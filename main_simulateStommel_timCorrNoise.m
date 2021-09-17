clear; clc;
%% Discretization
t0=0;
T=45;
dt=10^(-3);
% % for Rosenblatt process due to memory limitations
% T=10;
% dt=10^(-2);
dtval = strcat('dt',num2str(dt));
dtval = strrep(dtval,'.','K');
tspan = t0:dt:T;
M=10^4;
% % for Rosenblatt process
% M=10^3;

%% Choose type of noise
noise = 'BM';
%% Define parameters according to chosen distribution and generate sample paths of noise
switch noise
    case 'BM'
        %% Generate M Brownian motion (BM) paths
        rng(1)
        dW = sqrt(dt)*randn(M,length(tspan)-1);
%         W =cumsum([zeros(M,1) dW],2);
        dNoise = dW;
        noiseChoice = 'BM';
    case 'fBM'
        %% simulate increments of fraction Brownian motion (fBM) with Hurst parameter H
        %% inspired by http://www0.cs.ucl.ac.uk/staff/J.Shawe-Taylor/courses/ATML-1.pdf (last checked: 24.03.2020)
        H=0.9;
        C=covMatrixFBM(tspan,H);
        R = chol(C);
        rng(1)
        BH = [zeros(M,1) randn(M,length(tspan)-1)*R];
        dBH = diff(BH,1,2);
        dNoise = dBH;
        % nomenclature for saving figures
        Hval = strcat('H',num2str(H));
        Hval = strrep(Hval,'.','K');
        noiseChoice = strcat('fBM_',Hval);
    case 'color'
        %% Generate M Brownian motion (BM) paths
        rng(1)
        dW = sqrt(dt)*randn(M,length(tspan)-1);
%         W =cumsum([zeros(M,1) dW],2);
        y0 = 0;
        tau = 0.05;
        driftCoeff = @(x,t) -1/tau*x;
        diffusionCoeff = @(x,t) ones(size(x));
        OUP = zeros(M,length(tspan));
        for i=1:M
            OUP(i,:) = solveSDE_Euler(1, driftCoeff, diffusionCoeff, dt, tspan, dW(i,:), y0);
        end
        dOUP = diff(OUP,1,2);
        dNoise = dOUP;
        % nomenclature for saving figures
        tauval = strcat('tau',num2str(tau));
        tauval = strrep(tauval,'.','K');
        noiseChoice = strcat('color_',tauval);
    case 'Rosenblatt'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sparse implementation with cell arrays with only 1 dim and sparse matrices
    % concatenation immediately and no further cell array dimension
    H=0.9
    
    rng('default')
    dW = sqrt(dt)*randn(M,length(tspan)-1);
    % W =cumsum([zeros(M,1) dW],2);
    d_H = 1/(H+1)*(H/(2*(2*H-1)))^(-0.5);
    Hprime = (H+1)/2;
    % constant from Tudor(2008): "Analysis of the Rosenblatt process"
    c = sqrt((Hprime*(2*Hprime-1))/beta(2-2*Hprime,Hprime-0.5));
    
    dK = @(t,s) c*s^(0.5-Hprime)*(t-s).^(Hprime-1.5).*t.^(Hprime-0.5);
    y0 = zeros(M,1);
    RoPro = zeros(M,length(tspan));
    RoPro(:,1) = y0;
    R2 = sparse(length(tspan)-2,length(tspan)-2);
    C1 = cell(1,length(tspan)-2);
    for k=2:length(tspan)-1
        for l=2:length(tspan)-1
            v  = sparse(length(tspan)-2,1);
            if l==k
            C1{k-1} = [C1{k-1} v];
                continue;
            end
            v(max(k+1,l+1)-2:end) = cumsum(dK(tspan(max(k+1,l+1):end),tspan(k)).*dK(tspan(max(k+1,l+1):end),tspan(l))*dt)';
            C1{k-1} = [C1{k-1} v];
        end
    end
    for i=1:M
        for k=2:length(tspan)-1
            R2(:,k-1) = sum(C1{k-1}.*dW(i,1:end-1),2);
        end
        RoPro(i,3:end) = sum(R2.*dW(i,1:end-1),2)';
    end
    RoPro = d_H*RoPro;
    
    dRoPro = diff(RoPro,1,2);
    
    dNoise = dRoPro;
    % nomenclature for saving figures
    Hval = strcat('H',num2str(H));
    Hval = strrep(Hval,'.','K');
    noiseChoice = strcat('RoPro_',Hval);
end

%% Define ODE for parameter p as \dot(p) = \epsilon
p0 = 1.4;
% % for Rosenblatt process
% p0 = 1.0642;
epsilon = -0.01;

%% Define SDE
etaSquare = 7.5;
% determine initial value for SDE according to p0 such that it lies on critical manifold
y0 = fzero(@(x) -etaSquare*x.^3+15*x.^2-8.5*x+p0,1);
y0rep = y0*ones(M,1);
sigma = 0.01;
driftCoeff = @(p,x,t) p-x.*(1+etaSquare*(1-x).^2); % reduced Stommel Cessi model equation (64) from Kuehn (2013): "A mathematical framework for critical transitions: normal forms, variance, and applications"
diffusionCoeff = @(x,t) sigma*ones(size(x));

%% Solve Stommel Cessi SDE for chosen noise
[p, Y] = solveSDE1dim_EulerParameterized(driftCoeff, diffusionCoeff, dt, tspan, dNoise, y0rep,p0,epsilon);

%% Sample Var
varY_MC = var(Y);

%% for saving to .mat-file
Tval = strcat('T',num2str(T));
Tval = strrep(Tval,'.','K');
y0val = strcat('yIni',num2str(y0));
y0val = strrep(y0val,'.','K');
p0val = strcat('pIni',num2str(p0));
p0val = strrep(p0val,'.','K');
sigmaval = strcat('sigma',num2str(sigma));
sigmaval = strrep(sigmaval,'.','K');
epsval = strcat('eps',num2str(epsilon));
epsval = strrep(epsval,'.','K');

%% plot results
% sample variance evolution of Stommel Cessi SDE driven by chosen noise
fig1 = figure(1);
hold on;
from = 501; % corresponds to t=0.5
to = 2*10^3;
% from = 1001;
% % for Rosenblatt process
% from = 101;
% to = 0;
l3 = plot(p(from:end-to),varY_MC(from:end-to),'LineWidth',2);
ym = min(ylim);
yM = max(ylim);
pcrit1 = 11/9-1/sqrt(15); % from Kuehn (2013): "A mathematical framework for critical transitions: normal forms, variance, and applications"
xcrit1 = 1/15*(10+sqrt(15)); % from Kuehn (2013): "A mathematical framework for critical transitions: normal forms, variance, and applications"
l2 = plot([pcrit1 pcrit1],[ym yM],'k','LineWidth',2);
xlabel('Bifurcation parameter');
ylabel('Var (sample)');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend(l3,'Sample variance');
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
%% Save figure
% savefig(fig1,strcat('var_MC_Stommel_',noiseChoice,'_',Tval,y0val,p0val,sigmaval,epsval,dtval,'.fig'));
% saveas(fig1,strcat('var_MC_Stommel_',noiseChoice,'_',Tval,y0val,p0val,sigmaval,epsval,dtval,'.eps'),'epsc');

%% log-log-plot to discover inverse root scaling law
% sample variance evolution of Stommel SDE driven by standard BM
from = 501; % corresponds to t=0.5
to = 3000;
% % for Rosenblatt process
% from = 101;
% to = 150;
reg_MC = polyfit(log(p(from:end-to)-pcrit1),log(varY_MC(from:end-to)),1);
fig2 = figure(2);
hold on
l1 = plot(log(p(from:end-to)-pcrit1),log(varY_MC(from:end-to)),'LineWidth',2);
plot(log(p(from:end-to)-pcrit1),reg_MC(1)*log(p(from:end-to)-pcrit1)+reg_MC(2),'LineWidth',2);
xlabel('Logarithmic distance to bifurcation parameter');
ylabel('Logarithmic Var');
ax.FontSize = 12;
ax.Interpreter = 'latex';
str = ['slope= ',num2str(reg_MC(1))];
xm = min(xlim);
xM = max(xlim);
ym = min(ylim);
yM = max(ylim);
%% Comparable y-axis ranges
specRange = 2.5;
diff = specRange - (yM - ym);
ymNew = min(ylim)-diff/2;
yMNew = max(ylim)+diff/2;
ylim([ymNew yMNew]);
text(xm + (xM - xm)*.6, ymNew + (yMNew - ymNew)*.95, str,'FontSize',14)
grid on
set(gca,'FontSize',12);
%% Save figure
% savefig(fig2,strcat('loglog_var_MC_Stommel_',noiseChoice,'_',Tval,y0val,p0val,sigmaval,epsval,dtval,'.fig'));
% saveas(fig2,strcat('loglog_var_MC_Stommel_',noiseChoice,'_',Tval,y0val,p0val,sigmaval,epsval,dtval,'.eps'),'epsc');

%% Plot attracting critical manifold together with sample paths
fig3 = figure(3);
pcrit2 = 11/9+1/sqrt(15); % from Kuehn (2013): "A mathematical framework for critical transitions: normal forms, variance, and applications"
xcrit2 = 1/15*(10-sqrt(15)); % from Kuehn (2013): "A mathematical framework for critical transitions: normal forms, variance, and applications"
x_att1 = xcrit1:0.001:1.4;
x_rep = xcrit2:0.001:xcrit1;
x_att2 = 0:0.001:xcrit2;
h0 = @(x) x.*(1+etaSquare*(1-x).^2);
y_att1=h0(x_att1);
y_rep = h0(x_rep);
y_att2 = h0(x_att2);
plot(y_att1,x_att1,'k','LineWidth',2);
hold on
plot(y_rep,x_rep,'k--','LineWidth',2);
plot(y_att2,x_att2,'k','LineWidth',2);
plot(p,Y(1:min(M,100),:),'r');
xlabel('Bifurcation parameter');
ylabel('Salinity difference x(t)');
ax.FontSize = 12;
ax.Interpreter = 'latex';
plot(pcrit1,xcrit1,'o','Markersize',5,'MarkerFaceColor','b','MarkerEdgeColor','b');
set(gca,'FontSize',12);
% zoomed in commands based on
% https://de.mathworks.com/matlabcentral/answers/33779-zooming-a-portion-of-figure-in-a-figure,
% last checked: 11.02.21
% create a new pair of axes inside current figure
axes('position',[.55 .45 .25 .35])
box on % put box around new pair of axes
% for BM, fBM, colored noise
plot(y_att1(1:76),x_att1(1:76),'k','LineWidth',2);
hold on
plot(y_rep(440:end),x_rep(440:end),'k--','LineWidth',2);
plot(p(40000:end),Y(1:100,40000:end),'r');

% % for Rosenblatt process
% plot(y_att1(1:120),x_att1(1:120),'k','LineWidth',2);
% hold on
% plot(y_rep(380:end),x_rep(380:end),'k--','LineWidth',2);
% plot(p(1:end),Y(1:min(M,100),1:end),'r');
plot(pcrit1,xcrit1,'o','Markersize',5,'MarkerFaceColor','b','MarkerEdgeColor','b');
axis tight
%% Save figure
% savefig(fig3,strcat('pathsCMa_Stommel_',noiseChoice,'_',Tval,y0val,p0val,sigmaval,epsval,dtval,'.fig'));
% saveas(fig3,strcat('pathsCMa_Stommel_',noiseChoice,'_',Tval,y0val,p0val,sigmaval,epsval,dtval,'.eps'),'epsc');
