function [ Y ] = solveSDE_Euler(dim, driftCoeff, diffusionCoeff, dt, tspan, dW, Y0)
%calculates the stochastic Euler approximation based on
    % dim = dimension of problem (scalar)
    % drift = driftcoefficient as a function of (x,t) (pay attention to order)
    % diffusion = diffusion coefficient as a function of(x,t) (pay attention to order)
    % dt = time increment (scalar), assumption: time increment constant!
    % tspan = vector containing points in time
    % dW = the Brownian increments corresponding to the length of tspan
    % Y0 = starting value
% and gives back
    %the approximation as a dim x length(tspan) matrix with one discretisation time corresponding to
    %one column, that is 'time runs over columns'

    Y = zeros(dim,length(tspan));
    %Fill in initial value
    Y(:,1) = Y0;
    %solution for step size dt according to stochastic Euler
    %recursion Y_{n+1} = Y_n + a(t_n,Y_n)*dt_n + b(t_n,Y_n)*dW_n
    for i=1:length(tspan)-1
        Y(:,i+1) = Y(:,i) + driftCoeff(Y(:,i),tspan(i))*dt + diffusionCoeff(Y(:,i),tspan(i)).*dW(:,i);
    end
end