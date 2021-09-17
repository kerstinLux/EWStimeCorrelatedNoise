function [C] = covMatrixFBM(tspan,H)
% calculates the covariance matrix of a fractional Brownian motion (fBM)
% for time>0 (time=dt und nicht time=0)
% with Hurst parameter H
C = zeros(length(tspan)-1,length(tspan)-1); % initialize covariance matrix
for i=1:length(tspan)-1
    for j=1:length(tspan)-1
        C(i,j) = 0.5*(tspan(i+1).^(2*H)+tspan(j+1).^(2*H) - abs(tspan(j+1)-tspan(i+1)).^(2*H));
    end
end
end

