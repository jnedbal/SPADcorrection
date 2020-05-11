function [h, mu, sigma, tau, offset, exitFlag] = exgfit(param0)
% exgfit  Fit ExGaussian distribution to data
%
% [MU,SIGMA,TAU] = exgfit(X,S) fits the ExGaussian distribution to data in
% vector X using maximum likelihood and returns the fitted parameters MU,
% SIGMA, and TAU. The ExGaussian distribution is formed by the convolution
% of independent normal and exponential observations. MU and SIGMA denotes
% the mean and standard deviation of the normal component and TAU denotes
% the mean of the exponential component. S is a three-element vector of
% starting values for MU, SIGMA, and TAU when fitting the distribution to
% data. SIGMA must be positive and TAU at least zero. 
%
%   Example
%         n = 200;
%         mu = 500; sigma = 200; tau = 400;
%         x = randn(1,n)*sigma+mu-log(1-rand(1,n))*tau;
%         s = [100,100,100];
%         [mu,sigma,tau] = exgfit(x,s); 
%         % Below plots the results without Statistics Toolbox
%         ncdf = @(x) 0.5*(1+erf(x/sqrt(2)));
%         epdf = @(x,mu,sigma,tau) (1/tau).*exp((mu/tau)+(sigma^2/(2*tau^2))-...
%             (x/tau)).*ncdf( (x-mu-(sigma^2/tau))./sigma);
%         histogram(x,'Normalization','pdf','Facecolor',[.6,.6,.6]);
%         hold on
%         xx = linspace(min(x),max(x),1e3);
%         plot(xx,epdf(xx,mu,sigma,tau),':k','linew',2);
%         set(gca,'Color',[.98,.98,.98],'FontSize',15)
%         grid on
%         xlabel('X','FontSize',25)
%         ylabel('pdf','FontSize',25)
%
%
%   Original code by Tobias Johansson
%   Tobias Johansson (2020). exgfit - Fit ExGaussian distribution to data,
%   (https://www.mathworks.com/matlabcentral/fileexchange/70225-exgfit-fit-exgaussian-distribution-to-data),
%   MATLAB Central File Exchange. Retrieved April 27, 2020. 

% Define complementary error function
% ncdf = @(x) 0.5 * (1 + erf(x / sqrt(2)));
% Define the ExGaussian function
%epdf = @(x, h, mu, sigma, tau) ...
%        (h / tau) .* ...
%        exp((mu / tau) + (sigma ^ 2 / (2 * tau ^ 2)) - (x / tau)) .* ...
%        ncdf((x - mu - (sigma ^ 2 / tau)) ./ sigma);
% Define the merit function
% L = @(p) sum(-log(epdf(x, p(1), p(2), p(3))));
%L = @(x, p) sum(abs(epdf((1 : numel(x))', p(1), p(2), p(3)) - x / sum(x)));
% Run the minimization
% pfit = fmincon(L, s, [], [], [], [], [-Inf, 0 + eps, 0], []);
% Deal the outcome
options = optimset('MaxFunEvals', 400 * numel(param0), ...
                   'Display', 'off');
[pfit, ~, exitFlag] = fminsearch(@costFun, param0, options);
[h, mu, sigma, tau, offset] = ...
    deal(pfit(1), pfit(2), pfit(3), pfit(4), pfit(5));
end

% Define complementary error function
function y = ncdf(x)
    y = 0.5 * (1 + erf(x / sqrt(2)));
end

% Define the ExGaussian function
function y = epdf(h, mu, sigma, tau, offset)
    global pixIndex
    y = (h / tau) .* ...
        exp((mu / tau) + (sigma ^ 2 / (2 * tau ^ 2)) - (pixIndex / tau)) .* ...
        ncdf((pixIndex - mu - (sigma ^ 2 / tau)) ./ sigma) + ...
        offset;
end

function cost = costFun(p)
    global pixHist
    %h = p(1)*p(4)*exp(1);%/ncdf(( - (p(3) ^ 2 / p(4))) ./ p(3))
    cost = sum(abs(epdf(p(1), p(2), p(3), p(4), p(5)) - pixHist));
end