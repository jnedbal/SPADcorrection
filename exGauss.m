% Define the ExGaussian function
function y = exGauss(x, h, mu, sigma, tau, offset)
    y = (h ./ tau) .* ...
        exp((mu ./ tau) + (sigma .^ 2 ./ (2 * tau .^ 2)) - (x ./ tau)) .* ...
        ncdf((x - mu - (sigma .^ 2 ./ tau)) ./ sigma) + ...
        offset;
end


% Define complementary error function
function y = ncdf(x)
    y = 0.5 * (1 + erf(x / sqrt(2)));
end

