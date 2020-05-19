function fitIRF(XYZimage)
% fitIRF fits an exponentially modified Gaussian to the IRF
% A measurement of the instrument response function is done using the SPAD
% camera. Use a pulsed laser, quenched fluorescein, and accumulate many
% frames into a single histogram map.
%
% Syntax: fitIRF(XYZimage)
%
% Inputs:
%   The function requires a 3D array, where 1st and 2nd dimensions
%   represent th pixels od the SPAD array and the 3rd dimension is the time
%   domain with the histogram bins of the TDC.
%   The function also requires the global variable struct 'correction' with
%   the fields holding the parameters of the correction.
%
% Outputs:
%   The global variable struct 'correction' is populated with new fields:
%
%   correction.fit          Fit parameters of the exponentially modified
%                             Gaussian in each pixel
%   correction.fit.h          Exp-mod Gaussian amplitude
%   correction.fit.mu         Exp-mod Gaussian peak position
%   correction.fit.sigma      Exp-mod Gaussian standard deviation
%   correction.fit.tau        Exp-mod Gaussian exponential decay constant
%   correction.fit.offset     Exp-mod Gaussian exponential zero offset
%   correction.exitFlag       1 for successsful fit, 0 for failed fit
%   correction.goodfit        Matrix of good fits. This is generated based
%                               on the similarity of the results, not
%                               meeting fit stopping conditions. The value
%                               is 1, i.e. good fit, if sigma is between
%                               0.5x and 1.5x median sigma, tau is between
%                               0.5x and 1.5x median tau, and offset is not
%                               an  outlier according to 'isoutlier'.
%   correction.peak         The peak, maximum, rising- and fallin-edge
%                             positions, and full-width half-maximum of the
%                             fitted data. These values are calculated by
%                             iterrative fitting, as I could not find the
%                             analytical solution to the problem. It needs
%                             looking at in the future, as the solution
%                             must exist.
%   correction.peak.Pos       Exp-mod Gaussian peak position
%   correction.peak.Max       Exp-mod Gaussian peak maximum
%   correction.peak.risingHM  Exp-mod Gaussian peak rising edge position
%   correction.peak.fallingHM Exp-mod Gaussian peak falling edge position
%   correction.peak.FWHM      Exp-mod Gaussian peak full-width half maximum
%   correction.peak.XXinterp  Below are further fields that interpolate the
%                               above values, where the goodfit is not '1'
%   correction.peak.FWHMinterp
%   correction.peak.PosInterp
%   correction.peak.risingHMinterp
%   correction.peak.fallingHMinterp
%
% Examples:
%   fitIRF(XYZimage)
%
% Toolbox requirement: Optimization Toolbox
% Other m-files required: extractHistogramData, exgfit, exGauss
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: analyzeCDMmap, analyzeCDMhist, loadCDMdata

% Jakub Nedbal
% King's College London
% Aug 2018
% Last Revision: 08-May-2020 - Moved data import code into 'loadCDMdata'
% Revision: 02-Feb-2020 - Added support for repeated CDM measurements.
% Revision: 13-Aug-2018
%
% Copyright 2018-2020 Jakub Nedbal
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


%% Define global variables
global correction
global pixHist
global pixIndex


%% Create a waitbar
hwaitbar = waitbar(0, '', ...
                   'Name', 'Fitting IRF with exGaussian function...');
% Make it a big higher
hwaitbar.Position(4) = hwaitbar.Position(4) + 10;


%% Run a routine that extract all important aspects of the linearization
%  and the input histogram data
extractHistogramData;


%% Find the position of the peak in each pixel
fit.h = zeros(inputSize([1, 2]));
fit.mu = fit.h;
fit.sigma = fit.h;
fit.tau = fit.h;
fit.offset = fit.h;
fit.exitFlag = fit.h;
peak.Pos = nan(inputSize([1,2]));
peak.Max = peak.Pos;
peak.risingHM = peak.Pos;
peak.fallingHM = peak.Pos;

tic
numberPixels = double(numberPixels); %#ok<NODEF>
                                     % loaded by extractHistogramData 

%% Go through the image pixel by pixel
for i = 1 : numberPixels
    % update the waitbar
    if mod(i, 100) == 0
        if ishandle(hwaitbar)
            % Estimated time to completion
            eta = numberPixels / i * toc / 86400;
            if isinf(eta); eta = 0; end
            % Update the waitbar
            waitbar(i / numberPixels, ...
                    hwaitbar, ...
                    sprintf('%d of %d\n%s (%s ETA)', ...
                            i, ...
                            numberPixels, ...
                            datestr(toc / 86400, 'HH:MM:SS'), ...
                            datestr(eta, 'HH:MM:SS')));
        end
    end
    % Extract the histogram for pixel
    pixIndex = double(1 : (lastCalBin(i) - firstCalBin(i) + 1))';
    %pixHist = double(inputHist(firstCalBin(i) : lastCalBin(i), i));
    pixHist = double(XYZimage(pixIndex, i));
    % Find the maximum in each pixel
    [h0, mu0] = max(pixHist);
    sigma0 = 3;
    tau0 = 8;
    offset0 = median(pixHist);
    param0 = [h0 * tau0 * exp(1), mu0, sigma0, tau0, offset0];
    % Run the fit
    [fit.h(i), ...
        fit.mu(i), ...
        fit.sigma(i), ...
        fit.tau(i), ...
        fit.offset(i), ...
        fit.exitFlag(i)] = exgfit(param0);
    %fprintf('%d\n', i);
end

% Reshape outputs into matrices of the same size as the image sensor
%h = reshape(h, inputSize([1, 2]));
%mu = reshape(mu, inputSize([1, 2]));
%sigma = reshape(sigma, inputSize([1, 2]));
%tau = reshape(tau, inputSize([1, 2]));
%exitFlag = reshape(exitFlag, inputSize([1, 2]));

%% Calculate the well fitted bins
% These are ones that give results similar to the median
fit.goodfit = (fit.sigma < 1.5 * median(fit.sigma(:))) & ...
              (fit.sigma > 0.5 * median(fit.sigma(:))) & ...
              (~isoutlier(fit.offset)) & ...
              (fit.tau < 1.5 * median(fit.tau(:))) & ...
              (fit.tau > 0.5 * median(fit.tau(:)));
fit.goodfit([1, end], :) = false;

%% Calculate the fit parameters
tic
% Set options for the minimization function
options = optimoptions('fmincon', 'Display', 'none');
for i = find(fit.goodfit)'
    if mod(i, 100) == 0
        if ishandle(hwaitbar)
            % Estimated time to completion
            eta = numberPixels / i * toc / 86400;
            if isinf(eta); eta = 0; end
            % Update the waitbar
            waitbar(i / numberPixels, ...
                    hwaitbar, ...
                    sprintf('%d of %d\n%s (%s ETA)', ...
                            i, ...
                            numberPixels, ...
                            datestr(toc / 86400, 'HH:MM:SS'), ...
                            datestr(eta, 'HH:MM:SS')));
        end
    end
    % Find the position of the peak, and the the maximum of the peak
    % Minimization function
    mFun = @(x) (-(exGauss(x, ...
                           fit.h(i), ...
                           fit.mu(i), ...
                           fit.sigma(i), ...
                           fit.tau(i), ...
                           fit.offset(i))));
    % Run minimization routine to find the peak position
    peak.Pos(i) = fminsearch(mFun, fit.mu(i));
    % Evaluate the exGauss function to get the value of the peak
    peak.Max(i) = exGauss(peak.Pos(i), ...
                          fit.h(i), ...
                          fit.mu(i), ...
                          fit.sigma(i), ...
                          fit.tau(i), ...
                          fit.offset(i));
    % Create a minimization function for half-maximum
    mFun = @(x) abs(peak.Max(i) / 2 - ...
                    exGauss(x, fit.h(i), fit.mu(i), fit.sigma(i), ...
                               fit.tau(i), fit.offset(i)) + ...
                    fit.offset(i) / 2);
    % Run minimization routine to find the left half-maximum
    peak.risingHM(i) = fmincon(mFun, ...
                               peak.Pos(i) - fit.sigma(i), ...
                               [], [], [], [], ...
                               -Inf, ...
                               peak.Pos(i), ...
                               [], options);
    % Run minimization routine to find the right half-maximum
    peak.fallingHM(i) = fmincon(mFun, ...
                                peak.Pos(i) + fit.tau(i), ...
                                [], [], [], [], ...
                                peak.Pos(i), ...
                                Inf, ...
                                [], options);
end
% Calculate the full width at half-maximum
peak.FWHM = peak.fallingHM - peak.risingHM;

%% Run interpolations for the data that was not properly fitted
% interpolate the FWHM for points that are not fitted well
% Create a matrix of points
[X, Y] = meshgrid(1 : size(peak.FWHM, 2), 1 : size(peak.FWHM, 1));
Xv = X(~fit.goodfit);
Yv = Y(~fit.goodfit);
X = X(fit.goodfit);
Y = Y(fit.goodfit);

% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the full-width half maximum for the missing pixels
F = scatteredInterpolant(X, Y, peak.FWHM(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.FWHMinterp = peak.FWHM;
peak.FWHMinterp(~fit.goodfit) = F(Xv, Yv);



% Interpolate the peak position for points that are not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peak.Pos(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.PosInterp = peak.Pos;
peak.PosInterp(~fit.goodfit) = F(Xv, Yv);



% Interpolate the rising edge half-maximum for pixels not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peak.risingHM(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.risingHMinterp = peak.risingHM;
peak.risingHMinterp(~fit.goodfit) = F(Xv, Yv);



% Interpolate the falling edge half-maximum for pixels not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peak.fallingHM(fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peak.fallingHMinterp = peak.fallingHM;
peak.fallingHMinterp(~fit.goodfit) = F(Xv, Yv);

% Close the waitbar
delete(hwaitbar)

% Assign results to the main structure for future use
correction.IRF.fit = fit;
correction.IRF.peak = peak;




%% Run a routine to correct the IRF offset
correction.IRF.corrected = resampleHistogramPar(correction.IRF.raw);