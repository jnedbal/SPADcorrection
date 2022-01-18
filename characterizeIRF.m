function characterizeIRF(XYZimage)
% fitIRF fits an Gaussian to the IRF to find its peak
% A measurement of the instrument response function is done using the SPAD
% camera. Use a pulsed laser, quenched fluorescein, and accumulate many
% frames into a single histogram map.
%
% Syntax: fitIRF(XYZimage)
%
% Inputs:
%   The function requires a 3D array, where 1st and 2nd dimensions
%   represent the pixels od the SPAD array and the 3rd dimension is the
%   time domain with the histogram bins of the TDC.
%   The function also requires the global variable struct 'correction' with
%   the fields holding the parameters of the correction.
%
% Outputs:
%   The global variable struct 'correction' is populated with new fields:
%
%   correction.fit          Fit parameters of the Gaussian in each pixel
%   correction.IRF.fit.               Struct of IRF fit output parameters
%       .param.filter.model           Model for smoothing IRF data
%       .param.filter.kernel          Kernel size for smoothing IRF data
%   correction.IRF.fit.h              Gaussian amplitude (2D array)
%   correction.IRF.fit.mu             Gaussian peak position (2D array)
%   correction.IRF.fit.sigma          Gaussian standard dev. (2D array)
%   correction.IRF.fit.offset         Gaussian offset (2D array)
%   correction.IRF.fit.rsquare        R-squared of the fit (2D array)
%   correction.IRF.fit.exitFlag       Converged fit flag
%                                       1 for successsful fit
%                                       0 for failed fit
%   correction.IRF.fit.goodfit        Matrix of good fits. This is 
%                                       generated based on the similarity
%                                       of the results, not meeting fit 
%                                       stopping conditions. The value is
%                                       1, i.e. good fit, if sigma is not
%                                       an outlier, and r-squared is more
%                                       than 0.95
%   correction.IRF.fit.interp.h       Gaussian amplitude interpolated
%   correction.IRF.fit.interp.mu      Gaussian peak position interpolated
%   correction.IRF.fit.interp.sigma   Gaussian standard dev. interpolated
%   correction.IRF.fit.interp.offset  Gaussian offset interpolated
%   correction.IRF.peak             The peak, maximum, rising- and falling-
%                                     edge positions, and full-width half-
%                                     maximum of the fitted data.
%   correction.IRF.peak.Pos           Gaussian peak position
%   correction.IRF.peak.Max           Gaussian peak maximum
%   correction.IRF.peak.FWHM          Gaussian peak full-width at half max
%
%
% Examples:
%   fitIRF(XYZimage)
%
% Toolbox requirement: none
% Other m-files required: extractHistogramData
% Subfunctions: fitGauss
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: analyzeCDMmap, analyzeCDMhist, loadCDMdata

% Jakub Nedbal
% King's College London
% Apr 2021
% Last Revision: 15-Apr-2021 - Created a new function to fit Gaussians
%                              rather than exponentially modified Gaussians
%                              Fixed the problem with inequal bin sizes
% 
%
% Copyright 2018-2021 Jakub Nedbal
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
global currPeak
global fitRange

%% Define analysis parameters
% see help smoothdata for options
Fit.param.filter.model = 'gaussian';
Fit.param.filter.kernel = 9;

%% Create a waitbar
hwaitbar = waitbar(0, '', 'Name', 'Finding IRF parameters...');
% Make it a big higher
hwaitbar.Position(4) = hwaitbar.Position(4) + 10;

%% Run a routine that extract all important aspects of the linearization
%  and the input histogram data
extractHistogramData;

%% Smoothen the data to find the peaks
IRFsmooth = smoothdata(double(XYZimage), 1, ...
                       Fit.param.filter.model, Fit.param.filter.kernel);

%% Find the position of the peak in each pixel
Fit.h = zeros(inputSize([1, 2]));
Fit.mu = Fit.h;
Fit.sigma = Fit.h;
Fit.offset = Fit.h;
%Fit.gof = struct('sse', Fit.h, ...
%                 'rsquare', Fit.h, ...
%                 'dfe', Fit.h, ...
%                 'adjrsquare', Fit.h, ...
%                 'rmse', Fit.h);
Fit.rsquare = Fit.h;
Fit.exitflag = false(inputSize([1, 2]));
Peak.Pos = NaN(inputSize([1, 2]));
Peak.Max = Peak.Pos;
%Peak.risingHM = Peak.Pos;
%Peak.fallingHM = Peak.Pos;

%% Find the peaks
[peak, index] = max(IRFsmooth);
% Setup a zoom of relevant fit data
fitRangeFull = (-49 : 49)';
fitRange = fitRangeFull;
% Extract the peaks
peaks = zeros(numel(fitRange), size(IRFsmooth, 2));
% Fixed fit parameters
B = 0;      % Peak Position
C = 6.42;   % Sigma = C/sqrt(2)

tic
numberPixels = double(numberPixels); %#ok<NODEF>
options = optimset('Display', 'off');
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
    fitRange = fitRangeFull;
    binIndex = index(i) + fitRange;
    % If the peak is too close to the start, negative indices of the IRF
    % will be called and result in an error, for this purpose, only
    % positive indices can be selected.
    validIndex = binIndex > 0;
    currPeak = IRFsmooth(binIndex(validIndex), i);
    fitRange = fitRange(validIndex);
    % Now do the extrapolation. If the current peak is in full, just leave
    % it as is

    
    %currPeak = IRFsmooth(binIndex, i);
    peaks(validIndex, i) = currPeak;
    %D = mean(currPeak([1 : 4, end - (0 : 3)])); % Offset
    D = mean(currPeak(end - (0 : 7)));          % Offset
    A = peak(i) - D;                            % Amplitude
    % Perform the Gaussian fit
    %[f1, gof] = fit(fitRange, currPeak, gaussEqn, 'Start', [A, B, C, D]);
    [fp1, res, exitflag] = fminsearch(@fitGauss, [A, B, C, D], options);
    % Deal the results
    Fit.h(i) = fp1(1);
    Fit.mu(i) = index(i) + fp1(2);
    Fit.sigma(i) = fp1(3) / sqrt(2);
    Fit.offset(i) = fp1(4);
    % Calculate R-squared
    Fit.rsquare(i) = 1 - res / sum((mean(currPeak) - currPeak) .^ 2);
    % Store the exit flag of the minimization
    Fit.exitflag(i) = exitflag == 1;
    % If the selection was trunkated because the IRF was too early in the
    % range, use to offset to simulate the missing values
    if any(~validIndex)
        % Fit the measured values cropped at the start
        peaks(~validIndex, i) = ...
            fp1(1) * exp(-(fitRangeFull(~validIndex) - fp1(2)) .^ 2 / ...
            (fp1(3) ^ 2)) + fp1(4);
    end
        
end

%% Find goodfit
% Assume that adjusted r^2 less than 0.95 is not good
Fit.goodfit = Fit.rsquare > 0.95;
% Only keep those fits that converged well
Fit.goodfit = Fit.goodfit & Fit.exitflag;
% Find the well fitted indices
gf = find(Fit.goodfit);
% Assume that any peak with sigma being an outlier is not good
Fit.goodfit(gf(isoutlier(Fit.sigma(Fit.goodfit)))) = false;

%% Interpolate values
% interpolate values in pixels that are not fitted well
% Create a matrix of points
[X, Y] = meshgrid(1 : inputSize(2), 1 : inputSize(1));
Xv = X(~Fit.goodfit);
Yv = Y(~Fit.goodfit);
X = X(Fit.goodfit);
Y = Y(Fit.goodfit);

% Data to interpolate
params = {'sigma', 'h', 'mu', 'offset'};
% Interpolate one dataset after the other

for i = 1 : numel(params)
    % Set the interpolated value
    Fit.interp.(params{i}) = Fit.(params{i});
    % Create an interpolant object that can be called
    % Use the natural neighbor interpolation method and 
    % nearest neightbor extrapolation method
    % to estimate the full-width half maximum for the missing pixels
    F = scatteredInterpolant(X, Y, Fit.(params{i})(Fit.goodfit), ...
                             'natural', 'nearest');

    % Create a new matrix, which contains the estimated and the interpolated
    % values of the full-width half-maximum
    Fit.interp.(params{i})(~Fit.goodfit) = F(Xv, Yv);
end

%% Calculate the peak parameters
Peak.FWHM = sqrt(8 * log(2)) * Fit.interp.sigma;
Peak.Pos = Fit.interp.mu;
Peak.Max = Fit.interp.h;
%Peak.risingHM = Peak.Pos - Peak.FWHM / 2;
%Peak.fallingHM = Peak.Pos + Peak.FWHM / 2;

% Close the waitbar
delete(hwaitbar)

%% Store the calculated figures in the correction struct

% Assign results to the main structure for future use
correction.IRF.fit = Fit;
correction.IRF.peak = Peak;
% Add the average bin width from the good pixels

% Calculate the global bin width used in the linearized data
correction.globalBinWidth.corr = ...
    mean(correction.avgBinWidth(correction.IRF.fit.goodfit));
end

function z = fitGauss(p)
global currPeak
global fitRange
%cx = p(1);
%wx = p(2);
%amp = p(3);

zx = p(1) * exp(-(fitRange - p(2)) .^ 2 / (p(3) ^ 2)) + p(4) - currPeak;

z = sum(zx.^2);
end
