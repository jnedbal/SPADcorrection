function analyzeCDM
% ANALYZECDM processes a code density map from a SPAD array
% Function that processes measured code density map (CDM) from a SPAD array
% sensors and infers the bin size for each time bin and pixel.
%
% Syntax: analyzeCDM
%
% Inputs:
%   The function does not need any input. It required a global variable
%   struct 'correction' with the field 'CDM'. This field must hold a 3D
%   code-density map of measurements from the SPAD camera.
%
% Outputs:
%   The global variable struct 'correction' is populated with new fields:
%
%   correction.medianPhotons        Median of the number of photons in 
%                                     non-zero bins for each pixel
%   correction.lastBin              2D array of indices of last valid 
%                                     bins for each pixel histogram
%   correction.firstBin             2D array of indices of first valid 
%                                     bins for each pixel histogram
%   correction.calibratedBins       3D array of booleans, where 1 stand for
%                                     the calibrated bins for each pixel
%   correction.repPeriod            Experimental laser repetition rate
%                                     expressed in picoseconds
%   correction.nrBins               Number of active bins in CDM
%                                     measurement (2D array)
%   correction.avgBinWidth          Average bin width of active pixels in
%                                     CDM measurement (2D array)
%   correction.globalBinWidth.raw   Linearized bin width in picoseconds
%                                     calculated from all pixels
%   correction.avgPhotons           Average number of photons in a time bin
%                                     for all pixels stored in a 2D array
%   correction.photonsPerPicosecond Average photon count per picosecond for
%                                     all pixels stored in a 2D array
%   correction.binWidth             Actual calibrated bin width for all 
%                                     pixels stored in a 3D array
%   correction.INL                  3D array of integral nonlinearity for
%                                     all calibrated bins and pixels
%   correction.DNL                  3D array of differential nonlinearity
%                                     for all calibrated bins and pixels
%
% Examples:
%   analyzeCDM
%
% Toolbox requirement: none
% Other m-files required: loadCDMdata
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: analyzeCDMmap, analyzeCDMhist, loadCDMdata

% Jakub Nedbal
% King's College London
% Aug 2018
% Last Revision: 15-Apr-2020 - Changed analysis to support global bin width
% Revision: 08-May-2020 - Moved data import code into 'loadCDMdata'
% Revision: 02-Feb-2020 - Added support for repeated CDM measurements.
% Revision: 13-Aug-2018
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


%% Create a waitbar
cn = 0;
hwaitbar = waitbar(cn, '', 'Name', 'Processing the data ...');


%% The one bin is particularly high, which due to zeros from, when the TDCs
% haven't triggered. Get rid of it.
% Find the highest bin
[~, in] = max(sum(correction.CDM, [1 2]));
% Make that bin zero
correction.CDM(:, :, in) = 0;

% Check if the input data has been flipped, if the zero bin is towards the
% beginning, flip the matric along the 3rd dimension
if (size(correction.CDM, 3) - in) > (in - 1)
    correction.CDM = flip(correction.CDM, 3);
end


%% Convert the CDM into doubles
correction.CDM = double(correction.CDM);


%% Find which bins contain nonzero photons, first average the pixels.
% Because there might be some noise, define nonzero as more than 0.1-times
% the maximum photons
% Calculate the average value for each bin across all pixels
avgCDM = mean(correction.CDM, [1, 2]);
% Find the bins that exceed 0.1-times the average value
nonzeroBins = avgCDM > max(avgCDM) * 0.1;
% For each pixel find the median of photons in the nonzero bins
correction.medianPhotons = median(correction.CDM(:, :, nonzeroBins), 3);


%% Find the bins, which have less than a half of the median photon count
% This routine can be improved to use interpolation and find the cut-off
% with a sub-bin precision. But this process is questionable, as each bin
% has different bin width. Therefore, It is all done rounded to single
% bins.
% First create the map of bin values less than half the median
zeroMap = correction.CDM < correction.medianPhotons / 2;
% Then create a lookup table with bin indices
% first get the size of the CDM
CDMsize = size(correction.CDM);
% Create a matrix the size of CDM, but filled with bin indices
LUT = repmat(reshape(1 : CDMsize(3), 1, 1, CDMsize(3)), ...
             [CDMsize(1), CDMsize(2), 1]);
% Keep only the nonzero bin values in the LUT
% Make the zero values the maximum for uint16
LUT(zeroMap) = NaN;
% Extract the first and last bins
correction.lastBin = max(LUT, [], 3);
correction.firstBin = min(LUT, [], 3);

% Create map of calibrated bins, which is ignoring the first and last 
% 20 bins in each histogram. These are always a bit rubbish.
ignoreBins = 20;
correction.calibratedBins = LUT >= correction.firstBin + ignoreBins & ...
                            LUT <= correction.lastBin - ignoreBins + 2;


%% Convert to repetition period in picoseconds [ps]
correction.repPeriod = 1e6 / correction.repRate;

%% Calculate the number of bins with CDM in each pixels
correction.nrBins = correction.lastBin - correction.firstBin + 1;

%% Calculate the average bin width in picoseconds [ps]
% ******  WARNING   ******
% This is probably not the right way to do it as there is jitter to
% consider. We would need to also subtract the FWHM of the instrument
% response function. But I don't have this measurement at the moment.
% *  Change: May 2020  *
% The code stays the same here, but there is another function 'fitIRF'
% called later that corrects for the SPAD jitter later. The results below
% will be eventually discarded and and replaced with the corrected results.
% However, for now they are required for the function fitIRF to work.
% ******  WARNING  ******
correction.avgBinWidth = correction.repPeriod ./ correction.nrBins;
% Calculate the global bin width used in the linearized data
correction.globalBinWidth.raw = mean(correction.avgBinWidth, 'all');

%% Calculate the average number of photons per bin for each pixel
tmpCDM = correction.CDM;
tmpCDM(zeroMap) = 0;
correction.avgPhotons = sum(tmpCDM, 3) ./ correction.nrBins;


%% Calculate the average number of photons per picosecond
% This number gives us the number of photons accumulated in each picosecond
%correction.photonsPerPicosecond = correction.avgPhotons ./ ...
%                                  correction.avgBinWidth;
correction.photonsPerPicosecond = sum(tmpCDM, 3) / correction.repPeriod;


%% Calculate the actual bin width for each pixel
% Create a matrix size of CDM that stacks the photons per picosecond
%photonsPerPicosecond = ...
%    repmat(correction.photonsPerPicosecond, [1 1 CDMsize(3)]);
% Make the photons per picosecond infinite where its outside of the
% calibrated range. This makes the bin width 0 for those uncalibrated bins.
%photonsPerPicosecond(~correction.calibratedBins) = Inf;
% Calculate the bin width in.
%correction.binWidth = correction.CDM ./ photonsPerPicosecond;
tmpCDM(~correction.calibratedBins) = NaN;
correction.binWidth = tmpCDM ./ correction.photonsPerPicosecond;


%% Create idealized bin width array
%correction.idealBinWidth = ...
%    repmat(correction.avgBinWidth, [1, 1, size(correction.CDM, 3)]);
% Make the bin width 0 for bins outside of calibrated range.
%correction.idealBinWidth(~correction.calibratedBins) = 0;


%% Calculate integral non-linearity
% This is the accumulated difference of the actual and idealized TDC in 
% picoseconds
%correction.INL = cumsum(correction.binWidth, 3) - ...
%                 cumsum(correction.idealBinWidth, 3);
correction.INL = cumsum(correction.binWidth - correction.avgBinWidth, ...
                        3, 'omitnan');
correction.INL(~correction.calibratedBins) = NaN;


%% Calculate differential non-linearity
% This is the difference of the actual and idealized TDC in picoseconds
correction.DNL = correction.binWidth - correction.avgBinWidth;
%correction.DNL(~correction.calibratedBins) = NaN;


%% Close the waitbar
delete(hwaitbar)


