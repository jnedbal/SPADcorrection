function correctCDMbyIRF(saveOld)
% CORRECTCDMBYIRF function used the IRF calibration data to apply minor
% corrections to the CDM linearization data. Essentially, it corrects for
% the fact, that calibration of the bin size, achieved by providing a
% precise known clock during the CDM measurement is affected by the SPAD
% jitter. This routine takes the full-width hal-maximum of the IRF into
% account to correct the bin size.
%
% Syntax: analyzeCDM(saveOld)
%
% Inputs:
%   saveOld input is optional. If set to 1, it tells the function to keep
%   the non-corrected calibration in the 'correction' struct, in a field
%   called 'noIRF'. Generally, it is not recommended to use the input, as
%   it only takes up memory space. It is only for debugging purpose.
%
%   It requires a global variable 'correction', which carries the 
%   calibration data. The function will be modified by this function.
%
% Outputs:
%   The global variable struct 'correction' is updated with new values in
%   several of its fields:
%   
%   correction.avgBinWidth          The actual average bin width calculated
%                                     for each pixel based on the CDM
%                                     measurement.
%   correction.photonsPerPicosecond The average number of photons per 
%                                     picosecond for all pixels stored in a
%                                     2D array.
%   correction.binWidth             The actual calibrated bin width for
%                                     all pixels stored in a 3D array.
%   correction.idealBinWidth        The calibrated ideal bin width for
%                                     all pixels stored in a 3D array. It
%                                     is a 3D representation of avgBinWidth
%   correction.INL                  The integral nonlinearity for all
%                                     calibrated bins and pixels in a 3D
%                                     array.
%   correction.DNL                  The differential nonlinearity for
%                                     all calibrated bins and pixels in a
%                                     3D array.
%
%   correction.noIRF                Optional new field that holds the copy
%                                     of the fields before the IRF
%                                     correction. This is for debugging
%                                     only, when saveOld == 1;
%
% Examples:
%   correctCDMbyIRF
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: fitIRF, SPADcorrection

% Jakub Nedbal
% King's College London
% Aug 2020
% Last Revision: 11-May-2020 - Created this backup file
%
% Copyright 2020 Jakub Nedbal
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

% Create a waitbar
cn = 0;
hwaitbar = waitbar(cn, '', 'Name', 'Correcting CDM data ...');

% We can store the original information for future record, but it does take
% up memory and time, so it is conditional.
if nargin > 1 && saveOld
    % Store the original non-corrected value of average bin width
    correction.noIRF.avgBinWidth = ...
        correction.avgBinWidth;
    % Store the original non-corrected value
    correction.noIRF.photonsPerPicosecond = ...
        correction.photonsPerPicosecond;
    % Store the original non-corrected value
    correction.noIRF.binWidth = ...
        correction.binWidth;
    % Store the original non-corrected value
    correction.noIRF.idealBinWidth = ...
        correction.idealBinWidth;
    % Store the original non-corrected value
    correction.noIRF.idealBinWidth = ...
        correction.idealBinWidth;
    % Store the original non-corrected value
    correction.noIRF.INL = ...
        correction.INL;
    % Store the original non-corrected value
    correction.noIRF.DNL = ...
        correction.DNL;
end

% Calculate the average bin width in picoseconds [ps]
% Start by counting the bins in the rep rate period, corrected by the IRF
% full-width half-maximum
periodWidth = correction.lastBin - ...
              correction.firstBin - ...
              correction.IRF.peak.FWHMinterp;
% It could make sense to interpolate the curve around the first and last
% bin to get a sub-pixel resolution, but because of the unknown actual bin
% size, it is not possible to do accurately anyway. Therefore, I stick to
% the value rounded to integer bin counts. The resultant error should be
% less than 0.1% for a repetition period with 1000 bins. Therefore pretty
% negligible.
%

% Calculate the corrected average bin width taking the FWHM of the IRF into
% account.
correction.avgBinWidth = correction.repPeriod ./ periodWidth;

%% Calculate the average number of photons per bin for each pixel
%tmpCDM = correction.CDM;
%tmpCDM(zeroMap) = 0;
%correction.avgPhotons = sum(tmpCDM, 3) ./ ...
%                        (correction.lastBin - correction.firstBin + 1);

% Calculate the average number of photons per picosecond
% This number gives us the number of photons accumulated in each picosecond
% Calculate the corrected photons per picoseond.
correction.photonsPerPicosecond = correction.avgPhotons ./ ...
                                  correction.avgBinWidth;

% Calculate the actual bin width for each pixel
% Create a matrix size of CDM that stacks the photons per picosecond
photonsPerPicosecond = ...
    repmat(correction.photonsPerPicosecond, [1 1 size(correction.CDM, 3)]);
% Make the photons per picosecond infinite where its outside of the
% calibrated range. This makes the bin width 0 for those uncalibrated bins.
photonsPerPicosecond(~correction.calibratedBins) = Inf;
% Calculate the bin width in.
correction.binWidth = correction.CDM ./ photonsPerPicosecond;

% Create idealized bin width array
correction.idealBinWidth = ...
    repmat(correction.avgBinWidth, [1, 1, size(correction.CDM, 3)]);
% Make the bin width 0 for bins outside of calibrated range.
correction.idealBinWidth(~correction.calibratedBins) = 0;

% % Get the map of calibrated bins
% linearize.calibratedBins = LUT >= linearize.firstBin + ignoreBins & ...
%                            LUT <= linearize.lastBin - ignoreBins;

% Calculate integral non-linearity
% This is the accumulated difference of the actual and idealized TDC in 
% picoseconds
correction.INL = cumsum(correction.binWidth, 3) - ...
                     cumsum(correction.idealBinWidth, 3);
correction.INL(~correction.calibratedBins) = NaN;

% Calculate differential non-linearity
% This is the difference of the actual and idealized TDC in picoseconds
correction.DNL = correction.binWidth - ...
                     correction.idealBinWidth;
correction.DNL(~correction.calibratedBins) = NaN;


% % Update the waitbar the last time
% if ishandle(hwaitbar)
%     waitbar(1, hwaitbar, 'Saving the linearized data ...');
% end
% 
% % save the calibration data in a file
% save(linearize.files.binCorrection, 'linearize', '-v7.3')

% Close the waitbar
delete(hwaitbar)


