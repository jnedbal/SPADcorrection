% SPADCORRECTION is a pipeline which calls scripts to analyse SPADs and
% create a dataset for their correction.
%
% It requires a dataset of code density maps (CDM) and a dataset with
% instrument response function (IRF).
%
% Syntax: SPADcorrection
%
% Inputs:
%   The function does not need any input. It will ask questions at the
%   start about location of files and measurement parameters.
%
% Outputs:
%   MAT-file containing struct 'correction' is produced. There is a second
%   mat file produced, which contains the minimum amount of data to perform
%   the correction. The complete dataset is stored in a file with the
%   extension ".FULL.MAT", whereas the minimal version just with the
%   extension ".MAT".
%   The struct 'correction' contains the following fields:
%
%   correction.files                Files of the source datasets and 
%                                     location of the stored data
%   correction.repRate              Experimental laser repetition rate
%                                     expressed in MHz
%   correction.CDM                  Input raw density map with the first
%                                     bin containing many zeros removed
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
%   correction.avgPhotons           Average number of photons in a time bin
%                                     for all pixels stored in a 2D array
%   correction.photonsPerPicosecond Average photon count per picosecond for
%                                     all pixels stored in a 2D array
%   correction.binWidth             Actual calibrated bin width for all 
%                                     pixels stored in a 3D array
%   correction.idealBinWidth        Calibrated ideal bin width for all
%                                     pixels stored in a 3D array
%   correction.INL                  3D array of integral nonlinearity for
%                                     all calibrated bins and pixels
%   correction.DNL                  3D array of differential nonlinearity
%                                     for all calibrated bins and pixels
%   correction.IRF                  Parameters of the IRF calculated from
%                                     fitting  of exponentially modified
%                                     Gaussian function to the raw data
%   correction.IRF.raw                Raw IRF measurement data as loaded
%                                       from the experimental data file
%   correction.IRF.linear             Linearized IRF measurement based on
%                                       the CDM linearization routine
%   correction.noIRF                (Optional) This field is produced if
%                                     needed during the IRF correction
%                                     stage to store the non-corrected data
%
%   The fields in the minial version of the corretion struct are:
%   'calibratedBins', 'binWidth', 'avgBinWidth', and 'IRF.peak.PosInterp'
%
%   Four figures are produced, if required:
%     * Interactive figure that graphically shows the performance of the
%       sensor, the TDC histograms, the INL, DNL, and Fourier transform.
%     * Histogram of TDC bin widths
%     * Histogram of number of bins in a laser period
%     * Example histogram of bin width calibration
%
% Examples:
%   SPADcorrection
%
% Toolbox requirement: Optimization Toolbox, Parallel Processing Toolbox
%                      (optional)
% Other m-files required: loadCDMdata, analyzeCDM, fitIRF, correctIRFbyCDM,
%                         analyzeCDMmap, analyzeCDMhist, exGauss, exgfit,
%                         resampleHistogramPar, syntheticPhotons (MEX-file)
% Subfunctions: none
%
% See also: analyzeCDM, fitIRF

% Jakub Nedbal
% King's College London
% May 2020
% Last Revision: 11-May-2020 - First version of the file
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

%% Remove the global variable. If it stays in the memory, it will corrupt
%  the IRF correction routine in this script
clear global correction

%% Define global variables
global correction


%% Get the inputs from the user
% Start by asking a few questions at the beginning
% This way, the code can be left running unattended

% Open a dialog box to select files with code density map data
[file, path] = ...
    uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
               '*.*', 'All Files (*.*)'}, ...
              'Select files with code-density map measurements', ...
              'MultiSelect', 'on');
% Check if nothing has been returned
if isequal(file, 0) || isequal(path, 0)
    return
end
CDMfiles = strcat(path, file);


% Locate the IRF file
[file, path] = uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                          '*.*', 'All Files (*.*)'}, ...
                         'Select IRF file');

% Keep the file empty if there is no IRF to process
if isequal(file, 0) || isequal(path, 0)
    correction.files.IRF = [];
else
    % Store the IRF file name for later
    correction.files.IRF = fullfile(path, file);
end


% Open a list selection box with graphical format
list = {'EPS', 'PDF', 'PNG'};
in = listdlg('ListString', list, ...
             'Name', 'Image Format', ...
             'ListSize', [300, 300], ...
             'PromptString', {'Select image format(s)', ...
                              'for graphical outputs:'});
% Don't do anything if empty result was returned
if isempty(in)
    correction.files.graphics = [];
else
    correction.files.graphics = list{in};
end


% Ask the user for the laser repetition rate
correction.repRate = inputdlg('Laser Repetition Rate [MHz]', ...
                             'Linearize SPADs', ...
                             1, ...
                             {'20'});
% Stop the function if Cancel pressed
if isempty(correction.repRate)
    return
end
% Convert the returned string of repetition rate to numbers
correction.repRate = str2double(correction.repRate{1});


% Ask the user for the place where to store the correction data
[file, path] = uiputfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                         '*.*', 'All Files (*.*)'}, ...
                         'Save the correction data', ...
                         'binCorrection.mat');
if isequal(file, 0) || isequal(path, 0)
    return
end
correction.files.binCorrection = fullfile(path, file);


%% Load the CDM data and store it into the correction global variable
loadCDMdata(CDMfiles);


%% Analyze the code density map data. Load all the selected files and 
%  analyze the widths of individual histogram bins in every pixel
analyzeCDM;


%% Run IRF correction, if the instrument function is loaded and linearized
if ~isempty(correction.files.IRF)
    % Make a comment
    fprintf('Loading IRF file %s...\n', correction.files.IRF);
    % Load the IRF file
    load(correction.files.IRF, 'XYZimage')

    % Linearize the IRF data
    correction.IRF.raw = XYZimage;
    % Make a comment
    fprintf(['Linearizing the IRF. This may take a while,\n', ...
             'depending on the number of photons.\n']);
    % Run the linearization routine
    correction.IRF.linear = resampleHistogramPar(XYZimage);

    % Fit the exponentially modified Gaussian models to the IRF data
    fitIRF(correction.IRF.linear);
    % Correct the bin size by the IRF
    correctCDMbyIRF;
    % Create parameters for FLIM fitting by SlimCurve
    fitFLIMparam;
end


%% Plot the results, unless there is no need to
if ~isempty(correction.files.graphics)
    analyzeCDMmap(correction.files.graphics, ...
                  fileparts(correction.files.binCorrection));
    analyzeCDMhist(correction.files.graphics, ...
                   fileparts(correction.files.binCorrection));
end



%% Save the calibration data in a file
% SAve two files. One adding full. This is a large file with all the data
% contained within. Then there is another file, which contains the minimum
% amount of data for faster loading and less memory demand.
[path, fname, ext] = fileparts(correction.files.binCorrection);
correction.files.binCorrectionLong = fullfile(path, [fname, '.full' ext]);
% Make a comment
fprintf(['Saving results into file %s.\n', ...
         'This can take a while. Stand by ...\n'], ...
        correction.files.binCorrectionLong);
% save the results
save(correction.files.binCorrectionLong, 'correction', '-v7.3')
% Get rid of the redunndant variables
clear path; clear fname; clear ext

% Create a temporary struct that will be used to store correction in a
% smaller file with just the minimum set of data required for correction
S.correction.calibratedBins = correction.calibratedBins;
S.correction.binWidth = correction.binWidth;
S.correction.avgBinWidth = correction.avgBinWidth;
S.correction.IRF.peak.PosInterp = correction.IRF.peak.PosInterp;
S.correction.IRF.fit.goodfit = correction.IRF.fit.goodfit;
S.correction.fitFLIM = correction.fitFLIM;
% Make a comment
fprintf(['Saving results into file %s.\n', ...
         'This can take a while. Stand by ...\n'], ...
        correction.files.binCorrection);
% save the results
save(correction.files.binCorrection, '-struct', 'S', '-v7.3')
clear S;