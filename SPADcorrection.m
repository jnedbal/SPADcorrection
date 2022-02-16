% SPADCORRECTION is a pipeline which calls scripts to analyse SPADs and
% create a dataset for their correction.
%
% It requires a dataset of code density maps (CDM), a dataset with
% instrument response function (IRF), and a dataset with dark count map
% (DCR)
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
%   correction.lastBin              Indices of last valid bins for each
%                                     pixel histogram (2D array)
%   correction.firstBin             Indices of first valid bins for each
%                                     pixel histogram (2D array)
%   correction.calibratedBins       3D array of booleans where 1 stand for
%                                     the calibrated bins for each pixel
%   correction.repPeriod            Experimental laser repetition rate
%                                     expressed in picoseconds
%   correction.nrBins               Number of active bins in CDM
%                                     measurement (2D array)
%   correction.avgBinWidth          Average bin width of active pixels in
%                                     CDM measurement (2D array)
%   correction.globalBinWidth.raw   Linearized bin width in picoseconds
%                                     calculated from all pixels
%   correction.globalBinWidth.corr  Linearized bin width in picoseconds
%                                     calculated from good IRF fit pixels
%   correction.avgPhotons           Average number of photons in a time bin
%                                     for all pixels (2D array)
%   correction.photonsPerPicosecond Average photon count per picosecond for
%                                     all pixels (2D array)
%   correction.binWidth             Actual calibrated bin width for all 
%                                     pixels (3D array)
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
%   correction.IRF.fit.                Struct of IRF fit output parameters
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
%   correction.IRF.peak.corrPosInterp This is the peak position for the
%                                       resampled and timing-skew
%                                       corrected IRF
%   correction.IRF.corrected        Linearized and timing skew-corrected
%                                     IRF
%   correction.fitFLIM              A struct of parameters for
%                                     Levenberg-Marquardt fluorescence
%                                     lifetime decay fitting using SLIM
%                                     Curve
%   correction.fitFLIM.peak           Average bin index of IRF peak
%   correction.fitFLIM.start          The starting bin of decay data, the
%                                       IRF rising edge
%   correction.fitFLIM.fit_start      The starting of the fit, the IRF 
%                                       falling edge
%   correction.fitFLIM.fit_end        The last bin index with useful data
%                                       in all SPAD array pixels
%   correction.fitFLIM.data_start     The fitst bin with useful data in all
%                                       SPAD array pixels
%   correction.fitFLIM.expPrompt      The experimental IRF array for all
%                                       SPAD array pixels
%   correction.fitFLIM.timeBin        The bin index range for experimental
%                                       IRF s from all pixels
%   correction.fitFLIM.avgExpPrompt   The average experimental IRF
%   correction.fitFLIM.avgExpPromptFull The average experimental IRF
%                                       spanning the entire useful bin
%                                       range
%   correction.DCR                  A struct with dark count rate maps
%   correction.DRC.raw                Number of dark counts per pixel
%   correction.DRC.threshold          Vector of proportional cutoffs for
%                                       DCR map
%   correction.DRC.index              Sorted DCR value indices
%   correction.DRC.mapXX              Maps of pixels with DCR lower than
%                                       threshold
%
%
%
%   The fields in the minimal version of the corretion struct are:
%   'calibratedBins', 'binWidth', 'globalBinWidth', and 'IRF.peak.Pos'
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
% Last Revision: 14-Apr-2021 - Fix unequal bin width across the array,
%                              Switch to Gaussian IRF model
%                              Get rid os simulated exponentially modified
%                              Gaussian IRF.
% Revision: 11-May-2020 - First version of the file
%
% Copyright 2018-21 Jakub Nedbal
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
                         'Select IRF file(s)', ...
                         'MultiSelect', 'on');

% Keep the file empty if there is no IRF to process
if isequal(file, 0) || isequal(path, 0)
    correction.files.IRF = [];
else
    % Store the IRF file name(s) for later
    % correction.files.IRF = fullfile(path, file);
    correction.files.IRF = strcat(path, file);
end

% Check if correction against different IRF file is needed
ButtonName = ...
    questdlg('IRF Correction on Same IRF or Load a different one?', ...
             'IRF Correction', 'Same', 'Load New', 'Same');
switch ButtonName
    case 'Same'
        correction.files.IRFtest = [];
    case 'Load New'
        % Locate the IRF file
        [file, path] = uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                                  '*.*', 'All Files (*.*)'}, ...
                                 'Select Perormance Test IRF file(s)', ...
                                 'MultiSelect', 'on');

        % Keep the file empty if there is no IRF to process
        if isequal(file, 0) || isequal(path, 0)
            correction.files.IRFtest = [];
        else
            % Store the IRF file name(s) for later
            % correction.files.IRF = fullfile(path, file);
            correction.files.IRFtest = strcat(path, file);
            % load the IRF that will be cross-checked
            correction.IRF.rawTest = ...
                loadSPADdata(correction.files.IRFtest);
        end
    otherwise
        return
end

% Locate the DCR file
[file, path] = uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                          '*.*', 'All Files (*.*)'}, ...
                         'Select DCR file(s)', ...
                         'MultiSelect', 'on');

% Keep the file empty if there is no IRF to process
if isequal(file, 0) || isequal(path, 0)
    correction.files.DCR = [];
else
    % Store the IRF file name(s) for later
    % correction.files.IRF = fullfile(path, file);
    correction.files.DCR = strcat(path, file);
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
correction.repRate = inputdlg('CDM Pulse Repetition Rate [MHz]', ...
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
[correction.CDM, correction.files.CDMfiles] = loadSPADdata(CDMfiles);


%% Analyze the code density map data. Load all the selected files and 
%  analyze the widths of individual histogram bins in every pixel
analyzeCDM;


%% Run IRF correction, if the instrument function is loaded and linearized
if ~isempty(correction.files.IRF)
    % Make a comment
    % fprintf('Loading IRF file %s...\n', correction.files.IRF);
    % Load the IRF file
    % load(correction.files.IRF, 'XYZimage')

    % Linearize the IRF data
    correction.IRF.raw = loadSPADdata(correction.files.IRF);
    % Make a comment
    fprintf(['Linearizing the IRF. This may take a while,\n', ...
             'depending on the number of photons.\n']);
    % Run the linearization routine
    correction.IRF.linear = resampleHistogramPar(correction.IRF.raw);

    % Fit the exponentially modified Gaussian models to the IRF data
    %fitIRF(correction.IRF.linear);
    % Fit Gaussian models to the IRF data
    characterizeIRF(correction.IRF.linear);
    % Correct the bin size by the IRF
    %correctCDMbyIRF;
    % Create parameters for FLIM fitting by SLIM Curve
    fitFLIMparam;
end




%% Load the DCRmap, if the file is available
if ~isempty(correction.files.DCR)
    % Make a comment
    % fprintf('Loading IRF file %s...\n', correction.files.IRF);
    % Load the IRF file
    % load(correction.files.IRF, 'XYZimage')

    % Linearize the IRF data
    correction.DCR.raw = loadSPADdata(correction.files.DCR);
    % Mask only the calibrated bins
    correction.DCR.raw(~correction.calibratedBins) = 0;
    % Sum the number of photons
    correction.DCR.raw = sum(correction.DCR.raw, 3);
    % Set the threshold for DCR
    correction.DCR.threshold = [0.8, 0.85];
    % Sort the DCR values
    [~, correction.DCR.index] = sort(correction.DCR.raw(:));
    % Make a map of low DCR pixels
    %   0 : high DCR
    %   1 : low DCR
    for i = 1 : numel(sort(correction.DCR.threshold(:)))
        mapname = sprintf('map%d', 100 * correction.DCR.threshold(i));
        correction.DCR.(mapname) = true(size(correction.DCR.raw));
        correction.DCR.(mapname)(...
            correction.DCR.index(round(correction.DCR.threshold * ...
                                       numel(correction.DCR.index)) : ...
                                 end)) = false;
    end
end


%% Plot the results, unless there is no need to
if ~isempty(correction.files.graphics)
    analyzeCDMmap(correction.files.graphics, ...
                  fileparts(correction.files.binCorrection));
    analyzeCDMhist(correction.files.graphics, ...
                   fileparts(correction.files.binCorrection));
    analyzeIRFmap(correction.files.graphics, ...
                  fileparts(correction.files.binCorrection));
end



%% Save the calibration data in a file
% Save two files. One adding full. This is a large file with all the data
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
S.correction.globalBinWidth = correction.globalBinWidth;
S.correction.IRF.peak.Pos = correction.IRF.peak.Pos;
S.correction.IRF.fit.goodfit = correction.IRF.fit.goodfit;
S.correction.fitFLIM = correction.fitFLIM;
if isfield(correction, 'DCR')
    for i = 1 : numel(sort(correction.DCR.threshold(:)))
        mapname = sprintf('map%d', 100 * correction.DCR.threshold(i));
        S.correction.DCR.(mapname) = correction.DCR.(mapname);
    end
end
% Make a comment
fprintf(['Saving results into file %s.\n', ...
         'This can take a while. Stand by ...\n'], ...
        correction.files.binCorrection);
% save the results
save(correction.files.binCorrection, '-struct', 'S', '-v7.3')
clear S;
% Reset default figure settings
reset(0)
