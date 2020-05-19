function [fitResult, setting] = processFLIM(setting)
% PROCESSFLIM function used the IRF calibration data to apply minor
% corrections to the CDM linearization data. Essentially, it corrects for
% the fact, that calibration of the bin size, achieved by providing a
% precise known clock during the CDM measurement is affected by the SPAD
% jitter. This routine takes the full-width hal-maximum of the IRF into
% account to correct the bin size.
%
% Syntax: [fitResult, setting] = processFLIM(setting)
%
% Inputs:
%   setting     An optional struct with a selection of parameters. User can
%               provide as many or as little fields. The script will use
%               defaults and ask questions if needed.
%
%   setting.FLIMfile    This is a character array or a cell array of
%                       character arrays with link(s) to file(s) with raw
%                       FLIM transient data measured by the SPAD camera.
%                       The field can be left empty, which will open a
%                       dialog box to select the files for analysis.
%   setting.saveICS     The decays will be saved into an ICS file for
%                       analysis in TRI2. The file name(s) will be same as
%                       the FLIMfile input, except for the extension, which
%                       will be '.ICS'. Default TRUE.
%   setting.saveIRF     The prompt (IRF) will be saved into an ICS file for
%                       analysis in TRI2. The file name(s) will be same as 
%                       the FLIMfile input, except for the extension, which
%                       will be '.IRF.ICS'. The IRF can only be saved if 
%                       the average IRF are selected, i.e. the average 
%                       experimental IRF or the average simulated IRF are 
%                       used. The pixel-specific IRFs cannot be saved 
%                       because TRI2 cannot process them. Default TRUE.
%   setting.runLM       The input data will be fitted with a
%                       Levenberg-Marquardt routine using the SlimCurve
%                       FLIM fitting library. Default TRUE.
%   setting.saveLM      The fitted data in from the Levenberg-Marquardt 
%                       routing will be saved into a MAT file. The file 
%                       name(s) will be same as the FLIMfile input, except 
%                       for the extension, which will be '.FIT.MAT'.
%                       Default TRUE. 
%   setting.fitType     Choose the fitting model for the 
%                       Levenberg-Marquardt routine. The options are 
%                       'exp1', 'exp2', 'exp3', and 'expStretch', for
%                       single-exponential, double-exponential,
%                       triple-exponential, and stretched-exponential
%                       models, respectively. Default 'exp1'.
%   setting.promptType  Choose the type of prompt to use for the
%                       Levenberg-Marquardt fitting and/or to save into the
%                       IRF file. The options are 'avgSimPrompt', 
%                       'simPrompt', 'avgExpPrompt', and 'expPrompt'. These
%                       stand for average simulated IRF, pixel-specific 
%                       simulated IRF, average experimental IRF, and
%                       pixel-specific experimental IRF. 
%                       Default 'avgExpPrompt'.
%   setting.correctionFile  This is a character array with a link to '.MAT'
%                           file containing the 'correction' struct, which 
%                           contains all the information needed for the 
%                           SPAD data correction. This file if produced by
%                           the 'SPADcorrection' script. The 'correction'
%                           struct is kept in the global variable space and
%                           does not need to be repeatedly loaded to save
%                           time. The field can be left empty, which means
%                           that either the global variable 'correction' is
%                           used or a file load dialog will pop up to ask
%                           for the link to the file. Default ''.
%   setting.interactive     When set to TRUE, the input from the 'setting'
%                           struct and the hardwired default values will be
%                           used to run the script without showing the
%                           dialog window to help make the analysis
%                           choices. Default FALSE.
%   setting.overwriteFiles  Is set to TRUE, it will overwrite any files 
%                           that are being produced without giving a 
%                           warning. The script derives file names for 
%                           saving the IRF, the corrected transient data, 
%                           and/or the results of the Levenberg-Marquardt
%                           analysis into files with the same body, but 
%                           different extentions. If these produced files
%                           already exists, it will warn the user and ask 
%                           for an alternative filename. Default FALSE.
%
%   The script tries to use a global variable 'correction' with the
%   calibration data for the specific SPAD array. The array can be
%   overwritten by placing a link to the file with the calibration data
%   into the 'setting.calibrationFile' field. Alternatively, use 
%   'clear global correction', to remove the global variable and the script
%   will ask for the file from which it should load it.
%
% Outputs:
%   fitResult   A struct with the settings used in the run of this script.
%               It has a number of fields, described below. Not all fields
%               are produced each time. It depends on the 'setting' struct
%               field values. The number of elements in 'fitResult' is the
%               same as in 'setting.FLIMfile', i.e. one set of result for
%               each raw SPAD data file analyzed.
%
%   fitResult.fname     Full link to the source file with the raw data from
%                       the SPAD camera.
%   fitResult.transient The corrected FLIM data from the SPAD camera.
%   fitResult.lma_param Result of the fitting parameters as produced by 
%                       SlimCurve.
%   fitResult.lma_fit   Fitted decays produced by SlimCurve.
%   fitResult.offset    Additive constant to the fitted decay.
%   fitResult.chi2      Chi^2 (Chi-squared) error of the fit.
%   fitResult.A         Amplitude of the fitted decay. It can have 1, 2 or
%                       3 elements in the 3rd dimension depending on the
%                       number of exponential elements in the fit.
%   fitResult.tau       Lifetime of the fitted decay. It can have 1, 2 or
%                       3 elements in the 3rd dimension depending on the
%                       number of exponential elements in the fit.
%   fitResult.H         Stretching exponent in the stretched-exponential
%                       fit.
%   fitResult.prompt    The IRF used in the fitting.
%   fitResult.binWidth  Average TDC bin width used in the TRI2 fitting
%                       expressed in nanoseconds
%   fitResult.IRFfile   Full link to the produced IRF-containing file for
%                       use in TRI2.
%   fitResult.IRFfile   Full link to the produced ICS file with the
%                       corrected data from the SPAD.
%   fitResult.LMfile    Full link to the produced MAT file with the
%                       corrected data from the SPAD camera.
%
%   setting     A struct with the settings used in the analysis of the SPAD
%               output data. It has the same fields as the 'setting' struct
%               on the input of the function.
%   
% Examples:
%   processFLIM
%
% Other m-files required: exportICS2, mxSlimCurve, resampleHistogramPar,
%                         saveIRF
% Subfunctions: check_call, fit_call, irf_call, ok_call
% MAT-files required: Selected by 'setting.FLIMfile' field or by 
%                     having a global variable 'correction' loaded
%
% See also: SPADcorrection, mxSlimCurve

% Jakub Nedbal
% King's College London
% Aug 2020
% Last Revision: 19-May-2020 - Created this file. No interactive output yet
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



global correction
% Create a dummy variable to let MATLAB work with nested functions
XYZimage = [];

% Check if setting has not been provided
% Parse the setting, making sure that it is in the proper format
% Check saveICS exists, instruction to save an ICS file with the transient
% This is used to save files for use in TRI2.
if nargin == 0 || ~isfield(setting, 'saveICS')
    % Set an empty struct
    setting.saveICS = true;
end
% Check saveIRF exists, instruction to save an IRF file with the transient
% This is used to save files for use in TRI2
if ~isfield(setting, 'saveIRF')
    setting.saveIRF = true;
end
% Check runLM exists, instruction to run Levenberg-Marquardt fitting of the
% transient data
if ~isfield(setting, 'runLM')
    setting.runLM = true;
end
% Check saveLM exists, instruction to save the results on the 
% Levenberg-Marquardt fitting of the transient into a MATLAB.mat file
if ~isfield(setting, 'saveLM')
    setting.saveLM = true;
end
% Check fitType exists, type of fit used by the L-M fitting engine. Can be
% 'exp1', 'exp2', 'exp3', 'expStretch'
if ~isfield(setting, 'fitType')
    setting.fitType = 'exp1';
end
% Check promptType exists, type of prompts to use. Can be one of the
% following: 'avgSimPrompt', 'simPrompt', 'avgExpPrompt', 'expPrompt'.
if ~isfield(setting, 'promptType')
    setting.promptType = 'avgExpPrompt';
end
% Check correctionFile exists, the full path to the file containing the
% correction struct with all the necessary calibration data.
if ~isfield(setting, 'correctionFile')
    setting.correctionFile = '';
else
    % Check the file actually exists
    assert(exists(setting.correctionFile, 'file') == 2, ...
           'Correction file %s does not exist.', setting.correctionFile);
end
% Check FLIMfile exist, the full file names and paths to files containing 
% the raw SPAD data.
if ~isfield(setting, 'FLIMfile')
    setting.FLIMfile = '';
else
    % Check that it is a string or a cell array of strings
    switch class(setting.FLIMfile)
        case 'char'
            % Check the file actually exists
            assert(exists(setting.FLIMfile, 'file') == 2, ...
                   'FLIM data file %s does not exist.', setting.FLIMfile);
        case 'cell'
            % Check the files actually exists
            for i = 1 : numel(setting.FLIMfile)
                assert(exists(setting.FLIMfile{i}, 'file') == 2, ...
                       'FLIM data file %s does not exist.', ...
                       setting.FLIMfile{i});
            end
        otherwise
            error(['The type of filenames must be a character or ', ...
                   'cell array. You provided type %s.'], ...
                  class(setting.FLIMfile));
    end
end
% Check if interactive is present. It means that a setting window is
% produced. It will have the default items selected
if ~isfield(setting, 'interactive')
    setting.interactive = true;
end
% Check if overwriteFiles is present. It will force overwriting all files.
if ~isfield(setting, 'overwriteFiles')
    setting.overwriteFiles = false;
end

% Set fit parameters
% fitType
%   exp1:      Single exponential
%   exp2:      Double exponential
%   exp3:      Triple exponential
%   expStrech: Stretched exponential (fitting function GCI_stretchedexp)
setting.fitTypes = {'exp1', 'exp2', 'exp3', 'expStretch'};

% Set prompt
% promptType
%   avgSimPrompt: Average simulated prompt same for all pixels.
%   simPrompt:    Simulated prompt different for each pixels. It was
%                 obtained by fitting exponentially modified Gaussian to
%                 the measured prompts, or for poorly fitted pixels with
%                 high DCR by interpolation from neighboring pixels.
%   avgExpPrompt: Average experimental prompt same for all pixels.
%   expPrompt:    Experimental prompt specific to each pixel.
setting.promptTypes = {'avgExpPrompt', 'expPrompt', ...
                       'avgSimPrompt', 'simPrompt'};

% If correction file is provided, load it.
if ~isempty(setting.correctionFile)
    load(setting.correctionFile);
end

% If the correction file link is not provided and correction struct is not
% in the global workspace, load it now.
if isempty(correction)
    uiwait(warndlg(['There is no correction struct. ', ...
                    'Choose a file to load it.'], ...
                    'fitFLIM: Warning'));
    % load the calibration file
    [file, path] = ...
        uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                   '*.*', 'All Files (*.*)'}, ...
                  'Load calibration file');
    % Check if nothing has been returned
    if isequal(file, 0) || isequal(path, 0)
        return
    end
    setting.correctionFile = strcat(path, file);
    load(setting.correctionFile);
end

%% Get the input from the user
% Start by asking a few questions at the beginning
% This way, the code can be left running unattended

if isempty(setting.FLIMfile)
    % Open a dialog box to select files with code density map data
    [file, path] = ...
        uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                   '*.*', 'All Files (*.*)'}, ...
                  'Select files with FLIM measurements', ...
                  'MultiSelect', 'on');
    % Check if nothing has been returned
    if isequal(file, 0) || isequal(path, 0)
        return
    end
    % Combine the path and file to get the full link
    setting.FLIMfile = strcat(path, file);
    % It there was just one file selected and char was returned, turn it
    % into a cell array
    if ischar(setting.FLIMfile)
        setting.FLIMfile = {setting.FLIMfile};
    end
end


if setting.interactive
    %% Create a figure of a setting dialog questionnaire
    h.f = figure('Units', 'Pixels', 'Position', [200, 200, 340, 470], ...
                 'Toolbar', 'None', 'Menu', 'None', 'Name', 'fitFLIM', ...
                 'NumberTitle', 'off');
    % Center the figure
    movegui(h.f, 'center')
    % Create checkboxes, text boxes and radio buttons
    % Save ICS file checkbox
    h.c(1) = uicontrol('Style', 'Checkbox', 'Units', 'Pixels', ...
                       'Value', setting.saveICS, ...
                       'Position', [10, 430, 320, 22], ...
                       'String', 'Save into ICS file', ...
                       'Callback', @check_call);
    % Save IRF checkbox
    h.c(2) = uicontrol('Style', 'Checkbox', 'Units', 'Pixels', ...
                       'Value', setting.saveIRF, ...
                       'Position', [30, 400, 300, 22], ...
                       'String', 'Save accompanying IRF', ...
                       'Callback', @check_call, 'Enable', 'off');
    % Save IRF checkbox comment
    h.t(1) = uicontrol('Style', 'Text', 'Units', 'Pixels', ...
                       'Position', [47, 370, 320, 22], ...
                       'String', '(Requires average IRF)', ...
                       'HorizontalAlignment', 'left', 'Enable', 'off');
    % Run L-M fitting checkbox
    h.c(3) = uicontrol('Style', 'Checkbox', 'Units', 'Pixels', ...
                       'Value', setting.runLM, ...
                       'Position', [10, 340, 320, 22], ...
                       'String', 'Perform L-M fitting', ...
                       'Callback', @check_call);
    % Save L-M fitting checkbox
    h.c(4) = uicontrol('Style', 'Checkbox', 'Units', 'Pixels', ...
                       'Value', setting.saveLM, ...
                       'Position', [30, 310, 320, 22], ...
                       'String', 'Save Fit', ...
                       'Callback', @check_call);
    % Radio buttons for the fit type
    % Single-exponential fit radio button
    h.o(1) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 280, 300, 22], ...
                       'String', 'Single exponential', ...
                       'Callback', @fit_call, ...
                       'Value', strcmp(setting.fitTypes{1}, setting.fitType));
    % Double-exponential fit radio button
    h.o(2) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 250, 300, 22], ...
                       'String', 'Double exponential', ...
                       'Callback', @fit_call, ...
                       'Value', strcmp(setting.fitTypes{2}, setting.fitType));
    % Triple-exponential fit radio button
    h.o(3) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 220, 300, 22], ...
                       'String', 'Triple exponential', ...
                       'Callback', @fit_call, ...
                       'Value', strcmp(setting.fitTypes{3}, setting.fitType));
    % Stretched-exponential fit radio button
    h.o(4) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 190, 300, 22], ...
                       'String', 'Stretched exponential', ...
                       'Callback', @fit_call, ...
                       'Value', strcmp(setting.fitTypes{4}, setting.fitType));
    % Prompt text
    h.t(2) = uicontrol('Style', 'Text', 'Units', 'Pixels', ...
                       'Position', [10, 160, 320, 22], ...
                       'String', 'IRF (prompt) type:', ...
                       'HorizontalAlignment', 'left');
    % Radio buttons for the type of IRF
    % Average experimental IRF radio button
    h.i(1) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 130, 300, 22], ...
                       'String', 'Average experimental IRF', ...
                       'Callback', @irf_call, ...
                       'Value', strcmp(setting.promptTypes{1}, ...
                                       setting.promptType));
    % Pixel-specific experimental IRF radio button
    h.i(2) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 100, 320, 22], ...
                       'String', 'Pixel-specific experimental IRF', ...
                       'Callback', @irf_call, ...
                       'Value', strcmp(setting.promptTypes{2}, ...
                                       setting.promptType));
    % Average simulated IRF radio button
    h.i(3) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 70, 300, 22], ...
                       'String', 'Average simulated IRF', ...
                       'Callback', @irf_call, ...
                       'Value', strcmp(setting.promptTypes{3}, ...
                                       setting.promptType));
    % Pixel-specific simulated IRF radio button
    h.i(4) = uicontrol('Style', 'Radiobutton', 'Units', 'Pixels', ...
                       'Position', [30, 40, 300, 22], ...
                       'String', 'Pixel-specific simulated IRF', ...
                       'Callback', @irf_call, ...
                       'Value', strcmp(setting.promptTypes{4}, ...
                                       setting.promptType));
    % Create OK pushbutton   
    h.p = uicontrol('Style', 'Pushbutton', 'Units', 'Pixels', ...
                    'Position', [140, 10, 60, 20], 'String', 'OK', ...
                    'Callback', @ok_call);

end
    % Pushbutton callback nested function
    function ok_call(varargin)
        % Store the save ICS setting if selected
        setting.saveICS = get(h.c(1), 'Value');
        % Store the save IRF setting if selected and enabled
        setting.saveIRF = h.c(2).Value && isequal(h.c(2).Enable, 'on');
        % Store the run fit setting if selected
        setting.runLM = get(h.c(3), 'Value');
        % Store the fit type setting if selected and enabled
        setting.saveLM = h.c(4).Value && isequal(h.c(4).Enable, 'on');
        % Store the fit type setting
        setting.fitType = setting.fitTypes{[h.o.Value] == 1};
        % Store the prompt type setting
        setting.promptType = setting.promptTypes{[h.i.Value] == 1};
        % Report that the OK button was pressed
        setting.OK = true;
        delete(h.f)
    end

    % Check box click callback nested function
    function check_call(varargin)
        % Check which check box was clicked
        switch find(h.c == varargin{1})
            case 1
                % Save ICS button
                % If average prompt is selected, enable the save IRF button
                % and the checkbox value is positive
                if contains(setting.promptType, 'avg') && ...
                        varargin{1}.Value
                    h.c(2).Enable = 'on';
                    h.t(1).Enable = 'on';
                else
                    h.c(2).Enable = 'off';
                    h.t(1).Enable = 'off';
                end
                if varargin{1}.Value
                    if ~h.c(3).Value
                        set(h.i([1 3]), 'Enable', 'on')
                    end
                else
                end
            case 2
                % Save IRF button
            case 3
                % Run L-M analysis
                % Enable radio buttons for fitting
                if h.c(3).Value
                    set([h.o, h.i, h.c(4)], 'Enable', 'on')
                else
                    if h.c(1).Value
                        set([h.o, h.i([2 4]), h.c(4)], 'Enable', 'off')
                    else
                        set([h.o, h.i, h.c(4)], 'Enable', 'off')
                    end
                end
            case 4
                % Save L-M fit button
        end
    end

    % Check box click callback nested function
    function fit_call(varargin)
        % Make sure the button pressed is selected
        set(varargin{1}, 'Value', 1)
        % Set other radio button false
        set(h.o(h.o ~= varargin{1}), 'Value', 0)
    end

    % Check box click callback nested function
    function irf_call(varargin)
        % Make sure the button pressed is selected
        set(varargin{1}, 'Value', 1)
        % Set other radio button false
        set(h.i(h.i ~= varargin{1}), 'Value', 0)
        % Store the selected radio button
        if contains(setting.promptTypes{[h.i.Value] == 1}, 'avg') && ...
                h.c(1).Value
            % Disable the Save IRF button
            h.c(2).Enable = 'on';
            h.t(1).Enable = 'on';
        else
            % Enable the Save IRF button
            h.c(2).Enable = 'off';
            h.t(1).Enable = 'off';
        end
    end

if setting.interactive
    % Wait for the window to close
    uiwait(h.f)

    % Check if OK was pressed
    if ~setting.OK
        return
    end
end

%% Start the analysis of the SPAD data file-by-file

fitResult = struct([]);
for i = 1 : numel(setting.FLIMfile)
    % Store the file name
    fitResult(i).fname = setting.FLIMfile{i};
    % Load the data
    load(setting.FLIMfile{i})
    % Run the linearization and skew correction routine
    corrHist = resampleHistogramPar(XYZimage);
    % Get the image size
    imSize = size(XYZimage);
    % Check if L-M fitting is required
    if setting.runLM
        % Reshape the data for decay analysis
        transient = reshape(corrHist, imSize(1) * imSize(2), imSize(3))';
        % Select only the part of the curves that should be fitted
        transient = double(transient(correction.fitFLIM.start : ...
                                     correction.fitFLIM.fit_end, :));
    	% Store the transient for later
        fitResult(i).transient = transient;
        % Prepare the prompt. It might be a single prompt for all pixels or
        % a prompt specific to each pixel.
        prompt = correction.fitFLIM.(setting.promptType);
        % The prompt treated depending on the number of dimensions
        if ndims(prompt) == 3
            prompt = reshape(prompt, ...
                             imSize(1) * imSize(2), ...
                             size(prompt, 3))';
        else
            prompt = prompt(:);
        end
        % The bin width, i.e. the linearized TDC time increment 
        x_inc = correction.avgBinWidth(:);
        % The fit start index
        fit_start = correction.fitFLIM.fit_start - ...
                    correction.fitFLIM.start;
        % Convert the fit type into a numerical value
        fit_type = find(contains(setting.fitTypes, setting.fitType));
        % Make a comment
        fprintf('Running the LMA fit on %s...\n', setting.FLIMfile{i})
        % Run the fit
        [fitResult(i).lma_param, ~, fitResult(i).lma_fit] = ...
            mxSlimCurve(transient, ...
                        prompt, ...
                        x_inc, ...
                        fit_start, ...
                        fit_type);
        % Reorganize the LMA fittted curves a rectangular 3D matrix
        % representing the SPAD camera pixels
        fitResult(i).lma_fit = ...
            reshape(fitResult(i).lma_fit', ...
                    imSize(1), ...
                    imSize(2), ...
                    size(fitResult(i).lma_fit, 1));
        % Reorganize the LMA fit results into a rectangular 3D matrix
        % representing the SPAD camera pixels
        fitResult(i).lma_param = ...
            reshape(fitResult(i).lma_param', ...
                    imSize(1), ...
                    imSize(2), ...
                    size(fitResult(i).lma_param, 1));
        % Distribute the results into a more legible format
        % Two of the results are same regardless of the fit type
        % offset
        fitResult(i).offset = squeeze(fitResult(i).lma_param(:, :, 1));
        % chi^2
        fitResult(i).chi2 = squeeze(fitResult(i).lma_param(:, :, end));
        % The other results differ according to the fit type
        switch setting.fitType
            case 'exp1'
                % Single exponential fitting
                % amplitude
                fitResult(i).A = squeeze(fitResult(i).lma_param(:, :, 2));
                % lifetime [ps]
                fitResult(i).tau = squeeze(fitResult(i).lma_param(:,:, 3));
            case 'exp2'
                % Double exponential fitting
                % amplitude 1 and 2
                fitResult(i).A = fitResult(i).lma_param(:, :, [2, 4]);
                % lifetime 1 and 2 [ps]
                fitResult(i).tau = fitResult(i).lma_param(:, :, [3, 5]);
            case 'exp3'
                % Triple exponential fitting
                % amplitude 1, 2, and 3
                fitResult(i).A = fitResult(i).lma_param(:, :, [2, 4, 6]);
                % lifetime 1, 2, and 3 [ps]
                fitResult(i).tau = fitResult(i).lma_param(:, :, [3, 5, 7]);
            case 'expStretch'
                % Stretched exponential fitting
                % amplitude
                fitResult(i).A = squeeze(fitResult(i).lma_param(2));
                % lifetime [ps]
                fitResult(i).tau = squeeze(fitResult(i).lma_param(3));
                % stretching exponent
                fitResult(i).H = squeeze(fitResult(i).lma_param(4));
        end
        % Clean up the transient shape for saving
        fitResult(i).transient = ...
            reshape(fitResult(i).transient', ...
                    imSize(1), ...
                    imSize(2), ...
                    size(fitResult(i).transient, 1));
    end

    %% Save the IRF if this was requested
    if setting.saveIRF
        % Prepare the prompt. It will be a single prompt for all pixel made
        % either by experimental prompt averaging or by an exponentially
        % modified Gaussian model
        fitResult(i).prompt = ...
            correction.fitFLIM.([setting.promptType 'Full']) ...
                (1 : correction.fitFLIM.fit_end);
        % The prompt needs to be three dimensional
        fitResult(i).prompt = reshape(fitResult(i).prompt, ...
                                      [1 1 numel(fitResult(i).prompt)]);
        % binWidth is the average bin width across the range, expressed in
        % nanoseconds
        fitResult(i).binWidth = mean(correction.avgBinWidth(:)) * 1e-3;
        % Create the file name to save the data
        [path, file] = fileparts(setting.FLIMfile{i});
        fitResult(i).IRFfile = fullfile(path, [file, '.irf.ics']);
        % Check that the file doesn't exist
        if exist(fitResult(i).IRFfile, 'file') == 2 && ...
                ~setting.overwriteFiles
            % Give a warning that the file already exists
            uiwait(warndlg({sprintf('File %s already exists. ', ...
                                    fitResult(i).IRFfile), ...
                            'Choose a filename to save the IRF.'}, ...
                           'fitFLIM: Warning'));
            % Ask for a new file name or confirm the original file
            [file, path] = ...
                uiputfile({'*.irf.ics', 'MAT-files (*.irf.ics)'; ...
                           '*.*',   'All Files (*.*)'}, ...
                          'Save as', fitResult(i).IRFfile);
            % If box is closed or Cancel pressed, do not save
            if ischar(file)
                % Set the new filename
                fitResult(i).IRFfile = strcat(path, file);
                % Save the file
                saveIRF(fitResult(i).prompt, ...
                        [1 1], ...
                        fitResult(i).IRFfile, ...
                        fitResult(i).binWidth);
            else
                % File is 'none', because Cancel was pressed
                fitResult(i).IRFfile = 'none';
            end
        else
            % If the file doesn't exist yet, save it directly
            saveIRF(fitResult(i).prompt, ...
                    [1 1], ...
                    fitResult(i).IRFfile, ...
                    fitResult(i).binWidth);
        end
        % Clean up the prompt shape for saving
        fitResult(i).prompt = fitResult(i).prompt(:);
    end

    %% Save the ICS file for TRI2 analysis if this was requested
    if setting.saveICS
        % Reshape the data for decay analysis
        % Select only the part of the curves that should be fitted
        transient = double(corrHist(:, :, 1 : correction.fitFLIM.fit_end));
        % Next, check if the range has not been already been calculated, 
        % e.g. the IRF has not been saved
        if ~setting.saveIRF
            % binWidth is the average bin width across the range expressed
            % in nanoseconds
            fitResult(i).binWidth = mean(correction.avgBinWidth(:)) * 1e-3;
        end
        % Range is the time range of the prompt expressed in seconds
        range = fitResult(i).binWidth * size(transient, 3) * 1e-9;
        % Create the file name to save the data
        [path, file] = fileparts(setting.FLIMfile{i});
        fitResult(i).ICSfile = fullfile(path, [file, '.ics']);
        % Check that the file doesn't exist
        if exist(fitResult(i).ICSfile, 'file') == 2 && ...
                ~setting.overwriteFiles
            % Give a warning that the file already exists
            uiwait(warndlg({sprintf('File %s already exists. ', ...
                                    fitResult(i).ICSfile), ...
                            'Choose filename to save data for TRI2.'}, ...
                           'fitFLIM: Warning'));
            % Ask for a new file name or confirm the original file
            [file, path] = uiputfile({'*.ics', 'MAT-files (*.ics)'; ...
                                      '*.*',   'All Files (*.*)'}, ...
                                     'Save as', fitResult(i).ICSfile);
            % If box is closed or Cancel pressed, do not save
            if ischar(file)
                % Set the new filename
                fitResult(i).ICSfile = strcat(path, file);
                % Save the file
                exportICS2(transient, fitResult(i).ICSfile, range);
            else
                % File is 'none', because Cancel was pressed
                fitResult(i).ICSfile = 'none';
            end
        else
            % If the file doesn't exist yet, save it directly
            exportICS2(transient, fitResult(i).ICSfile, range);
        end
    end

    %% Check if the results should be saved
    if setting.saveLM
        % Reorder the transient so that it is in a form of a pixel array of
        % TDC histograms
        
        % Create the file name to save the data
        [path, file] = fileparts(setting.FLIMfile{i});
        fitResult(i).LMfile = fullfile(path, [file, '.fit.mat']);
        % Check that the file doesn't exist
        if exist(fitResult(i).LMfile, 'file') == 2 && ...
                ~setting.overwriteFiles
            % Give a warning that the file already exists
            uiwait(warndlg({sprintf('File %s already exists. ', ...
                                    fitResult(i).LMfile), ...
                            'Choose a filename to save the data.'}, ...
                           'fitFLIM: Warning'));
            % Ask for a new file name or confirm the original file
            [file, path] = ...
                uiputfile({'*.fit.mat', 'MAT-files (*.fit.mat)'; ...
                           '*.*',   'All Files (*.*)'}, ...
                          'Save as', fitResult(i).LMfile);
            % If box is closed or Cancel pressed, do not save
            if ischar(file)
                % Set the new filename
                fitResult(i).LMfile = strcat(path, file);
                % Create a temporary struct for saving. This is to save
                % only the result of this iterration processing iterration
                S.fitResult = fitResult(i);
                S.setting = setting;
                % Save the file
                save(fitResult(i).LMfile, '-struct', 'S');
                clear S;
            else
                % File is 'none', because Cancel was pressed
                fitResult(i).LMfile = 'none';
            end
        else
            % If the file doesn't exist yet, save it directly
            % Create a temporary struct for saving. This is to save only
            % the result of this iterration processing iterration
            S.fitResult = fitResult(i);
            S.setting = setting;
            % Save the file
            save(fitResult(i).LMfile, '-struct', 'S');
            clear S;
        end
    end
end
    
end
