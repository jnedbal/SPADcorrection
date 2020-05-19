global correction

% Set fit parameters
% fit_type
%   1: Single exponential
%   2: Double exponential
%   3: Triple exponential
%   4: Stretched exponential (fitting function GCI_stretchedexp)
fit_type = 1;

% Set prompt
% prompt_type
%   avgSimPrompt: Average simulated prompt same for all pixels.
%   simPrompt:    Simulated prompt different for each pixels. It was
%                 obtained by fitting exponentially modified Gaussian to
%                 the measured prompts, or for poorly fitted pixels with
%                 high DCR by interpolation from neighboring pixels.
%   avgExpPrompt: Average experimental prompt same for all pixels.
%   expPrompt:    Experimental prompt specific to each pixel.
prompt_type = 'simPrompt';

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
    corrFile = strcat(path, file);
    load(corrFile);
end

%% Get the input from the user
% Start by asking a few questions at the beginning
% This way, the code can be left running unattended

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
FLIMfiles = strcat(path, file);

fitResult = struct([]);
for i = 1 : numel(FLIMfiles)
    % Store the file name
    fitResult(i).fname = FLIMfiles{i};
    % Load the data
    load(FLIMfiles{i})
    % Run the linearization and skew correction routine
    corrHist = resampleHistogramPar(XYZimage);
    % Get the image size
    imSize = size(XYZimage);
    % Reshape the data for decay analysis
    corrHist = reshape(corrHist, imSize(1) * imSize(2), imSize(3))';
    % Select only the part of the curves that should be fitted
    transient = double(corrHist(correction.fitFLIM.start : ...
                                correction.fitFLIM.fit_end, :));
	% Store the transient for later
    fitResult(i).transient = transient;
    % Prepare the prompt. It might be a single prompt for all pixels or a
    % prompt specific to each pixel.
    prompt = correction.fitFLIM.(prompt_type);
    % The prompt needs to be treated depending on its number of dimensions
    if ndims(prompt) == 3
        prompt = reshape(prompt, imSize(1) * imSize(2), size(prompt, 3))';
    else
        prompt = prompt(:);
    end
    % The bin width, i.e. the linearized TDC time increment 
    x_inc = correction.avgBinWidth(:);
    % The fit start index
    fit_start = correction.fitFLIM.fit_start - correction.fitFLIM.start;
    % Make a comment
    fprintf('Running the LMA fit on %s...\n', FLIMfiles{i})
    % Run the fit
    [fitResult(i).lma_param, ~, fitResult(i).lma_fit] = ...
        mxSlimCurve(transient, ...
                    prompt, ...
                    x_inc, ...
                    fit_start, ...
                    fit_type);
    % Distribute the results into a more legible format
    switch fit_type
        case 1
            % Single exponential fitting
            fitResult(i).A = fitResult(i).lma_param(2);     % amplitude
            fitResult(i).tau = fitResult(i).lma_param(3);   % lifetime [ps]
            fitResult(i).offset = fitResult(i).lma_param(1);% offset
            fitResult(i).chi2 = fitResult(i).lma_param(4);  % chi^2
        case 2
            % Double exponential fitting
            fitResult(i).A(1) = fitResult(i).lma_param(2);  % amplitude 1
            fitResult(i).tau(1) = fitResult(i).lma_param(3);% lifetime [ps]
            fitResult(i).A(2) = fitResult(i).lma_param(4);  % amplitude 2
            fitResult(i).tau(2) = fitResult(i).lma_param(5);% lifetime [ps]
            fitResult(i).offset = fitResult(i).lma_param(1);% offset
            fitResult(i).chi2 = fitResult(i).lma_param(6);  % chi^2
    end
end
    
    %