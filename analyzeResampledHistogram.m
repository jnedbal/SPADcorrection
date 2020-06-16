function analyzeResampledHistogram(XYZimage, graphics)
% ANALYZERESAMPLEDHISTOGRAM produces a dynamic figure with map of CDM data
%
% Syntax: analyzeResampleHistogram(XYZimage, graphics)
%
% Inputs:
%   XYZimage can be one of three inputs (1) raw CDM measurement data array,
%            (2) a string with a link to a file containing the raw CDM
%            measurement data array, (3) nothing at all.
%
%   graphics can is a string or a cell array of strings with one or more of
%            the following 'png', 'eps' or 'pdf'. Depending on the choice,
%            it will save the output figures in one or more data formats.
%            If graphics is not supplied, or contains an empty string '',
%            no file will be saved.
%
% Outputs:
%   Image files 'SPADresampled.xxx' and 'SPAD.stDev.xxx', where xxx can be 
%   one or more of png, eps and pdf. The image shows a snapshot of the CDM 
%   map, code density map of the original data, resampled data and ideal
%   synthetic data. Then there is a graph of the 
%   integral nonlinearity, differential nonlinearity, and Fourier transform
%   of the differential nonlinearity.
%
% Examples:
%   analyzeResampleHistogram
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: 

% Jakub Nedbal
% King's College London
% Aug 2018
% Last Revision: 04-Jun-2020
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

% If the correction file link is not provided and correction struct is not
% in the global workspace, load it now.
if isempty(correction) || ~isfield(correction, 'avgPhotons')
    uiwait(warndlg(['Full correction struct is missing. ', ...
                    'Choose a file to load it.'], ...
                    'Process FLIM: Warning'));
    % load the calibration file
    [file, path] = ...
        uigetfile({'*.full.mat', 'MATLAB Files (*.full.mat)'; ...
                   '*.*', 'All Files (*.*)'}, ...
                  'Load calibration file');
    % Check if nothing has been returned
    if isequal(file, 0) || isequal(path, 0)
        return
    end
    correctionFile = strcat(path, file);
    load(correctionFile); %#ok<LOAD>
end

%% Input data can be a link to a file, a time-resolved image matrix with a 
%  code density map or empty.
% If there is no input data, make it empty
if nargin == 0
    XYZimage = [];
end

switch class(XYZimage)
    case 'char'
        % This looks like a link to a file
        % Test this file actually exists
        assert(exist(XYZimage, 'file') == 2, ...
               'File %s does not exist.', XYZimage)
        % Now that it looks like it exists, load it
        load(XYZimage); %#ok<LOAD>
    case 'uint32'
        % This looks like a code density map
        % leave it as it is
    otherwise
        % Looks like it is some gibberish and so we need ask for the file.
        
        % load the data file
        [file, path] = ...
            uigetfile({'*.mat', 'MATLAB Files (*.mat)'; ...
                       '*.*', 'All Files (*.*)'}, ...
                      'Load input data file');
        % Check if nothing has been returned
        if isequal(file, 0) || isequal(path, 0)
            return
        end
        inputFile = strcat(path, file);
        load(inputFile); %#ok<LOAD>
end

%% Check if there is any graphics input
if nargin < 2
    graphics = '';
end

%% Run the linearization
corrHist = resampleHistogramPar(XYZimage);

% convert values to doubles to make working on them possible
corrHist = double(corrHist);
XYZimage = double(XYZimage);

%% calculate the standard deviations of the raw and corrected CDM 
%  normalized by the square root of their sums

rawSTD = zeros(size(corrHist, 1), size(corrHist, 2));
corrSTD = zeros(size(corrHist, 1), size(corrHist, 2));
for i = 1 : size(rawSTD, 1)
    for j = 1 : size(rawSTD, 2)
        in = correction.firstBin(i, j) : correction.lastBin(i, j);
        % Get rid of the first and last few indices, which tend to be off
        in([1 : 3, end - 2 : end]) = [];
        rawSTD(i, j) = ...
            std(XYZimage(i, j, in)) / sqrt(sum(XYZimage(i, j, in)));
        in = find(corrHist(i, j, :), 1) : find(corrHist(i, j, :),1,'last');
        % Get rid of the first and last few indices
        in([1 : 3, end - 2 : end]) = [];
        corrSTD(i, j) = ...
            std(corrHist(i, j, in)) / sqrt(sum(corrHist(i, j, in)));
    end
end

%% Create a figure to plot the result
close all
h.fig = figure;
set(h.fig, 'Units', 'centimeters', ...
           'Position', [1 1 20 15], ...
           'PaperUnits', 'centimeters', ...
           'PaperPosition', [0 0 20 15], ...
           'PaperPositionMode', 'auto', ...
           'Color', 'w');

% Create axis for the image of the number of pixels per bin
h.ax.image = axes;
% Place the image of photon count in the axes
h.im = imagesc(correction.avgPhotons);
% Highlight pixel [1, 1]
hold on
h.line.pixel = plot([0 1 1 0 0] + 0.5, [0 0 1 1 0] + 0.5, ...
                    'w-', 'LineWidth', 2);
hold on
% Color map
colormap('default')
% Add a colorbar for scale
h.ax.CB = colorbar;
% label the graph
title('Median Number of Photons per Pixel per Bin')
xlabel('Horizonal Pixels')
ylabel('Vertical Pixels')
% Label the colorbar
ylabel(h.ax.CB, 'Median Photons per Bin', 'FontWeight', 'bold', ...
       'FontSize', 8)
% Set the images
set(h.ax.image, 'Units', 'centimeters', ...
                'Position', [1.5, 5, 8, 8], ...
                'Box', 'on', ...
                'XAxisLocation', 'top', ...
                'FontSize', 8);

% Create fake invisible axes over colorbar
h.ax.CBinvisible = axes;
% Set the axes so it matches the colorbar
set(h.ax.CBinvisible, 'Units', get(h.ax.CB, 'Units'), ...
                      'Position', get(h.ax.CB, 'Position'));

% Highlight the number of photons in pixel [1, 1] in the colorbar
h.line.bar = plot(h.ax.CBinvisible, ...
                 [0 1], [1 1] * correction.avgPhotons(1, 1), ...
                 'w-', 'LineWidth', 2);
set(h.ax.CBinvisible, 'YLim', get(h.ax.CB, 'YLim'), ...
                      'Visible', 'off');

% Plot the histograms
h.ax.hist = axes;
Xdata = squeeze(find(correction.calibratedBins(1, 1, :)));
% Get rid of the last few bins due to bin size change
Xdata(end - 1 : end) = [];
Xdata = [Xdata, ...
         (1 : numel(Xdata))' + find(corrHist(1, 1, :), 1) - 1, ...
         (1 : numel(Xdata))'];
Ydata = [squeeze(XYZimage(1, 1, Xdata(:, 1))), ...
         squeeze(corrHist(1, 1, Xdata(:, 2)))];
% Normalize the Ydata
Ydata = Ydata ./ mean(Ydata);
h.line.hist = plot(Xdata(:, 3), Ydata);

h.ax.hist.XLim = Xdata([1 end], 3);

title(h.ax.hist, 'Normalized Code Density Map (Pixel [1, 1])')
ylabel(h.ax.hist, 'Photon Probability Densisty')
xlabel(h.ax.hist, 'TDC Bin Index')
h.leg.hist = legend('Raw', 'Corrected', ...
                    'Location', 'SouthOutside', ...
                    'Orientation', 'Horizontal');

% Set the image
set(h.ax.hist, 'Units', 'centimeters', ...
               'Position', [13, 9, 6, 5], ...
               'Box', 'on', ...
               'FontSize', 8);

% Plot the Fourier transform of the differential nonlinearity
h.ax.FTcdm = axes;
% Get the power spectrum
Ps = abs(fft(Ydata - 1) / size(Ydata, 1));
% sampling frequency
Fs = 1 / (correction.avgBinWidth(1, 1)* 1e-12);
% Frequency domain
f = Fs * (0 : (size(Ydata, 1) / 2)) / size(Ydata, 1);
% Convert frequency to GHz
f = f * 1e-9;
% Get the single sided spectrum
Psss = Ps(1 : numel(f), :);
Psss(2 : end - 1, :) = 2 * Psss(2 : end - 1, :);
% Create the plot
h.line.FTcdm = plot(f, Psss);
h.ax.FTcdm.XLim = f([1 end]) + [-0.05, 0.05] * (f(end) - f(1));
h.ax.FTcdm.XTick = 0 : 2 : 12;
h.ax.FTcdm.XGrid = 'on';

title(h.ax.FTcdm, 'CDM Fourier Transform (Pixel [1, 1])')
ylabel(h.ax.FTcdm, 'Power Density |P(f)|')
xlabel(h.ax.FTcdm, 'Frequency [GHz]')

% Set the image
set(h.ax.FTcdm, 'Units', 'centimeters', ...
                'Position', [13, 1.5, 6, 5], ...
                'Box', 'on', ...
                'FontSize', 8);

% Create text labels for describing the data
X = 0.5;
Y = 2;
W = 2.0;
H = 0.5;
h.lab.Xvar = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'X:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.Yvar = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'Y:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+0.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.Pvar = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'Phot/Bin:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

X = 7.0;
h.lab.Xfix = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'X:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.Yfix = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'Y:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+0.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.Pfix = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'Phot/Bin:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.std = uicontrol(h.fig,  'Style', 'text', ...
                              'String', 'CDM Standard Deviation:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-0.5, 2.5*W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.raw = uicontrol(h.fig,  'Style', 'text', ...
                              'String', 'Raw:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.corr = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'Corrected:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-1.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

% h.lab.syn = uicontrol(h.fig,  'Style', 'text', ...
%                               'String', 'Synthetic:', ...
%                               'Units', 'centimeter', ...
%                               'Position', [X, Y-1.5, W, H], ...
%                               'BackgroundColor', 'w', ...
%                               'HorizontalAlignment', 'left');

X = 2.2;
h.txt.Xvar = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', '1', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

h.txt.Yvar = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', '1', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+0.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

avgP = num2str(round(correction.avgPhotons(1, 1)));
h.txt.Pvar = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', avgP, ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

X = 9;
h.txt.Xfix = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', '1', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

h.txt.Yfix = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', '1', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y+0.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

h.txt.Pfix = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', avgP, ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

rawstd = num2str(round(rawSTD(1, 1), 2, 'significant'));
h.txt.raw = uicontrol(h.fig, 'Style', 'edit', ...
                             'String', rawstd, ...
                             'Units', 'centimeter', ...
                             'Position', [X, Y-1, W, H], ...
                             'BackgroundColor', 'w', ...
                             'HorizontalAlignment', 'center', ...
                             'ForegroundColor', h.line.FTcdm(1).Color);

corrstd = num2str(round(corrSTD(1, 1), 2, 'significant'));
h.txt.corr = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', corrstd, ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-1.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center', ...
                             'ForegroundColor', h.line.FTcdm(2).Color);

% synstd = std(Ydata(:, 3));
% h.txt.syn = uicontrol(h.fig,  'Style', 'edit', ...
%                               'String', sprintf('%4.2g', synstd), ...
%                               'Units', 'centimeter', ...
%                               'Position', [X, Y-1.5, W, H], ...
%                               'BackgroundColor', 'w', ...
%                               'HorizontalAlignment', 'center');


% Give the figure responsiveness to mouse
set(h.fig, 'WindowButtonMotionFcn', @mouseMove, ...
           'WindowButtonDownFcn', @mouseClick);

% Check if graphics is a character, turn it into a cell
if ischar(graphics)
    if isempty(graphics)
        graphics = {};
    else
        graphics = {graphics};
    end
end



%% Produce a figure of standard deviations on the CDM.

% calculate the standard deviations of the raw and corrected CDM normalized
% by the square root of their sums

%% Add the figure for the standard deviations of the CDMs
h.fig(2) = figure;
h.fig(2).Units = 'centimeters';
h.fig(2).Position = [1 1 20 7];
h.fig(2).PaperUnits = 'centimeters';
h.fig(2).PaperPosition = [0 0 20 7];
h.fig(2).PaperPositionMode = 'auto';
h.fig(2).Color = [1, 1, 1];

% Create axis for the image of the number of pixels per bin
h.ax.corrSTD = axes;
% Place the image of photon count in the axes
h.imCorrSTD = surf(corrSTD);
% Color map
colormap('default')
% Add a colorbar for scale
h.ax.CBstd = colorbar;
% label the graph
title('Normalized Standard Deviation of Corrected CDM')
xlabel('Horizontal Pixels')
ylabel('Vertical Pixels')
% Label the colorbar
ylabel(h.ax.CBstd, 'Standard Deviation', 'FontWeight', 'bold', ...
       'FontSize', 8)
% Set the images
h.ax.corrSTD.Units = 'centimeters';
h.ax.corrSTD.Position = [2, 1, 5, 5];
h.ax.corrSTD.Box = 'on';
h.ax.corrSTD.XAxisLocation = 'top';
h.ax.corrSTD.FontSize = 8;
h.ax.corrSTD.ZScale = 'log';
h.ax.corrSTD.XLim = [0, size(corrSTD, 2)];
h.ax.corrSTD.XTick = [h.ax.corrSTD.XTick, h.ax.corrSTD.XLim(2)];
h.ax.corrSTD.YLim = [0, size(corrSTD, 1)];
h.ax.corrSTD.YTick = [h.ax.corrSTD.YTick, h.ax.corrSTD.YLim(2)];
h.ax.corrSTD.XLabel.Rotation = 20;
h.ax.corrSTD.XLabel.Position(1) = h.ax.corrSTD.XLim(2) * 0.5;
h.ax.corrSTD.XLabel.Position(2) = -h.ax.corrSTD.YLim(2) * 0.25;
h.ax.corrSTD.XLabel.Position(3) = h.ax.corrSTD.ZLim(1);
h.ax.corrSTD.XLabel.HorizontalAlignment = 'center';
h.ax.corrSTD.YLabel.Rotation = -30;
h.ax.corrSTD.YLabel.Position(3) = h.ax.corrSTD.ZLim(1);
h.ax.corrSTD.YLabel.Position(2) = h.ax.corrSTD.YLim(2) * 0.5;
h.ax.corrSTD.YLabel.Position(1) = -h.ax.corrSTD.XLim(2) * 0.25;
h.ax.corrSTD.YLabel.HorizontalAlignment = 'center';

h.imCorrSTD.EdgeColor = 'none';

% Create axis for the image of the number of pixels per bin
h.ax.rawSTD = axes;
% Place the image of photon count in the axes
h.imRawSTD = surf(rawSTD);
% Color map
colormap('default')
% Add a colorbar for scale
h.ax.CBstdR = colorbar;
% label the graph
title('Normalized Standard Deviation of Raw CDM')
xlabel('Horizontal Pixels')
ylabel('Vertical Pixels')
% Label the colorbar
ylabel(h.ax.CBstdR, 'Standard Deviation', 'FontWeight', 'bold', ...
       'FontSize', 8)
% Set the images
h.ax.rawSTD.Units = 'centimeters';
h.ax.rawSTD.Position = [12, 1, 5, 5];
h.ax.rawSTD.Box = 'on';
h.ax.rawSTD.XAxisLocation = 'top';
h.ax.rawSTD.FontSize = 8;
h.ax.rawSTD.ZScale = 'log';
h.ax.rawSTD.XLim = [0, size(corrSTD, 2)];
h.ax.rawSTD.XTick = [h.ax.rawSTD.XTick, h.ax.rawSTD.XLim(2)];
h.ax.rawSTD.YLim = [0, size(corrSTD, 1)];
h.ax.rawSTD.YTick = [h.ax.rawSTD.YTick, h.ax.rawSTD.YLim(2)];
h.ax.rawSTD.XLabel.Rotation = 20;
h.ax.rawSTD.XLabel.Position(1) = h.ax.rawSTD.XLim(2) * 0.5;
h.ax.rawSTD.XLabel.Position(2) = -h.ax.rawSTD.YLim(2) * 0.25;
h.ax.rawSTD.XLabel.Position(3) = h.ax.corrSTD.ZLim(1);
h.ax.rawSTD.XLabel.HorizontalAlignment = 'center';
h.ax.rawSTD.YLabel.Rotation = -30;
h.ax.rawSTD.YLabel.Position(3) = h.ax.corrSTD.ZLim(1);
h.ax.rawSTD.YLabel.Position(2) = h.ax.rawSTD.YLim(2) * 0.5;
h.ax.rawSTD.YLabel.Position(1) = -h.ax.rawSTD.XLim(2) * 0.25;
h.ax.rawSTD.YLabel.HorizontalAlignment = 'center';

h.imRawSTD.EdgeColor = 'none';
% Make same color limits on both graphs
zlim = [min(h.ax.rawSTD.ZLim(1), h.ax.corrSTD.ZLim(1)), ...
        max(h.ax.rawSTD.ZLim(2), h.ax.corrSTD.ZLim(2))];
h.ax.rawSTD.ZLim = zlim;
h.ax.corrSTD.ZLim = zlim;


%% SAve the figures

% A cell with with graphics formats
formats = {'png', '-dpng'; 'eps', '-depsc'; 'pdf', '-dpdf'};

fnames = {'SPADresampled.%s', 'SPADstDev.%s'};
for j = 1 : numel(fnames)
    % Go through all requested formats
    for i = 1 : numel(graphics)
        % Check if the requested format is one of the allowed ones
        ftype = find(strcmpi(formats(:, 1), graphics{i}));
        % If it is not found throw an error
        assert(~isempty(ftype), ...
               'Graphics type ''%s'' is not recognized.', graphics{i})
        % Save the file with the image
        print(h.fig(j), ...
              sprintf(fnames{j}, formats{ftype, 1}),...
              formats{ftype, 2})
    end
end


%% Mose callback functions
function mouseMove(~, ~)

    % Get the position of the mouse cursor
    C = round(get(h.ax.image, 'CurrentPoint'));
    C = C(1, [1 2]);
    % Get the pixel indices for the array
    xl = get(h.ax.image, 'XLim');
    yl = get(h.ax.image, 'YLim');

    % If the cursor is outside of the pixel array, exit the function
    if  C(1) < xl(1) || C(1) > xl(2) || C(2) < yl(1) || C(2) > yl(2)
        return
    end

    % Update the X pixel coordinate
    set(h.txt.Xvar, 'String', C(1))
    % Update the Y pixel coordinate
    set(h.txt.Yvar, 'String', C(2))

    % Get the median photons per bin value
    img = get(h.im, 'CData');
    % Write the median number of photons per pixel per bin
    h.txt.Pvar.String = num2str(round(img(C(2), C(1))));
end

function mouseClick(~, ~)

    % Get the position of the mouse cursor
    C = round(get(h.ax.image, 'CurrentPoint'));
    C = C(1, [1 2]);
    % Get the pixel indices for the array
    xl = get(h.ax.image, 'XLim');
    yl = get(h.ax.image, 'YLim');

    % If the cursor is outside of the pixel array, exit the function
    if  C(1) < xl(1) || C(1) > xl(2) || C(2) < yl(1) || C(2) > yl(2)
        return
    end

    % Update the X pixel coordinate
    set(h.txt.Xfix, 'String', C(1))
    % Update the Y pixel coordinate
    set(h.txt.Yfix, 'String', C(2))

    % Get the median photons per bin value
    img = get(h.im, 'CData');
    % Write the median number of photons per pixel per bin
    h.txt.Pfix.String = num2str(round(img(C(2), C(1))));

    % Move the highlighted pixel to wherever the mouse has been clicked
    set(h.line.pixel, 'XData', [0 1 1 0 0] + C(1) - 0.5, ...
                      'YData', [0 0 1 1 0] + C(2) - 0.5);

    % Gauge on the colorbar according to the median number of photons in the
    % selected pixel
%     Xdata = squeeze(find(correction.calibratedBins(C(2), C(1), :)));
%     Ydata = [squeeze(inputHist(C(2), C(1), Xdata)), ...
%              squeeze(corrHist(C(2), C(1), Xdata))];
    Xdata = squeeze(find(correction.calibratedBins(C(2), C(1), :)));
    % Get rid of the last few bins due to bin size change
    Xdata(end - 2 : end) = [];
    Xdata = [Xdata, ...
             (0 : (numel(Xdata)-1))' + find(corrHist(C(2), C(1), :), 1),...
             (1 : numel(Xdata))'];
    Ydata = [squeeze(XYZimage(C(2), C(1), Xdata(:, 1))), ...
             squeeze(corrHist(C(2), C(1), Xdata(:, 2)))];
    % Normalize the Ydata
    Ydata = Ydata ./ mean(Ydata);

    set(h.line.hist(1), 'XData', Xdata(:, 3), 'YData', Ydata(:, 1));
    set(h.line.hist(2), 'XData', Xdata(:, 3), 'YData', Ydata(:, 2));
    %set(h.line.hist(3), 'XData', Xdata, 'YData', Ydata(:, 3));
    h.ax.hist.XLim = Xdata([1 end], 3);

    % Update the standard deviations of the CDMs
    h.txt.raw.String = ...
        num2str(round(rawSTD(C(2), C(1)), 2, 'significant'));
    h.txt.corr.String = ...
        num2str(round(corrSTD(C(2), C(1)), 2, 'significant'));
    %set(h.txt.syn, 'String', sprintf('%4.2g', std(Ydata(:, 3))));

    % Update the axes titles
    title(h.ax.hist, ...
        sprintf('Normalized Code Density Map (Pixel [%d, %d])', C))
    title(h.ax.FTcdm, sprintf('CDM Fourier Transform (Pixel [%d, %d])', C))

    % Update the Fourier transform
    Ps = abs(fft(Ydata - 1) / size(Ydata, 1));
    % sampling frequency
    Fs = 1 / (correction.avgBinWidth(C(2), C(1)) * 1e-12);
    % Frequency domain
    f = Fs * (0 : (size(Ydata, 1) / 2)) / size(Ydata, 1);
    % Convert frequency to GHz
    f = f * 1e-9;
    % Get the single sided spectrum
    Psss = Ps(1 : numel(f), :);
    Psss(2 : end - 1, :) = 2 * Psss(2 : end - 1, :);
    % Create the plot
    set(h.line.FTcdm(1), 'XData', f(2:end), 'YData', Psss((2:end), 1));
    set(h.line.FTcdm(2), 'XData', f(2:end), 'YData', Psss((2:end), 2));
    % set(h.line.FTcdm(3), 'XData', f(2:end), 'YData', Psss((2:end), 3));

    % Update the graph limits
    h.ax.FTcdm.XLim = f([1 end]) + [-0.05, 0.05] * (f(end) - f(1));
    
    % Update the bar
    h.line.bar.YData = [1 1] * correction.avgPhotons(C(2), C(1));
    end

end