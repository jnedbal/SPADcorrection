function displayFLIM(fitResult)
% displayFLIM displays linearized annd corrected SPAD image data.
% It produces an interactive figure to display all the fit parameters, the
% microtime decays, fit residuals, IRF, and fit results.
%
% Syntax: fitIRF(fitResult)
%
% Inputs:
%   fitResult       Struct or a character array.
%                   * Struct fitResult is produced by the function 
%                     processFLIM.
%                   * Character array with the link to the (.fit.mat) file
%   Alternatively, no input has to be provided. Instead, a popup dialog box
%   asking for the file is going to appear.
%
% Outputs:
%   An interactive figure to play with.
%
% Examples:
%   displayFLIM(fitResult)
%
% Other m-files required: createPointer
% Subfunctions: edit_Callback, pop_Callback, keyhit, mouseclick, mousemove,
%               mouserelease
%               
% MAT-files required: none, but a link to the fit.mat file can be provided
%
% See also: processFLIM, SPADcorrection

% Jakub Nedbal
% King's College London
% May 2020
% Last Revision:  26-May-2020 - First working prototype
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



%% Process the input
% Check if nothing has been supplied
if nargin == 0
    % Open a dialog box to select files with the processed FLIM data
    [file, path] = ...
        uigetfile({'*.fit.mat', 'FLIM Data Files (*.fit.mat)'; ...
                   '*.*', 'All Files (*.*)'}, ...
                  'Select file with processed FLIM measurements');
    % Check if nothing has been returned
    if isequal(file, 0) || isequal(path, 0)
        return
    end
    fitResult = fullfile(path, file);
end

% Check if the input is a struct of a character array with the link to the
% source file containing the corrected FLIM data
switch class(fitResult)
    case 'char'
        % Character array, presumably containing the link to the file
        assert(exist(fitResult, 'file') == 2, ...
               'The file provided ''%s'' does not exist.', fitResult);
        % Load the file
        load(fitResult, 'fitResult');
    case 'struct'
        % Struct has been provided. Assume that it is OK
    otherwise
        error(['Expecting character array or struct, ', ...
               'received ''%s'' instead.'], class(fitResult));
end

%% Create titles for lifetime fits
% Number of fit parameters
display.nrNames = size(fitResult.A, 3) + ...
                  size(fitResult.tau, 3) + ...
                  size(fitResult.offset, 3) + ...
                  size(fitResult.chi2, 3);
% For stretched exponential, there is also stretch exponent H
if isfield(fitResult, 'H')
    display.nrNames = display.nrNames + 1;
end
display.order = zeros(1, display.nrNames);
display.names = cell(size(display.order));
display.graphNames = cell(size(display.order));
% Offset
display.names{1} = 'Z:';
display.graphNames{1} = 'Offset Z';
display.order(1) = numel(display.order) - 1;
% Chi^2
display.names{end} = [char([967, 178]), ':'];
display.order(end) = numel(display.order);
display.graphNames{end} = sprintf('Chi-squared %s', char([967, 178]));
% Amplitudes and lifetimes
if numel(display.names) == 4 || numel(display.names) == 5
    % Single-exponential situation
    display.names{2} = 'A:';
    display.order(2) = 1;
    display.graphNames{2} = 'Amplitude A';
    display.names{3} = [char(964), ':'];
    display.order(3) = 2;
    display.graphNames{3} = ...
        sprintf('Fluorescence Lifetime %s [ps]', char(964));
else
    % Double- or triple-exponential situation
    for j = 0 : 2 : numel(display.names) - 3
        display.names{j + 2} = ['A' char(8321 + j / 2) ':'];
        display.order(j + 2) = 1 + j;
        display.graphNames{j + 2} = ['Amplitude A', char(8321 + j / 2)];
        display.names{j + 3} = [char([964, 8321 + j / 2]) ':'];
        display.order(j + 3) = 2 + j;
        display.graphNames{j + 3} = ...
            sprintf('Fluorescence Lifetime %s [ps]', ...
                    char([964, 8321 + j / 2]));
    end
end
% Stretched exponential model
if isfield(fitResult, 'H')
    display.names{end - 1} = 'H:';
    display.order(end - 1) = numel(display.order) - 2;
    display.graphNames{end - 1} = 'Stretch Exponent H';
end
% Create and index of the names in the order
[~, display.index] = sort(display.order);
% Create the selected index for the parameter that is being plotted
display.selIndex = 3;

%% Set the raw image as the default to display
display.rawInterp = 'lma_param';

%% Work out the time histogram bins and limits
display.Xbins.exp = fitResult.binWidth * ...
                    (0 : (fitResult.fitFLIM.fit_end - 1));
display.Xlim.exp = display.Xbins.exp([1, end]);
display.Xbins.fit = fitResult.binWidth * ...
           (fitResult.fitFLIM.start : fitResult.fitFLIM.fit_end);
display.Xlim.fit = display.Xbins.fit([1, end]);
display.Xbins.fitResid = fitResult.binWidth * ...
           (fitResult.fitFLIM.fit_start : fitResult.fitFLIM.fit_end);
display.Xlim.fitResid = display.Xbins.fitResid([1, end]);


%% Work out the images to display
display.img.intensity = flipud(sum(fitResult.transient, 3));
display.img.fit = flipud(squeeze(fitResult.(display.rawInterp)(:, :, 3)));
display.img.resid = ...
    double(fitResult.transient(:, :, fitResult.fitFLIM.fit_start : end));
display.img.resid = display.img.resid - ...
    fitResult.lma_fit(:, :, (end - size(display.img.resid, 3) + 1) : end);
display.img.resid = display.img.resid ./ ...
    abs(sqrt(fitResult.lma_fit(:, :, ...
                       (end - size(display.img.resid, 3) + 1 : end))));
display.img.xlim = [0, size(fitResult.offset, 2)] + 0.5;
display.img.ylim = [0, size(fitResult.offset, 1)] + 0.5;
display.img.xpixLim = [0, size(fitResult.offset, 2) - 1];
display.img.ypixLim = [0, size(fitResult.offset, 1) - 1];

%% Work out the limits on the scales for each fit parameter
display.goodfit.mask = repmat(fitResult.fitFLIM.goodfit, ...
                              1, 1, size(fitResult.(display.rawInterp),3));
% Only keep the parameters from good pixels
display.goodfit.(display.rawInterp) = ...
    fitResult.(display.rawInterp)(display.goodfit.mask);
% Reshape the parameters to have a parameter per column
display.goodfit.(display.rawInterp) = ...
    reshape(display.goodfit.(display.rawInterp), ...
            numel(display.goodfit.(display.rawInterp)) / ...
                size(fitResult.(display.rawInterp), 3), ...
            size(fitResult.(display.rawInterp), 3));
display.goodfit.median = median(display.goodfit.(display.rawInterp));
display.goodfit.std = std(display.goodfit.(display.rawInterp));
display.goodfit.minLim = display.goodfit.median - 3 * display.goodfit.std;
display.goodfit.maxLim = display.goodfit.median + 3 * display.goodfit.std;

%% Display the image
close all
% Get screen size
sSize = get(0, 'ScreenSize');
% Prepare figure size
fSize = [(sSize(3) - sSize(4)) / 2, sSize(4) * 0.1, sSize(4) * [1, 0.8]];
% Prepare axes size
aSize = [fSize(3) * 0.07, fSize(3) * 0.06, ...
         min([0.4 * fSize([4, 4]); 0.8 * fSize([3, 3])])];
% Create a figure
h.fig = figure('Units', 'pixels', 'Position', fSize, ...
               'Name', 'displayFLIM', 'NumberTitle', false);
% Center the figure
movegui(h.fig, 'center')
% Create a pointer for the figure
h.fig.PointerShapeCData = createPointer;
h.fig.PointerShapeHotSpot = [16 16];

% Add some text elements
pos = [aSize(1) + aSize(3) * 1.025, aSize(2), 30, 20];
h.text.lab.Xpix = uicontrol('Style', 'text', ...
                            'String', 'X:', ...
                            'Position', pos + [0, 30, 0, 0]);
h.text.lab.Ypix = uicontrol('Style', 'text', ...
                            'String', 'Y:', ...
                            'Position', pos + [0, 0, 0, 0]);
h.text.val.Xpix = uicontrol('Style', 'edit', ...
                            'String', '', ...
                            'Position', pos + [pos(3), 30, 100, 0], ...
                            'BackgroundColor', [1 1 1], ...
                            'Callback', @edit_Callback, ...
                            'Tag', 'Xpix');
h.text.val.Ypix = uicontrol('Style', 'edit', ...
                            'String', '', ...
                            'Position', pos + [pos(3), 0, 100, 0], ...
                            'BackgroundColor', [1 1 1], ...
                            'Callback', @edit_Callback, ...
                            'Tag', 'Ypix');


%% Create starting position for lifetime fits
pos = pos + [0, 60, 0, 0];

for j = 1 : numel(display.names)
    % i is the counting index
    tpos = pos + [0, 30 * (numel(display.names) - display.order(j)), 0, 0];
    h.text.lab.fit(j) = uicontrol('Style', 'text', ...
                                  'String', display.names(j), ...
                                  'Position', tpos);
    tpos = tpos + [tpos(3), 0, 100, 0];
    h.text.val.fit(j) = uicontrol('Style', 'text', ...
                                  'String', num2str(j), ...
                                  'Position', tpos, ...
                                  'BackgroundColor', [1 1 1]);
end

%% Place four axes there
% Axes for intensity image
h.axes(1) = axes('Units', 'pixels', 'Position', aSize);
colormap(h.axes(1), 'gray')

% Axes for lifetime fit
aSize = repmat(aSize, 6, 1);
aSize(2, 1) = aSize(1, 1) + 0.57 * fSize(3);
aSize(2, 4) = aSize(1, 4) * 0.6;
h.axes(2) = axes('Units', 'pixels', ...
                 'Position', aSize(2, :), ...
                 'XLim', display.Xlim.exp);
% Axes for residuals
aSize(3, :) = aSize(2, :);
aSize(3, 2) = aSize(3, 2) + aSize(3, 4) * 7 / 6;
aSize(3, 4) = aSize(3, 4) * 0.5;
h.axes(3) = axes('Units', 'pixels', ...
                 'Position', aSize(3, :), ...
                 'XLim', display.Xlim.exp);

% Axes for lifetime image
h.axes(4) = axes('Units', 'pixels');
aSize(4, :) = aSize(3, :);
aSize(4, 4) = aSize(4, 4) / 0.5 / 0.6;
aSize(4, 2) = aSize(4, 2) + aSize(4, 4) * 0.5 - 20;
h.axes(4).Position = aSize(4, :);

% Axes for lifetime histogram
h.axes(6) = axes('Units', 'pixels');
aSize(6, 2) = aSize(4, 2) + aSize(4, 4) * 0.2;
aSize(6, 4) = aSize(4, 4) * 0.8;
h.axes(6).Position = aSize(6, :);
% Add x-axis label
xlabel(display.graphNames{display.selIndex});


%% Draw the maximum projection (intensity)
axes(h.axes(1))
h.image(1) = imagesc(display.img.intensity);
hold on
% Flip the Y axis, add the box
h.axes(1).YDir = 'reverse';
h.axes(1).Box = 'on';
h.axes(1).XAxisLocation = 'top';

% Draw two four dummy lines later used as a cross hair
for j = 1 : 4
    h.line.ax1(j) = plot([0 0], [0 0], 'g-');
end

% Draw the decay transient
axes(h.axes(2))
hold on
% draw dummy lines on axes 2
% This will be the decay
h.line.ax2(1) = plot(display.Xbins.exp, zeros(size(display.Xbins.exp)));

% Draw a line for start of the transient data
h.line.ax2(3) = ...
    plot([1, 1] * fitResult.fitFLIM.start * fitResult.binWidth, ...
         h.axes(2).YLim);
% Draw a line for start of the fit
h.line.ax2(4) = ...
    plot([1, 1] * fitResult.fitFLIM.fit_start * fitResult.binWidth, ...
    h.axes(2).YLim, 'Color', h.line.ax2(3).Color);
% This will be the fit
h.line.ax2(2) = plot(display.Xbins.fit, zeros(size(display.Xbins.fit)), ...
                     'LineStyle', '--', 'LineWidth', 2);
% Set the X limits
h.axes(2).XLim = display.Xlim.exp;
% Set the box around
h.axes(2).Box = 'on';
% Add the X- and Y-axis labels
xlabel('Microtime [ns]')
ylabel('Photon Count')

%% Draw a dummy line of residuals
axes(h.axes(3))
% This will be the line of residuals
h.line.ax3(1) = plot(display.Xbins.fitResid, ...
                     zeros(size(display.Xbins.fitResid)), '.');
hold on
% Draw a line for start of the transient data
h.line.ax3(2) = plot(h.line.ax2(3).XData, h.axes(3).YLim);
% Draw a line for start of the fit
h.line.ax3(3) = plot(h.line.ax2(4).XData, h.axes(3).YLim, ...
                     'Color', h.line.ax3(2).Color);
% Give it a Y axis label
ylabel('Residue')
% Draw the IRF on the right axis
yyaxis right
% if there is one common IRF for all pixels, plot it
h.line.ax3(4) = plot(fitResult.binWidth * fitResult.fitFLIM.timeBin, ...
                     zeros(size(fitResult.fitFLIM.timeBin)));
if isvector(fitResult.fitFLIM.prompt)
    h.line.ax3(4).YData = fitResult.fitFLIM.prompt;
end
% Set the X limits
h.axes(2).XLim = display.Xlim.exp;
% Add the Y-axis label
ylabel('IRF')
yticks([])


% Draw the lifetime map
axes(h.axes(4))
h.image(2) = ...
    imagesc(fitResult.(display.rawInterp)(:, :, display.selIndex));
hold on
cmap = parula(256);
cmap(end, :) = [1 0 0]; % Make highest value red
cmap(1, :) = [0 0 0];   % Make lowest value black
colormap(h.axes(4), cmap)
h.axes(5) = colorbar('westoutside');
h.axes(4).Position = aSize(4, :);
h.axes(5).Units = 'Pixels';
aSize(5, :) = h.axes(5).Position;
aSize(5, 1) = aSize(4, 1) - aSize(4, 3) * 0.3;
% pos2(1) = pos1(1) - pos1(3) * 0.3;
aSize(5, 2) = aSize(4, 2) + aSize(4, 4) * 0.15;
% pos2(2) = pos1(2) + pos1(4) * 0.15;
aSize(5, 4) = aSize(4, 4) * 0.7;
h.axes(5).Position = aSize(5, :);
% Label the colorbar
set(h.axes(5).Label, 'String', display.graphNames{display.selIndex})

% Invisible axes for lines over colorbar
h.axes(7) = axes('Units', 'pixels');
h.axes(7).Position = h.axes(5).Position;
% Plot the lower bound line
h.line.ax7(1) = plot([0, 1], [0 0], 'Color', [0, 0, 0], 'LineWidth', 2);
hold on
% Plot the upper bound line
h.line.ax7(2) = plot([0, 1], [2 2], 'Color', [1, 0, 0], 'LineWidth', 2);
% Plot the sliding value line
h.line.ax7(3) = plot([0, 1], [1 1], 'Color', [1, 1, 1], 'LineWidth', 2);
h.axes(7).XLim = [0, 1];
h.axes(7).YLim = h.axes(5).YLim;
h.axes(7).Visible = 'off';

% Change the colormap to grayscale
% Flip the Y axis, add the box
h.axes(4).YDir = 'reverse';
h.axes(4).Box = 'on';
h.axes(4).XAxisLocation = 'top';
h.axes(4).XLim = display.img.xlim;
h.axes(4).YLim = display.img.ylim;



% Draw two four dummy lines later used as a cross hair
for j = 1 : 4
    h.line.ax4(j) = plot([0 0], [0 0], 'g-');
end

%% Add a popupmenu (combo box) to choose linear or logarithmic plot
pos = [aSize(2, [1, 2]) + aSize(2, [3, 4]) - [80, 40], 70, 30];
h.pop(1) = uicontrol('Style', 'popupmenu', ...
                     'String', {'Lin', 'Log'}, ...
                     'Position', pos, ...
                     'Callback', @pop_Callback);
%% Add a popupmenu (combo box) for choosing the fit parameter
pos = [aSize(4, 1) + aSize(4, 3) - 70, aSize(4, 2) - 35, 70, 30];
h.pop(2) = uicontrol('Style', 'popupmenu', ...
                     'String', cellfun(@(x) x(1 : end - 1), ...
                                       display.names(display.index), ...
                                       'UniformOutput', false), ...
                     'Position', pos, ...
                     'Value', display.order(display.selIndex), ...
                     'Callback', @pop_Callback);
%% Add a popupmenu (combo box) for choosing the fit parameter
pos = pos + [-170, 0, 80, 0];
h.pop(3) = uicontrol('Style', 'popupmenu', ...
                     'String', {'Raw', 'Interpolated'}, ...
                     'Position', pos, ...
                     'Value', 1, ...
                     'Callback', @pop_Callback);
% If the lmaInterp field does not exist in the fitResult struct, make the
% above popummenu disabled
if ~isfield(fitResult, 'lmaInterp')
    h.pop(3).Enable = 'off';
end


%% Add lifetime limits
pos = [aSize(5, 1) + 0.5 * aSize(5, 3) - 35, aSize(4, 2), 70, 26];
h.spin(1) = uicontrol('Style', 'edit', ...
                      'Position', pos, ...
                      'Callback', @edit_Callback, ...
                      'Tag', 'lowTau');
pos(2) = aSize(4, 2) + aSize(4, 4) - 26;
h.spin(2) = uicontrol('Style', 'edit', ...
                    'Position', pos, ...
                    'Callback', @edit_Callback,...
                    'Tag', 'highTau');
% Set the default limits
display.CLim = [display.goodfit.minLim(3), display.goodfit.maxLim(3)];

% Set the histogram limits. They are 15 % more on each side
display.HLim = display.CLim + [-0.15, 0.15] * diff(display.CLim);
% Choose the data for the histogram
display.Hdata = fitResult.(display.rawInterp)(:, :, display.selIndex);

% Only select the data from within the histogram
display.Hdata = display.Hdata(display.Hdata >= display.HLim(1) & ...
                              display.Hdata <= display.HLim(2));
% histogram bins
display.Hbins = linspace(display.HLim(1), display.HLim(2), 100);
% Plot the histogram
axes(h.axes(6))
h.hist = histogram(display.Hdata, display.Hbins);
hold on
% plot the axes of the histogram
h.line.ax6(1) = plot(display.CLim(1) * [1, 1], h.axes(6).YLim);
h.line.ax6(2) = plot(display.CLim(2) * [1, 1], h.axes(6).YLim, ...
                     'Color', h.line.ax6(1).Color);
% Set the axis limits
h.axes(6).XLim = display.HLim;

%set(h.spin, 'SliderStep', [0.05, 0.05], ...
%            'Min', tauLim(1), ...
%            'Max', tauLim(2))
% round all spinner numbers to 3 significant digits
round3 = @(x) round(x, 3, 'significant');
% Add the spinner values
h.spin(1).String = num2str(round3(display.CLim(1)));
h.spin(2).String = num2str(round3(display.CLim(2)));
%set(h.spin(2), 'Value', tauLim(2), ...
%               'String', num2str(round(tauLim(2), 3, 'significant')))

% Update the lifetime graph
edit_Callback(h.spin(1))

% Set all the font sizes to be same for all figures
set(h.axes, 'FontSize', 10);

%% Make the figure responsive to the mouse and keyboard
% Call 'mousemove' function when mouse moves over the window
h.fig.WindowButtonMotionFcn = @mousemove;

% Call 'mouseclick' function when mouse is clicked over the window
h.fig.WindowButtonDownFcn = @mouseclick;

% Call 'mouserelease' function when mouse is clicked over the window
h.fig.WindowButtonUpFcn = @mouserelease;

% Call 'keyhit' function when key is pressed
h.fig.KeyPressFcn = @keyhit;

% Create flags for mouse actions
% Flag for stopping the hairpin movement
display.stopped = false;
% Flag for noting a click
display.clicked = false;
% Index of the line being moved
display.movedLine = 0;
    
    % Popupmenu callback - lin/log scale
    function pop_Callback(handle, ~)
        % Figure out which 
        switch handle
            case h.pop(1)   % Lin/log popup has been called
                switch h.pop(1).String{h.pop(1).Value}
                    % Asking for linear scale
                    case 'Lin'
                        h.axes(2).YScale = 'linear';
                    % Asking for a logarithmic scale
                    case 'Log'
                        h.axes(2).YScale = 'log';
                end
            case {h.pop(2), h.pop(3)}   % Type of data to be displayed
                % Check whether raw or interpolated data should be shown
                if h.pop(3).Value == 2
                    display.rawInterp = 'lmaInterp';
                else
                    display.rawInterp = 'lma_param';
                end
                % Find the index of the data to be displayed
                display.selIndex = display.index(h.pop(2).Value);

                % Update the image in the graph
                h.image(2).CData = flipud(fitResult.(display.rawInterp) ...
                    (:, :, display.selIndex));

                % Update the limits on the range
                display.CLim = ...
                    [display.goodfit.minLim(display.selIndex), ...
                     display.goodfit.maxLim(display.selIndex)];

                % Add the spinner values
                h.spin(1).String = num2str(round3(display.CLim(1)));
                h.spin(2).String = num2str(round3(display.CLim(2)));

                % Update the lifetime graph
                edit_Callback(h.spin(1))

                % Update the colormap axis name
                set(h.axes(5).Label, ...
                    'String', display.graphNames{display.selIndex});

                % Update the graph
                updateGraphs;
        end
    end

    % Function to check that the mousepointer is inside axes
    function index = inAxes
        % Index of existing axes 
        index = false(1, 6);
        % Check if the mouse is in either the 1st or 4th axes
        % Retrieve the coordinate of the mouse in 1st, 4th or 6th axes
        for i = [1 4 6]
            coord = get(h.axes(i), 'CurrentPoint');
            % Keep only the useful bit
            coord = coord(1, [1 2]);
            % Check if pointer is outside of the limits
            outX = any(diff([h.axes(i).XLim(1), ...
                             coord(1), ...
                             h.axes(i).XLim(2)]) < 0);
            outY = any(diff([h.axes(i).YLim(1), ...
                             coord(2), ...
                             h.axes(i).YLim(2)]) < 0);
            if ~any([outX, outY])
                index(i) = true;
            end
        end
    end

    % Mouse button up callback function
    function mouserelease(~, ~)
        % Flag that the mouse button is not held pressed
        display.clicked = false;

        % Check if the release is done after moving the line
        if display.movedLine > 0
            % Call the edit callback button to update the histogram
            edit_Callback(h.spin(display.movedLine));
            % Make sure the moved line is forgotten after button release
            display.movedLine = 0;
        end
    end

    % Mouse button down callback function
    function mouseclick(~, ~)
        % Check if the mouse is in either the 1st, 4th or 6th axes        
        index = inAxes;

        % End the function if it is outside in any of the directions
        if ~any(index)
            return
        end

        % If in the 1st or 4th axes
        if index(1) || index(4)
            display.stopped = ~display.stopped;
        end

        % If in the 6 axes (histogram)
        if index(6)
            display.clicked = true;
        end
        

        mousemove(false, false);
    end

    function mousemove(~, ~)
        % Check if the mouse is in either the 1st, 4th or 6th axes        
        index = inAxes;

        % If in the 1st or 4th axes
        if index(1) || index(4)
            % Set a camera focus pointer
            h.fig.Pointer = 'custom';
            % If the motion is not stopped, update the graph
            if ~display.stopped
                % Update the coordinates
                display.img.coord = get(h.axes(index), 'CurrentPoint');
                % Update the graph
                updateGraphs
            end
            % end the function
            return
        elseif index(6)
            % Get the coordinates in the hsitogram axes
            display.hist.coord = get(h.axes(index), 'CurrentPoint');
            % Get the distance and the index of the closest vertical bar
            [M, I] = ...
                min(abs(display.hist.coord(1) - display.CLim));
            %if display.movedLine == 0
            % Remember which line is being moved
            display.movedLine = I;
            %end
            % Check if the minimum is not more than 0.5% of the histogram
            % limits distance
            if M < diff(display.CLim) * 0.01
                h.fig.Pointer = 'cross';
            else
                % Make sure the pointer is an arrow
                h.fig.Pointer = 'arrow';
            end
            if display.clicked && display.movedLine > 0
                % Update the graph limit
                display.CLim(display.movedLine) = display.hist.coord(1);
                % Update the position of the line
                h.line.ax6(display.movedLine).XData = ...
                    display.hist.coord(1) * [1 1];
                % Update the number in the colorbar limit
                h.spin(display.movedLine).String = ...
                    num2str(round3(display.hist.coord(1)));
            end
            % end the function
            return
        end

        % Make sure the pointer is an arrow
        h.fig.Pointer = 'arrow';
    end

    % Function called when a key is hit
    function keyhit(~, event)
        switch event.Key
            case 'rightarrow'
                % Increase the existing X value by one
                xpixval = str2double(h.text.val.Xpix.String) + 1;
                % Check if the value is within the limits of the image
                if xpixval <= display.img.xpixLim(2)
                    % Update the edit box value
                    h.text.val.Xpix.String = num2str(xpixval);
                    % Update the graph
                    edit_Callback(h.text.val.Xpix);
                end
            case 'leftarrow'
                % Decrease the existing X value by one
                xpixval = str2double(h.text.val.Xpix.String) - 1;
                % Check if the value is within the limits of the image
                if xpixval >= display.img.xpixLim(1)
                    % Update the edit box value
                    h.text.val.Xpix.String = num2str(xpixval);
                    % Update the graph
                    edit_Callback(h.text.val.Xpix);
                end
            case 'uparrow'
                % Increase the existing Y value by one
                ypixval = str2double(h.text.val.Ypix.String) - 1;
                % Check if the value is within the limits of the image
                if ypixval >= display.img.ypixLim(1)
                    % Update the edit box value
                    h.text.val.Ypix.String = num2str(ypixval);
                    % Update the graph
                    edit_Callback(h.text.val.Ypix);
                end
            case 'downarrow'
                % Decrease the existing Y value by one
                ypixval = str2double(h.text.val.Ypix.String) + 1;
                % Check if the value is within the limits of the image
                if ypixval <= display.img.ypixLim(2)
                    % Update the edit box value
                    h.text.val.Ypix.String = num2str(ypixval);
                    % Update the graph
                    edit_Callback(h.text.val.Ypix);
                end
        end
    end

    
    % Editable text boxes feedback
    function edit_Callback(handle, ~)
        % First check that the value makes sense, that it is a number
        value = str2double(handle.String);
        if isnan(value)
            return
        end
        % Run the code according to the edit box pressed
        switch handle.Tag
            case 'Xpix'
                % Check that the value is within the range
                if value < 0 || value > display.img.xpixLim(2)
                    return
                end
                % Update the coordinate
                display.img.coord(1) = value + 1;
                % Update the graphs
                updateGraphs;
            case 'Ypix'
                % Check that the value is within the range
                if value < 0 || value > display.img.ypixLim(2)
                    return
                end
                % Update the coordinate
                display.img.coord(3) = value + 1;
                % Update the graphs
                updateGraphs;
            case {'lowTau', 'highTau'}
                % First check that the value is a number
                value = [str2double(h.spin(1).String), ...
                         str2double(h.spin(2).String)];
                if any(isnan(value))
                    % Revert the value to the original
                    h.spin(1).String = num2str(round3(display.CLim(1)));
                    h.spin(2).String = num2str(round3(display.CLim(2)));
                    return
                end

                % Check which of the two values have changed
                if value(1) == display.CLim(1) && ...
                        value(2) == display.CLim(2)
                    % If none has changed, just leave
                    return
                end

                % Check if value(1) is smaller than value(2)
                if value(1) < value(2)
                    % Update the limits according to the values in the box
                    display.CLim = round3(value);
                    h.axes(4).CLim = display.CLim;
                end

                % Add the rounded text values bax to the editable boxes
                h.spin(1).String = num2str(round3(display.CLim(1)));
                h.spin(2).String = num2str(round3(display.CLim(2)));

                % Update the histogram graph
                % Set the histogram limits. They are 15 % more on each side
                display.HLim = ...
                    display.CLim + [-0.15, 0.15] * diff(display.CLim);
                % Choose the data for the histogram
                display.Hdata = ...
                    fitResult.(display.rawInterp)(:, :, display.selIndex);

                % Only select the data from within the histogram
                display.Hdata = ...
                    display.Hdata(display.Hdata >= display.HLim(1) & ...
                                  display.Hdata <= display.HLim(2));
                % histogram bins
                display.Hbins = ...
                    linspace(display.HLim(1), display.HLim(2), 100);
                % Shrink the lines so they do not affect the y limit scale
                h.line.ax6(1).YData = [0, 0];
                h.line.ax6(2).YData = [0, 0];
                % Update the axes limits
                h.axes(6).XLim = display.HLim;
                % Update the histogram
                h.hist.Data = display.Hdata;
                h.hist.BinEdges = display.Hbins;
                % plot the axes of the histogram
                h.line.ax6(1).XData = display.CLim(1) * [1, 1];
                h.line.ax6(1).YData = h.axes(6).YLim;
                h.line.ax6(2).XData = display.CLim(2) * [1, 1];
                h.line.ax6(2).YData = h.axes(6).YLim;
                % Set the axis limits
                h.axes(6).XLim = display.HLim;
                %h.hist = histogram(data);
                %h.axes(6).XLim = value;
                set(h.axes(6).XLabel, ...
                    'String', display.graphNames{display.selIndex});

                % Update the lines over the Colorbar
                h.line.ax7(1).YData = display.CLim(1) * [1, 1];
                h.line.ax7(2).YData = display.CLim(2) * [1, 1];
                h.axes(7).YLim = display.CLim;
        end
    end

    % Function that updates the graph on mouse move or X, Y pixel box
    % change
    function updateGraphs
        % Draw the cursor
        ydiff = 0.05 * diff(display.img.ylim);
        h.line.ax1(1).XData = display.img.coord([1 1]);
        h.line.ax1(1).YData = ...
            [display.img.coord(3) + ydiff, display.img.ylim(2)];
        h.line.ax1(2).XData = display.img.coord([1 1]);
        h.line.ax1(2).YData = ...
            [display.img.ylim(1), display.img.coord(3) - ydiff];
        xdiff = 0.05 * diff(display.img.xlim);
        h.line.ax1(3).XData = ...
            [display.img.coord(1) + xdiff, display.img.xlim(2)];
        h.line.ax1(3).YData = display.img.coord([3 3]);
        h.line.ax1(4).XData = ...
            [display.img.xlim(1), display.img.coord(1) - xdiff];
        h.line.ax1(4).YData = display.img.coord([3 3]);

        for i = 1 : 4
            h.line.ax4(i).XData = h.line.ax1(i).XData;
            h.line.ax4(i).YData = h.line.ax1(i).YData;
        end
        % Round the pixel number, subtract 1 to start from zero
        C = round(display.img.coord);

        % Display the X and Y pixel positions
        h.text.val.Xpix.String = num2str(C(1) - 1);
        h.text.val.Ypix.String = num2str(C(3) - 1);

        % Calculate the image row and column index from mouse position
        rowIn = display.img.ypixLim(2) + 2 - C(3);
        colIn = C(1);
        % Display the fit values. This can be either interpolated or raw.
        % When interpolated, show both
        % The raw value is always good, when goodfit is true or when raw
        % data is selected
        if fitResult.fitFLIM.goodfit(rowIn, colIn) || ...
                isequal(display.rawInterp, 'lma_param')
            for i = 1 : numel(h.text.val.fit)
                % Just show the raw fitted value
                h.text.val.fit(i).String = ...
                    num2str(round3(fitResult.lma_param ...
                                   (rowIn, colIn, i)));
            end
        else
            % Show both the interpolated and the non-interpolated value
            for i = 1 : numel(h.text.val.fit)
                % Show the interpolated fitted value followed by the
                % raw value in parentheses
                h.text.val.fit(i).String = ...
                    sprintf('%g (%g)', ...
                            round3(fitResult.lmaInterp ...
                                           (rowIn, colIn, i)), ...
                            round3(fitResult.lma_param ...
                                           (rowIn, colIn, i)));
                
            end
        end

        % Draw the decay and fit
        set(h.line.ax2([3, 4]), 'YData', [0, 0])
        h.line.ax2(1).YData = squeeze(fitResult.transient(rowIn, colIn,:));
        yfitData = squeeze(fitResult.lma_fit(rowIn, colIn, :));
        yfitData(yfitData < max(1, double(min(h.line.ax2(1).YData)))) =NaN;
        h.line.ax2(2).YData = yfitData;
        %h.axes(2).YLimMode = 'auto';
        %h.axes(2).YLim(1) = max(1, min(h.line.ax2(1).YData));
        set(h.line.ax2([3, 4]), 'YData', h.axes(2).YLim)

        %set(h.axes(2), 'XLim', display.Xlim.exp)

        % Draw the residuals
        set(h.line.ax3([2, 3]), 'YData', [0, 0])
        h.line.ax3(1).YData = squeeze(display.img.resid(rowIn, colIn, :));
        set(h.line.ax3([2, 3]), 'YData', h.axes(3).YAxis(1).Limits)
        % Draw the IRF
        if ndims(fitResult.fitFLIM.prompt) == 3
            h.line.ax3(4).YData = fitResult.fitFLIM.prompt(rowIn, colIn,:);
        end

        % Draw the colorbar line
        h.line.ax7(3).YData = fitResult.(display.rawInterp)(rowIn, colIn, display.selIndex) * [1, 1];
    end
end





