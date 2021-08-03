function analyzeCDMmap(graphics, folder)
% ANALYZECDMMAP produces a dynamic figure with map of CDM data. Use the
% mouse to explore the behaviour of individual pixels of the SPAD camera.
%
% Syntax: analyzeCDMmap(graphics)
%
% Inputs:
%   The routine requires a global variable struct 'correction' in the
%   workspace to operate. This struct is produced by running the script
%   SPADcorrection.
%
%   graphics (optional) Can be a string or a cell of strings containing any 
%            of the following 'png', 'eps', 'pdf'. The figure is then saved
%            in the given formats.
%   folder   (optional) The directory into which the resultant images 
%            should be saved.
%
% Outputs:
%   An interactive figure that responds to mouse cursor position and clicks
%   An image file 'folder/SPADlinearization.xxx', where xxx can be one or
%   more of png, eps and pdf. The image shows a snapshot of the CDM map,
%   integral nonlinearity, differential nonlinearity, and Fourier transform
%   of the differential nonlinearity.
%
% Examples:
%   analyzeCDMmap('PNG', pwd)
%
% Toolbox requirement: none
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: SPADcorrection, analyzeCDM, analyzeCDMhist

% Jakub Nedbal
% King's College London
% Aug 2018
% Last revision: 15-Apr-2021 - Changed support for same bin width
% Revision: 08-May-2020 - Tidied up the code, added folder argument
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
global h


%% If no folder is supplied, use the existing one
if nargin < 2
    folder = pwd;
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

%% Create axis for the image of the number of pixels per bin
h.ax.image = axes;
% Place the image of photon count in the axes
h.im = imagesc(correction.avgPhotons);
% Highlight pixel [1, 1]
hold on
h.line.pixel = plot([0 1 1 0 0] + 0.5, [0 0 1 1 0] + 0.5, ...
                    'r-', 'LineWidth', 2);
hold on
% Color map
colormap('default')
% Add a colorbar for scale
h.ax.CB = colorbar;
% label the graph
title('Median Number of Photons per Pixel per Bin')
xlabel('Horizonal Pixels')
ylabel('Vertical Pixels')

%% Label the colorbar
ylabel(h.ax.CB, 'Photons per pixel', 'FontWeight', 'bold', 'FontSize', 8)
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
                 'r-', 'LineWidth', 2);
set(h.ax.CBinvisible, 'YLim', get(h.ax.CB, 'YLim'), ...
                      'Visible', 'off');


%% Plot the integral nonlinearity
h.ax.INL = axes;
h.line.INL = plot(squeeze(correction.INL(1, 1, :)));
title(h.ax.INL, 'Integral Nonlinearity (Pixel [1, 1])')
ylabel(h.ax.INL, 'Error [ps]')

% Set the image
set(h.ax.INL, 'Units', 'centimeters', ...
              'Position', [13, 11, 6, 3], ...
              'Box', 'on', ...
              'FontSize', 8);


%% Plot the differential nonlinearity
h.ax.DNL = axes;
h.line.DNL = plot(squeeze(correction.DNL(1, 1, :)));
title(h.ax.DNL, 'Differential Nonlinearity (Pixel [1, 1])')
ylabel(h.ax.DNL, 'Error [ps]')
xlabel(h.ax.DNL, 'Histogram Bin')

% Set the image
set(h.ax.DNL, 'Units', 'centimeters', ...
              'Position', [13, 6.5, 6, 3], ...
              'Box', 'on', ...
              'FontSize', 8);

%% Plot the Fourier transform of the differential nonlinearity
h.ax.FT = axes;
% Extract the DNL
DNL = squeeze(correction.DNL(4, 4, :));
% Get rid of the NaN values
DNL(isnan(DNL)) = [];
% Get the power spectrum
Ps = abs(fft(DNL) / numel(DNL));
% sampling frequency
Fs = 1 / (correction.avgBinWidth(1, 1)* 1e-12);
% Frequency domain
f = Fs * (0 : (numel(DNL) / 2)) / numel(DNL);
% Convert frequency to GHz
f = f * 1e-9;
% Get the single sided spectrum
Psss = Ps(1 : numel(f));
Psss(2 : end - 1) = 2 * Psss(2 : end - 1);
% Create the plot
h.line.FT = plot(f, Psss);

title(h.ax.FT, 'DNL Fourier Transform (Pixel [1, 1])')
ylabel(h.ax.FT, 'Power Density |P(f)|')
xlabel(h.ax.FT, 'Frequency [GHz]')

% Set the image
set(h.ax.FT, 'Units', 'centimeters', ...
              'Position', [13, 1.5, 6, 3], ...
              'Box', 'on', ...
              'FontSize', 8);


%% Create text labels for describing the data
X = 0.5;
Y = 2;
W = 1.5;
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

X = 7.5;
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

h.lab.Dstd = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'DNL std:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-0.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

h.lab.Istd = uicontrol(h.fig, 'Style', 'text', ...
                              'String', 'INL std:', ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'left');

X = 2;
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

avgP = num2str(correction.avgPhotons(1, 1));
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

Dstd = std(correction.DNL(1, 1, ~isnan(correction.DNL(1, 1, :))));
h.txt.Dstd = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', sprintf('%4.2g ps', Dstd), ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-0.5, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');

Istd = std(correction.INL(1, 1, ~isnan(correction.INL(1, 1, :))));
h.txt.Istd = uicontrol(h.fig, 'Style', 'edit', ...
                              'String', sprintf('%4.2g ps', Istd), ...
                              'Units', 'centimeter', ...
                              'Position', [X, Y-1, W, H], ...
                              'BackgroundColor', 'w', ...
                              'HorizontalAlignment', 'center');


% Give the figure responsiveness to mouse
set(h.fig, 'WindowButtonMotionFcn', @mouseMove, ...
           'WindowButtonDownFcn', @mouseClick);

% Check if graphics is a character, turn it into a cell
if ischar(graphics)
    graphics = {graphics};
end

% A cell with with graphics formats
formats = {'png', '-dpng'; 'eps', '-depsc'; 'pdf', '-dpdf'};

% Go through all requested formats
for i = 1 : numel(graphics)
    % Check if the requested format is one of the allowed ones
    ftype = find(strcmpi(formats(:, 1), graphics{i}));
    % If it is not found throw an error
    assert(~isempty(ftype), 'Graphics type ''%s'' is not recognized.', ...
                            graphics{i})
    % Save the file with the image
    print(h.fig, ...
          fullfile(folder, ...
                   sprintf('SPADlinearization.%s', formats{ftype, 1})),...
          formats{ftype, 2})
end


% end of analyzeCDMmap
end

function mouseMove(~, ~)
global h

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
set(h.txt.Pvar, 'String', num2str(img(C(2), C(1))))
end

function mouseClick(~, ~)
global h
global correction

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
set(h.txt.Pfix, 'String', num2str(img(C(2), C(1))))

% Move the highlighted pixel to wherever the mouse has been clicked
set(h.line.pixel, 'XData', [0 1 1 0 0] + C(1) - 0.5, ...
                  'YData', [0 0 1 1 0] + C(2) - 0.5);

% Gauge on the colorbar according to the median number of photons in the
% selected pixel
set(h.line.bar, 'YData', [1 1] * img(C(2), C(1)));

% Update the DNL and INL stadnard deviation
Dstd = correction.DNL(C(2), C(1), :);
Dstd = std(Dstd(~isnan(Dstd)));
Istd = correction.INL(C(2), C(1), :);
Istd = std(Istd(~isnan(Istd)));
set(h.txt.Dstd, 'String', sprintf('%4.2g ps', Dstd));
set(h.txt.Istd, 'String', sprintf('%4.2g ps', Istd));

% Update the axes titles
title(h.ax.INL, sprintf('Integral Nonlinearity (Pixel [%d, %d])', C))
title(h.ax.DNL, sprintf('Differential Nonlinearity (Pixel [%d, %d])', C))

% Update the INL and DNL graphs
set(h.line.INL, 'YData', squeeze(correction.INL(C(2), C(1), :)));
set(h.line.DNL, 'YData', squeeze(correction.DNL(C(2), C(1), :)));

% Update the Fourier transform
% Extract the DNL
DNL = squeeze(correction.DNL(C(2), C(1), :));
% Get rid of the NaN values
DNL(isnan(DNL)) = [];
% Get the power spectrum
Ps = abs(fft(DNL) / numel(DNL));
% sampling frequency
Fs = 1 / (correction.avgBinWidth(C(2), C(1))* 1e-12);
% Frequency domain
f = Fs * (0 : (numel(DNL) / 2)) / numel(DNL);
% Convert frequency to GHz
f = f * 1e-9;
% Get the single sided spectrum
Psss = Ps(1 : numel(f));
Psss(2 : end - 1) = 2 * Psss(2 : end - 1);
% Create the plot
set(h.line.FT, 'XData', f, ...
               'YData', Psss);
end