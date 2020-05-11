function analyzeCDMhist(graphics, folder)
% ANALYZECDMHIST produces histograms of pixel TDC characterisics
%
% Syntax: analyzeCDMhist(graphics, folder)
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
%   Image files '/path/to/TDCbinWidth.xxx', 'path/to/TDCrange.xxx', and
%   'path/to/TDCcalibration.xxx' where xxx can be one or more of png, eps,
%   and pdf. The images show histograms of average bin width per pixel, the
%   number of TDC bins in the period of the calibration clock, and an
%   example of the TDC histogram for one pixel, highlighting the
%   calibration clock period.
%
% Examples:
%   analyzeCDMhist('PNG', pwd)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: SPADcorrection, analyzeCDM, analyzeCDMmap

% Jakub Nedbal
% King's College London
% Aug 2018
% Last Revision: 08-May-2020 - Tidied up the code, added folder argument
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


%% Set figure properties
set(groot, 'defaultFigureUnits', 'centimeters', ...
           'defaultFigurePosition', [1, 1, 10, 10], ...
           'defaultFigureColor', 'w', ...
           'defaultFigurePaperUnits', 'centimeters', ...
           'defaultFigurePaperPosition', [0 0 10 10], ...
           'defaultFigurePaperPositionMode', 'auto')
% Set axes properties
set(groot, 'defaultAxesUnits', 'centimeters', ...
           'defaultAxesPosition', [1.5, 1.5, 7.5, 7.5])


%% Create a figure of average bin width histogram
% The different bin widths
X = floor(min(correction.avgBinWidth(:))) : ...
    0.1 : ...
    ceil(max(correction.avgBinWidth(:)));

% Create a new figure
h.fig(2) = figure;
% Calculate a histogram of average bin widths
Y = histc(correction.avgBinWidth(:), X);
% Plot a histogram of bin widths
h.bar = bar(X, Y);
% Set the bar parameters.
set(h.bar, 'BarWidth', 1, ...
           'EdgeColor', 'k', ...
           'LineWidth', 1, ...
           'FaceColor', [0.5 0.5 0.5])
% Label graph
xlabel('Bin Width [ps]')
ylabel('Number of Pixels')
title('TDC Bin Width')


%% Create a figure of TDC bin number per period
% The different bin widths
X = unique(correction.lastBin(:) - correction.firstBin(:));

% Create a new figure
h.fig(3) = figure;
% Calculate a histogram of average bin widths
Y = histc(correction.lastBin(:) - correction.firstBin(:), X);
% Plot a histogram of bin widths
h.bar = bar(X, Y);
% Set the bar parameters.
set(h.bar, 'BarWidth', 1, ...
           'EdgeColor', 'k', ...
           'LineWidth', 1, ...
           'FaceColor', [0.5 0.5 0.5])
% Label graph
xlabel(sprintf('Number of TDC Bins per %d ps', correction.repPeriod))
ylabel('Number of Pixels')
title('Calibrated TDC Bins')


%% Create a figure of mean vs variance plot
% Create a new figure
h.fig(4) = figure;
% Take random samples from the CDM map and plot their mean vs variance.

% Calculate a histogram of average bin widths
% The different bin widths
X = 1 : size(correction.calibratedBins, 3);
Y = squeeze(correction.CDM(1, 1, :));
% Plot a histogram of bin widths
h.line.pulseWidth(1) = plot(X, Y);
hold on
% Draw an extend line
X = double([correction.firstBin(1, 1) * [1, 1, 1], ...
            correction.lastBin(1, 1) * [1, 1, 1]]);
% out half of the median of calibrated bins
Mhalf = double(median(Y(Y > correction.medianPhotons(1, 1) / 2)) / 2);
Y = Mhalf + Mhalf * [-0.1 0.1 0 0 0.1, -0.1];
h.line.pulseWidth(2) = plot(X, Y, 'k', 'LineWidth', 2);
h.txt.pulseWidth = text(mean(X), 0.9 * Mhalf, ...
                        sprintf('\\it{T}\\rm = %g ns', ...
                                0.001 * correction.repPeriod), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'top');

% Label graph
xlabel('Histogram Bin')
ylabel('Number of Photons')
title('Absolute TDC Bin Width Calibration')



%% Save the figures
if nargin < 1
    return
end

% Check if graphics is a character, turn it into a cell
if ischar(graphics)
    graphics = {graphics};
end

% A cell with with graphics formats
formats = {'png', '-dpng'; 'eps', '-depsc'; 'pdf', '-dpdf'};

% Prepare filenames
fnames = {'TDCbinWidth', 'TDCrange', 'TDCcalibration'};
% Go through each figure
for j = 1 : numel(fnames)
    % Go through all requested formats
    for i = 1 : numel(graphics)
        % Check if the requested format is one of the allowed ones
        ftype = find(strcmpi(formats(:, 1), graphics{i}));
        % If it is not found throw an error
        assert(~isempty(ftype), ...
               'Graphics type ''%s'' is not recognized.', ...
               graphics{i})
        % Save the file with the image
        print(h.fig(j + 1), ...
              fullfile(folder, ...
                       sprintf('%s.%s', fnames{j}, formats{ftype, 1})), ...
              formats{ftype, 2})
    end
end
