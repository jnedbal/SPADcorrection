function analyzeIRFmap(graphics, folder)
% analyzeIRFmap displays linearized annd corrected SPAD array instrument 
% response function 
%
% Syntax: analyzeIRFmap(graphics, folder)
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
%   Image files '/path/to/IRFcorrection.xxx' where xxx can be one or more 
%   of png, eps, and pdf. The images show peak position of the instrument
%   response function across the SPAD array and their projection along the
%   vertical and horizontal axes. There are two sets of these images. One
%   is for the linearized SPAD data and the other is for the linearized and
%   timing skew corrected SPAD data.
%
% Outputs:
%   An interactive figure to play with.
%
% Examples:
%   analyzeIRFmap('png')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: SPADcorrection

% Jakub Nedbal
% King's College London
% June 2020
% Last Revision:  03-June-2020 - First working prototype
% Last Revision:  15-Apr-2021  - Added support for global bin width
%
% Copyright 2020-21 Jakub Nedbal
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
global currPeak
global fitRange



%% If no folder is supplied, use the existing one
if nargin < 2
    folder = pwd;
end


%% STart preparing the figure parameters
sSize = get(0, 'ScreenSize');
% Prepare figure size
fSize = [(sSize(3) - sSize(4)) / 2, sSize(4) * 0.1, sSize(4) * [1, 0.8]];
% Prepare axes size
aSize = [fSize(3) * 0.1, fSize(3) * 0.4, ...
         min([0.4 * fSize([4, 4]); 0.8 * fSize([3, 3])])];
aSize = repmat(aSize, 6, 1);
aSize(2, :) = aSize(2, :) .* [1, 0.25, 0.5, 0.5];
aSize(3, :) = aSize(2, :) + [aSize(2, 3), 0, 0, 0];

% Create a figure
h.fig = figure('Units', 'pixels', 'Position', fSize, ...
               'Name', 'displayFLIM', 'NumberTitle', false);
% Center the figure
movegui(h.fig, 'center')

%% Create axes to show the uncorrected peak position
h.ax(1) = axes('Units', 'pixels', 'Position', aSize(1, :));

% Create surface data
rawSurfData = ...
    flipud(correction.IRF.peak.Pos .* correction.avgBinWidth);
% Plot the peak position of the IRF
h.ax1.surf = surf(rawSurfData);
% Remove the obscuring lines
h.ax1.surf.EdgeColor = 'none';
% Set the axis limits to be tight
h.ax(1).XLim = [0, size(correction.IRF.peak.Pos, 2) - 1];
h.ax(1).XTick = [h.ax(1).XTick, h.ax(1).XLim(2)];
h.ax(1).YLim = [0, size(correction.IRF.peak.Pos, 1) - 1];
h.ax(1).YTick = [h.ax(1).YTick, h.ax(1).YLim(2)];
% Add box around the graph
h.ax(1).Box = 'on';

% Label the axes
h.ax1.zlabel = zlabel('IRF Peak Position [ps]');
h.ax1.xlabel = xlabel('Pixel Index', 'Rotation', 20);
h.ax1.xlabel.Position(1) = h.ax(1).XLim(2) * 0.35;
h.ax1.ylabel = ylabel('Pixel Index', 'Rotation', -31);
h.ax1.ylabel.Position(2) = h.ax(1).YLim(2) * 0.35;
h.ax1.title = title('Linearization Correction');

%% Create axes to show the horizontal peak position projection
h.ax(2) = axes('Units', 'pixels', 'Position', aSize(2, :));
% Plot the horizontal peak projection
h.ax2.line(1) = plot(0 : h.ax(1).YLim(2), mean(rawSurfData, 2),  ...
                     'LineWidth', 2);
hold on
spread = [mean(rawSurfData, 2) + std(rawSurfData, [], 2), ...
          mean(rawSurfData, 2) - std(rawSurfData, [], 2)]';
% Plot the errors of the horizontal peak position
h.ax2.line(2) = plot(0 : h.ax(1).YLim(2), spread(1, :),  ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax2.line(3) = plot(0 : h.ax(1).YLim(2), spread(2, :), ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax2.area = fill([0 : h.ax(1).YLim(2), h.ax(1).YLim(2) : -1 : 0], ...
                  [spread(1, :), fliplr(spread(2, :))], [0.7, 0.7, 0.7]);
uistack(h.ax2.area, 'bottom')
h.ax(2).XDir = 'reverse';
h.ax(2).XLim = h.ax(1).YLim;
h.ax(2).XTick = h.ax(1).YTick;
h.ax(2).YLim = h.ax(1).ZLim;
h.ax(2).YTick = h.ax(1).ZTick;
h.ax(2).YGrid = 'on';
ylabel(h.ax1.zlabel.String);
title({'Vertical', 'Projection'})
xlabel('Pixel Index')

%% Create axes to show the vertical peak position projection
h.ax(3) = axes('Units', 'pixels', 'Position', aSize(3, :));
% Plot the horizontal peak projection
h.ax3.line(1) = plot(0 : h.ax(1).XLim(2), mean(rawSurfData, 1),  ...
                     'LineWidth', 2);
hold on
spread = [mean(rawSurfData, 1) + std(rawSurfData, [], 1); ...
          mean(rawSurfData, 1) - std(rawSurfData, [], 1)];
% Plot the errors of the horizontal peak position
h.ax3.line(2) = plot(0 : h.ax(1).XLim(2), spread(1, :),  ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax3.line(3) = plot(0 : h.ax(1).XLim(2), spread(2, :), ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax3.area = fill([0 : h.ax(1).XLim(2), h.ax(1).XLim(2) : -1 : 0], ...
                  [spread(1, :), fliplr(spread(2, :))], [0.7, 0.7, 0.7]);
uistack(h.ax3.area, 'bottom')
h.ax(3).XLim = h.ax(1).XLim;
h.ax(3).XTick = h.ax(1).XTick;
h.ax(3).YLim = h.ax(1).ZLim;
h.ax(3).YTick = h.ax(2).YTick;
h.ax(3).YTickLabel = cell(1, numel(h.ax(3).YTick));
h.ax(3).YGrid = 'on';
title({'Horizontal', 'Projection'})
xlabel('Pixel Index')
h.leg(1) = legend(h.ax3.line([1, 2]), {'Average', 'StdDev'});

% Stretch and orient the legend
h.leg(1).Units = 'pixel';
h.leg(1).Orientation = 'horizontal';
h.leg(1).Position(1) = aSize(2, 1) + aSize(2, 3) - h.leg(1).Position(3) / 2;
h.leg(1).Position(2) = aSize(2, 2) + 10;

%% Now correct the IRF and find the peak positions
% Create a waitbar
hwaitbar = waitbar(0, '', ...
                   'Name', 'Finding IRF parameters...');
% Make it a big higher
hwaitbar.Position(4) = hwaitbar.Position(4) + 10;

% Run the linearization routine
corrIRF = resampleHistogramPar(correction.IRF.raw);
corrIRF = reshape(corrIRF, ...
                  size(corrIRF, 1) * size(corrIRF, 2), size(corrIRF, 3))';
corrIRF = double(corrIRF);
% Create a matrix of of peak positions
peakPos = correction.IRF.peak.Pos;
% Create a matrix of bin indices
%pixIndex = (1 : size(corrIRF, 2));
% Fit parameters
%fitParam = zeros(1, 5);
% Guess the first peak position
%[~, fitParam(2)] = max(corrIRF(find(correction.IRF.fit.goodfit, 1), :));
% total number of pixels
numberPixels = size(corrIRF, 2);

%% Smoothen the data to find the peaks
IRFsmooth = smoothdata(corrIRF, 1, ...
                       correction.IRF.fit.param.filter.model, ...
                       correction.IRF.fit.param.filter.kernel);


%% Find the peaks
[peak, index] = max(IRFsmooth);
% Setup a zoom of relevant fit data
fitRange = (-49 : 49)';
% Fixed fit parameters
B = 0;      % Peak Position

% % Set options for the minimization function
% for i = find(correction.IRF.fit.goodfit)'
%     if mod(i, 100) == 0
%         if ishandle(hwaitbar)
%             % Estimated time to completion
%             eta = numberPixels / i * toc / 86400;
%             if isinf(eta); eta = 0; end
%             % Update the waitbar
%             waitbar(i / numberPixels, ...
%                     hwaitbar, ...
%                     sprintf('%d of %d\n%s (%s ETA)', ...
%                             i, ...
%                             numberPixels, ...
%                             datestr(toc / 86400, 'HH:MM:SS'), ...
%                             datestr(eta, 'HH:MM:SS')));
%         end
%     end
%     % Fit the linearized IRF
%     
%     pixHist = corrIRF(i, :);
%     % These are the initial guesses for the fit
%     param0 = [correction.IRF.fit.h(i), ...
%               fitParam(2), ...
%               correction.IRF.fit.sigma(i), ...
%               correction.IRF.fit.tau(i), ...
%               correction.IRF.fit.offset(i)];
%     % Run the fit
%     [fitParam(1), ...
%         fitParam(2), ...
%         fitParam(3), ...
%         fitParam(4), ...
%         fitParam(5)] = exgfit(param0);
%     % Find the position of the peak, and the the maximum of the peak
%     % Minimization function
%     mFun = @(x) (-(exGauss(x, ...
%                            fitParam(1), ...
%                            fitParam(2), ...
%                            fitParam(3), ...
%                            fitParam(4), ...
%                            fitParam(5))));
%     % Run minimization routine to find the peak position
%     peakPos(i) = fminsearch(mFun, fitParam(2));
% end

options = optimset('Display', 'off');
for i = find(correction.IRF.fit.goodfit)'
    % update the waitbar
    if mod(i, 100) == 0
        if ishandle(hwaitbar)
            % Estimated time to completion
            eta = numberPixels / i * toc / 86400;
            if isinf(eta); eta = 0; end
            % Update the waitbar
            waitbar(i / numberPixels, ...
                    hwaitbar, ...
                    sprintf('%d of %d\n%s (%s ETA)', ...
                            i, ...
                            numberPixels, ...
                            datestr(toc / 86400, 'HH:MM:SS'), ...
                            datestr(eta, 'HH:MM:SS')));
        end
    end
    binIndex = index(i) + fitRange;
    currPeak = IRFsmooth(binIndex, i);

    C = correction.IRF.fit.sigma(i) * sqrt(2);  % Sigma * sqrt(2)
    D = correction.IRF.fit.offset(i);           % Offset
    A = peak(i) - D;                            % Amplitude
    % Perform the Gaussian fit
    %[f1, gof] = fit(fitRange, currPeak, gaussEqn, 'Start', [A, B, C, D]);
    fp1 = fminsearch(@fitGauss, [A, B, C, D], options);
    % Deal the results
    peakPos(i) = index(i) + fp1(2);
end

%% Run interpolations for the data that was not properly fitted
% interpolate the FWHM for points that are not fitted well
% Create a matrix of points
[X, Y] = meshgrid(1 : size(peakPos, 2), 1 : size(peakPos, 1));
Xv = X(~correction.IRF.fit.goodfit);
Yv = Y(~correction.IRF.fit.goodfit);
X = X(correction.IRF.fit.goodfit);
Y = Y(correction.IRF.fit.goodfit);

% Interpolate the peak position for points that are not fitted well
% Create an interpolant object that can be called
% Use the natural neighbor interpolation method and 
% nearest neightbor extrapolation method
% to estimate the peak position for the missing pixels
F = scatteredInterpolant(X, Y, peakPos(correction.IRF.fit.goodfit), ...
                         'natural', 'nearest');

% Create a new matrix, which contains the estimated and the interpolated
% values of the full-width half-maximum
peakPosInterp = peakPos;
peakPosInterp(~correction.IRF.fit.goodfit) = F(Xv, Yv);


% Close the waitbar
delete(hwaitbar)


%% Create axes to show the uncorrected peak position
aSize(4, 1) = fSize(3) * 0.6;
aSize(5, :) = aSize(2, :);
aSize(5, 1) = aSize(4, 1);
aSize(6, :) = aSize(5, :) + [aSize(5, 3), 0, 0, 0];

h.ax(4) = axes('Units', 'pixels', 'Position', aSize(4, :));

% Create surface data
corrSurfData = flipud(peakPosInterp .* correction.globalBinWidth.corr);
% Plot the peak position of the IRF
h.ax4.surf = surf(corrSurfData);
% Remove the obscuring lines
h.ax4.surf.EdgeColor = 'none';
% Set the axis limits to be tight
h.ax(4).XLim = [0, size(peakPosInterp, 2) - 1];
h.ax(4).XTick = [h.ax(4).XTick, h.ax(4).XLim(2)];
h.ax(4).YLim = [0, size(peakPosInterp, 1) - 1];
h.ax(4).YTick = [h.ax(4).YTick, h.ax(4).YLim(2)];
% Add box around the graph
h.ax(4).Box = 'on';

% Label the axes
h.ax4.zlabel = zlabel('IRF Peak Position [ps]');
h.ax4.xlabel = xlabel('Pixel Index', 'Rotation', 20);
h.ax4.xlabel.Position(1) = h.ax(4).XLim(2) * 0.35;
h.ax4.ylabel = ylabel('Pixel Index', 'Rotation', -31);
h.ax4.ylabel.Position(2) = h.ax(4).YLim(2) * 0.35;
h.ax4.title = title('Timing Skew Correction');

%% Create axes to show the horizontal peak position projection
h.ax(5) = axes('Units', 'pixels', 'Position', aSize(5, :));
% Plot the horizontal peak projection
h.ax5.line(1) = plot(0 : h.ax(4).YLim(2), mean(corrSurfData, 2),  ...
                     'LineWidth', 2);
hold on
spread = [mean(corrSurfData, 2) + std(corrSurfData, [], 2), ...
          mean(corrSurfData, 2) - std(corrSurfData, [], 2)]';
% Plot the errors of the horizontal peak position
h.ax5.line(2) = plot(0 : h.ax(4).YLim(2), spread(1, :),  ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax5.line(3) = plot(0 : h.ax(4).YLim(2), spread(2, :), ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax5.area = fill([0 : h.ax(4).YLim(2), h.ax(4).YLim(2) : -1 : 0], ...
                  [spread(1, :), fliplr(spread(2, :))], [0.7, 0.7, 0.7]);
uistack(h.ax5.area, 'bottom')
h.ax(5).XDir = 'reverse';
h.ax(5).XLim = h.ax(4).YLim;
h.ax(5).XTick = h.ax(4).YTick;
h.ax(5).YLim = h.ax(4).ZLim;
h.ax(5).YTick = h.ax(4).ZTick;
ylabel(h.ax4.zlabel.String);
h.ax(5).YGrid = 'on';
title({'Vertical', 'Projection'})
xlabel('Pixel Index')

%% Create axes to show the vertical peak position projection
h.ax(6) = axes('Units', 'pixels', 'Position', aSize(6, :));
% Plot the horizontal peak projection
h.ax6.line(1) = plot(0 : h.ax(4).XLim(2), mean(corrSurfData, 1),  ...
                     'LineWidth', 2);
hold on
spread = [mean(corrSurfData, 1) + std(corrSurfData, [], 1); ...
          mean(corrSurfData, 1) - std(corrSurfData, [], 1)];
% Plot the errors of the horizontal peak position
h.ax6.line(2) = plot(0 : h.ax(4).XLim(2), spread(1, :),  ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax6.line(3) = plot(0 : h.ax(4).XLim(2), spread(2, :), ...
                     'Color', [0.3, 0.3, 0.3]);
h.ax6.area = fill([0 : h.ax(4).XLim(2), h.ax(4).XLim(2) : -1 : 0], ...
                  [spread(1, :), fliplr(spread(2, :))], [0.7, 0.7, 0.7]);
uistack(h.ax6.area, 'bottom')
h.ax(6).XLim = h.ax(4).XLim;
h.ax(6).XTick = h.ax(4).XTick;
h.ax(6).YLim = h.ax(4).ZLim;
h.ax(6).YTick = h.ax(5).YTick;
h.ax(6).YTickLabel = cell(1, numel(h.ax(6).YTick));
h.ax(6).YGrid = 'on';
title({'Horizontal', 'Projection'})
xlabel('Pixel Index')
h.leg(2) = legend(h.ax6.line([1, 2]), {'Average', 'StdDev'});

% Stretch and orient the legend
h.leg(2).Units = 'pixel';
h.leg(2).Orientation = 'horizontal';
h.leg(2).Position(1) = aSize(5, 1) + aSize(5, 3) - h.leg(2).Position(3) / 2;
h.leg(2).Position(2) = aSize(5, 2) + 10;



%% Store the corrected peak position in the correction matrix
correction.IRF.peak.corrPosInterp = peakPosInterp;

%% Store the image
% Check if graphics is a character, turn it into a cell
if ischar(graphics)
    graphics = {graphics};
end

% A cell with with graphics formats
formats = {'png', '-dpng'; 'eps', '-depsc'; 'pdf', '-dpdf'};

% Make the font sizes the same
set(h.ax([2 3 5 6]), 'FontSize', h.ax(1).FontSize);
% Get rid of the overlapping zeros in the two adjacent axes
h.ax(3).XTickLabel{1} = '';
h.ax(6).XTickLabel{1} = '';
% Get the page size correct for proper PDF export
h.fig.Units = 'centimeters';
h.fig.PaperUnits = h.fig.Units;
h.fig.PaperPosition = [0, 0, h.fig.Position([3, 4])];
h.fig.PaperSize = h.fig.Position([3, 4]);
h.fig.PaperPositionMode = 'auto';

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
                   sprintf('IRFcorrection.%s', formats{ftype, 1})),...
          formats{ftype, 2})
end

saveas(h.fig, fullfile(folder, 'IRFcorrection.fig'))
end

function z = fitGauss(p)
global currPeak
global fitRange
%cx = p(1);
%wx = p(2);
%amp = p(3);

zx = p(1) * exp(-(fitRange - p(2)) .^ 2 / (p(3) ^ 2)) + p(4) - currPeak;

z = sum(zx.^2);
end