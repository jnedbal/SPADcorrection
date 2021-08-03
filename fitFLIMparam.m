function fitFLIMparam

global correction
%% Create input datasets for slimcurve fluorescence lifetime fitting
% This will require the 'prompt' (IRF), a synthetic prompt, 'fit_start', 
% which is the position of the rising edge, 'start', which is the position
% of the start of the fitting. 'fit_end' is the end of the fit. 

% Make a comment
fprintf('Running IRF correction. It may take a while. Stand by...\n')
% First start by shifting the IRF
correction.IRF.corrected = resampleHistogramPar(correction.IRF.raw);

% Find the peak
[~, correction.fitFLIM.peak] = ...
    max(sum(correction.IRF.corrected .* ...
            uint32(correction.IRF.fit.goodfit), [1, 2]));

% Find the starting of the fit, the rising edge
% 02-Jun-2020: I added the factor of three to expand the start of the fit.
% While I learned from Paul Barber that the start should be the position of
% the steepest gradient on the rising edge, in practice, it should be
% probably more like the start of the rise. Adding factor of x3 solved it.
% Now I have better fits at the start of the decay that before, when there
% was no multiplicative factor.
% 15-Apr-2021: I changed from using falling and rising edges to FWHM. 
correction.fitFLIM.start = ...
    round(mean(correction.fitFLIM.peak - ...
               1.5 * correction.IRF.peak.FWHM(:)));

% Find the starting of the fit index
correction.fitFLIM.fit_start = ...
    round(mean(correction.fitFLIM.peak + ...
               0.5 * correction.IRF.peak.FWHM(:)));

% Find the ending edge
% The code below does not work if the rep rate during the CDM measurement
% is not the same as the reprate in the IRF experiment
% The temporary solution is:
correction.fitFLIM.fit_end = ...
    floor(correction.fitFLIM.start + ...
    min(correction.lastBin(:) - ...
        correction.firstBin(:) - ...
        correction.IRF.peak.Pos(:)));
% % To get the ending edge, find the peak with the lowest peak position and
% % subtract it from the width of the CDM
% peakRange = correction.lastBin - ...
%             correction.firstBin - ...
%             correction.IRF.peak.PosInterp;
% % Find the lowest peak position in each row
% [m, inR] = min(peakRange);
% % Find the column index of the lowest peak position
% [~, inC] = min(m);
% % Find the row index of the lowest peak position
% inR = inR(inC);
% % Store the end position
% correction.fitFLIM.fit_end = floor(peakRange(inR, inC));

% Added 02-Jun-2020
% Find the start of consistent data. This is to account for the timing skew
peakRange = correction.IRF.peak.Pos(:);
correction.fitFLIM.data_start = ceil(max(peakRange) - min(peakRange));

% Store the experimental prompt
% Find the index of the prompt
promptStart = round(mean(correction.fitFLIM.peak - ...
                         correction.IRF.peak.FWHM(:)));

% Find the last index of the prompt
promptEnd = ...
        round(mean(correction.fitFLIM.peak + ...
                   1.5 * correction.IRF.peak.FWHM(:)));

% Store the experimental prompt
correction.fitFLIM.expPrompt = ...
    double(correction.IRF.corrected(:, :, promptStart : promptEnd));

% Normalize the prompt
correction.fitFLIM.expPrompt = correction.fitFLIM.expPrompt ./ ...
                               sum(correction.fitFLIM.expPrompt, 3);

% Create model prompt
% first create the time domain
correction.fitFLIM.timeBin = promptStart : promptEnd;
% timeBin = reshape(correction.fitFLIM.timeBin, ...
%                   1, 1, numel(correction.fitFLIM.timeBin));
%
% Calculate the average mu for the exGauss function
%correction.fitFLIM.mu = ...
%    correction.fitFLIM.peak + ...
%    mean(correction.IRF.fit.mu(correction.IRF.fit.goodfit) - ...
%         correction.IRF.peak.PosInterp(correction.IRF.fit.goodfit));

% Simulate the prompt for each pixel
%correction.fitFLIM.simPrompt = exGauss(timeBin, ...
%                                       correction.IRF.fit.h, ...
%                                       correction.fitFLIM.mu, ...
%                                       correction.IRF.fit.sigma, ...
%                                       correction.IRF.fit.tau, ...
%                                       correction.IRF.fit.offset);
% Normalize the prompt
%correction.fitFLIM.simPrompt = correction.fitFLIM.simPrompt ./ ...
%                               sum(correction.fitFLIM.simPrompt, 3);

%% Create an average experimental prompt from well-fit pixels
% First, select only the good pixel data
correction.fitFLIM.avgExpPrompt = ...
    correction.fitFLIM.expPrompt .* correction.IRF.fit.goodfit;
% Sum over the well fitted pixels
correction.fitFLIM.avgExpPrompt = ...
    squeeze(sum(correction.fitFLIM.avgExpPrompt, [1 2], 'omitnan'));

% Normalize the experimental prompt
correction.fitFLIM.avgExpPrompt = ...
    correction.fitFLIM.avgExpPrompt / sum(correction.fitFLIM.avgExpPrompt);

%% Create an average experimental prompt that spans the whole TDC histogram
%  range
correction.fitFLIM.avgExpPromptFull = ...
    double(correction.IRF.corrected) .* correction.IRF.fit.goodfit;

% Sum over the well fitted pixels
correction.fitFLIM.avgExpPromptFull = ...
    squeeze(sum(correction.fitFLIM.avgExpPromptFull, [1 2], 'omitnan'));

% Normalize to the full range of uint16 if larger than that
if max(correction.fitFLIM.avgExpPromptFull) > intmax('uint16')
    correction.fitFLIM.avgExpPromptFull = ...
        correction.fitFLIM.avgExpPromptFull / ...
        max(correction.fitFLIM.avgExpPromptFull) * ...
            double(intmax('uint16'));
end

% %% Create a simulated average prompt
% %  Calculate the average parameters
% h = mean(correction.IRF.fit.h(correction.IRF.fit.goodfit));
% sigma = mean(correction.IRF.fit.sigma(correction.IRF.fit.goodfit));
% tau = mean(correction.IRF.fit.tau(correction.IRF.fit.goodfit));
% offset = mean(correction.IRF.fit.offset(correction.IRF.fit.goodfit));
% 
% % Simulate the average prompt
% correction.fitFLIM.avgSimPrompt = exGauss(timeBin, ...
%                                           h, ...
%                                           correction.fitFLIM.mu, ...
%                                           sigma, ...
%                                           tau, ...
%                                           offset);
% % Normalize the prompt
% correction.fitFLIM.avgSimPrompt = ...
%     squeeze(correction.fitFLIM.avgSimPrompt ./ ...
%             sum(correction.fitFLIM.avgSimPrompt, 3));
% 
% %% Create an average simulated prompt that spans the whole TDC histogram
% %  range
% % Create a long time bin
% timeBin = 1 : size(correction.calibratedBins, 3);
% % Simulate the average prompt
% correction.fitFLIM.avgSimPromptFull = exGauss(timeBin, ...
%                                               h, ...
%                                               correction.fitFLIM.mu, ...
%                                               sigma, ...
%                                               tau, ...
%                                               offset);
% 
% % Normalize to the full range of uint16 if smaller than that
% if max(correction.fitFLIM.avgSimPromptFull) > intmax('uint16')
%     correction.fitFLIM.avgSimPromptFull = ...
%         correction.fitFLIM.avgSimPromptFull / ...
%         max(correction.fitFLIM.avgSimPromptFull) * intmax('uint16');
% end