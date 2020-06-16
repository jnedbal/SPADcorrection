function corrHist = resampleHistogramPar(XYZimage)
% RESAMPLEHISTOGRAMPAR linearizes SPAD TDC histograms
% Function that processes an existing histogram and resamples it to take
% into account the actual calibrated TDC bin widths. The functions uses a
% Parallel Computing Toolbox, if deemed helpful to speed up the processing
% speed. The heavy number-crunching is done in a compiled C-code (MEX).
%
% Syntax: corrHist = resampleHistogramPar(XYZimage)
%
% Inputs:
%   The function requires input of a 3D array containing a histogram of
%   measured TDC codes. It also requires a global variable 'correction',
%   which carries the calibration data. The function will be modified by
%   this function.
%
%   inputHist is a 3d array of the format XxYxT, where X and Y are SPAD
%           array pixel coordinates and T are the TDC histogram values. It
%           must be same size as the matrices in 'correction'.
%
%
% Outputs:
%   corrHist is a 3D array of TDC histograms resampled according to the
%            calibrated bin widths.
%
%
% Examples:
%   corrHist = resampleHistogramPar(inputHist)
%
% Other m-files required: none
% Subfunctions: syntheticPhotons.c, extractHistogramData
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: analyzeCDM, SPADcorrection

% Jakub Nedbal
% King's College London
% Aug 2018
% Last revision: 11-May-2020 - Added support for skew correction
% Revision: 11-May-2020 - Tidied up the code
% Revision: 14-Apr-2020
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
global correction %#ok<NUSED>


%% Run a routine that extract all important aspects of the linearization
%  and the input histogram data
extractHistogramData;

%% Check if parallel processing toolbox is installed and if it makes sense 
%  using it
% First check if parallel processing toolbox is available
isParCompAvailable = ~isempty(ver('distcomp'));
% Set the minimum reasonable number of photons to bother with parallel
% computing (1e9 seems good on a modern PC)
minPhotons = 1e9;
% It makes sense to run parallel computing toolbox if there are more than
% 1B photons
isParCompHelpful = sum(XYZimage(:)) > minPhotons;

tic

% Run it in parallel if it is helpful and toolbox is available
if isParCompAvailable && isParCompHelpful
    % parfor loop has got a peculiar behavior, where it cannot work with
    % variables that have been loaded from a MAT file. The solution is to
    % explicitly make them copies of themselves. Weird, but it works
    firstCalBin = firstCalBin; %#ok<*ASGSL,*NODEF>
    lastCalBin = lastCalBin;
    realBinWidth = realBinWidth;
    linBinWidth = linBinWidth;
    peakPos = peakPos;
    % Divide the whole dataset into a number chunks equal to the maximum
    % number of threads available on the machine
    chunks = int32(maxNumCompThreads);
    % Each chunk is a block of a number of pixels
    block = numberPixels / chunks;
    % Preallocate a cell array of integer matrices
    % parfor does not allow updating the same array, hence the cell array
    cH = repmat({zeros(numberBins, numberPixels / chunks, 'int32')}, ...
                1, chunks);
    % Run a resampling code in four chunks to benefit from multicore CPU
    parfor i = 0 : (chunks - 1)
        % index of the block of input data
        index = (1 : block) + i * block;
        % Run the compiled C-code for the bin resampling
        cH{i + 1} = syntheticPhotons(XYZimage(:, index), ...
                                     firstCalBin(index), ...
                                     lastCalBin(index), ...
                                     realBinWidth(:, index), ...
                                     linBinWidth(index), ...
                                     peakPos(index)); ...
                                     %#ok<PFBNS>
    end
    % Concatenate the results into a same matrix
    corrHist = horzcat(cH{:});
else
    % Run the analysis single-threaded
    corrHist = syntheticPhotons(XYZimage, ...
                                firstCalBin, ...
                                lastCalBin, ...
                                realBinWidth, ...
                                linBinWidth, ...
                                peakPos);
end
% Reshape the matrix to get back to the size of the image array
corrHist = reshape(corrHist', inputSize);
toc



