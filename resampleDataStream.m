function corrHist = resampleDataStream(microtimeData)
% RESAMPLEDATASTREAM linearizes SPAD TDC histograms
% Function that processes an existing histogram and resamples it to take
% into account the actual calibrated TDC bin widths. The functions uses a
% Parallel Computing Toolbox, if deemed helpful to speed up the processing
% speed. The heavy number-crunching is done in a compiled C-code (MEX).
%
% Syntax: corrHist = resampleDataStream(dataStream)
%
% Inputs:
%   The function requires input of a 3D array containing a histogram of
%   measured TDC codes. It also requires a global variable 'correction',
%   which carries the calibration data. The function will be modified by
%   this function.
%
%   dataStream is a 3d array of the format XxYxT, where X and Y are SPAD
%           array pixel coordinates and T are the TDC measurements..
%
%
% Outputs:
%   corrHist is a 3D array of TDC histograms resampled according to the
%            calibrated bin widths.
%
%
% Examples:
%   corrHist = resampleHistogramPar(dataStream)
%
% Toolbox requirement: Parallel Processing Toolbox (optional)
% Other m-files required: none
% Subfunctions: syntheticPhotons.c, extractHistogramData
% MAT-files required: none, but global variable 'correction' is needed
%
% See also: analyzeCDM, SPADcorrection

% Jakub Nedbal
% King's College London
% Aug 2018
% Last revision: 15-Apr-2021 - Changed to use global bin width
% Revision: 11-May-2020 - Added support for skew correction
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
extractDataStreamData;

tic

% Run the analysis single-threaded
corrHist = syntheticOJIP(microtimeData, ...
                         firstCalBin, ...
                         lastCalBin, ...
                         realBinWidth, ...
                         linBinWidth, ...
                         peakPos);

% Reshape the matrix to get back to the size of the image array
%corrHist = reshape(corrHist', inputSize);
toc



