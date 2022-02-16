function exportICS3(utimeHist, fname, range, binWidth)
% EXPORTICS2 is a code for saving microtime histograms into ICS image
% format for use with the TRI2 software
% It requires having the dip_image library installed and loaded.
%
% Inputs:
%   utimeHist: Microtime histogram in the format YxXxT
%   fname:     Filename of the file into which the data should be saved
%   range:     The time range of the histogram in seconds, 
%              i.e. number of bins * bin width
%   binWidth:  duration of each bin in nanoseconds
%
%
% Examples:
%   correctCDMbyIRF
%
% Other m-files required: dipio_imagewriteics, dip_image, writeics
% Subfunctions: none
% MAT-files required: none
%
% See also: dipio_imagewriteics, dip_image, writeics
%
%% Jakub Nedbal
% King's College London
%
% Last Revision: 04-Feb-2022 - added support for writeics in newer DIPimage
%                              versions
%                14-May-2020
%
% Copyright 2020-22 Jakub Nedbal
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


% Range is the number of bins multiplied by seconds per bin
if nargin < 3
    % Make some guess if not available
    range = 512 * 40e-12;
end
% Bin width is the size of each TDC bin in the histogram
if ~exist('binWidth', 'var')
    % Bin with in nanoseconds
    binWidth = range * 1e9 / (size(utimeHist, 1) - 1);
end

% Check if DIPimage exists
if ~exist('writeics', 'file') && ~exist('dipio_imagewriteics', 'file')
    error('DIP Image is not installed')
end


%% Generate histograms in ICS format
physDims.dimensions = [binWidth,  0.296296, 0.296296];
physDims.origin = [0 0 0];
physDims.dimensionUnits = {'ns', 'microns', 'microns'};
physDims.intensity = 1;
physDims.offset = 0;
physDims.intensityUnit = 'relative';
history = ...
    {'microscope name	Caerdydd', ...
     'metadata format ver	1.0', ...
     'microscope	MicroFLiC 1.0', ...
     'manufacturer	King''s College London', ...
     'microscope built on	$WCNOW$', ...
     'atd_microscopy ver	$WCREV$ $WCDATE$', ...
     'atd_hardware ver	$HARDWARE_REV$ $HARDWARE_DATE$', ...
     'atd_libraries ver	$LIBRARIES_REV$ $LIBRARIES_DATE$', ...
     'icsviewer ver	$ICSVIEWER_REV$ $ICSVIEWER_DATE$', ...
     'experimenter	Jakub Nedbal', ...
     'study	DefaultStudy', ...
     'experiment	DefaultExperiment', ...
     sprintf('creation date	%s', datestr(now, 'HH:MM:SS dd\\mm\\yy')), ...
     'type	FluorescenceLifetime', ...
     'labels	t x y', ...
     sprintf('dimensions	%d %d %d', size(utimeHist)), ...
     'offsets	0 0 0', ...
     'units	s m m', ...
     sprintf('extents	%g %g %g', ...
             range, ...
             size(utimeHist, 1), ...
             size(utimeHist, 3))};


% save using older dipio_imagewriteics function
if ~exist('writeics', 'file')
    % reshuffle the data - this is needed for older dipio_imagewriteics
    % function
    if numel(size(utimeHist)) == 3
        utimeHist = permute(utimeHist, [1, 3, 2]);
    end

    % Export the file
    dipio_imagewriteics(dip_image(utimeHist, 'uint16'), ... % Data
                        fname, ...                     % Filename
                        'gray', ...                    % Photometrix
                        physDims, ...                  % Physical dimension
                        history, ...                   % Tags
                        16, ...                        % Significant bits
                        2, ...                         % ICS version
                        1)                             % Compression
else

    % reshuffle the data - this is needed for older dipio_imagewriteics
    % function
    if numel(size(utimeHist)) == 3
        utimeHist = permute(utimeHist, [3, 1, 2]);
        utimeHist = dip_image(utimeHist, 'uint16');
    else
        % Expand the dimensionality to match the requirements of TRI2
        uh = newim(1, numel(utimeHist), 1, 'uint16');
        uh(:, :, :) = dip_image(utimeHist, 'uint16');
        utimeHist = uh;
    end
    writeics(utimeHist, ...     % Data
             fname, ...         % Filename
             history, ...       % Tags
             16, ...            % Significant bits
             {'v2', ...         % ICS version
              'gzip', ...       % Compression
              'fast'});         % write directly
end