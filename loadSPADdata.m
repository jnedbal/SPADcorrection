function [SPADdata, inputFiles] = loadSPADdata(SPADfiles)
% LOADSPADDATA is a routine to import TCSPC datasets from SPAD camera.
% These datasets are histograms of photon arrival times produced by 
% the SPAD array. They are often stored in multiple measurements due to the
% length of the experiment and this routine loads all these files and sums
% them into a single histogram
%
% Syntax: [SPADdata, inputFiles] = loadSPADdata(SPADfiles)
%
% Inputs:
%   The function requires a string or a cell of strings as an input. These
%   lead to Matlab MAT-file(s) with the SPAD TCSPC data histograms.
%   SPADfiles can be:
%       (1) Cell array of strings with links to files containing the 
%           code density maps. Each file must contain a variable called
%           XYZimage, which is the SPAD histogram data
%       (2) It can be a regular expression of all file names with repeats
%           of SPAD histogram data.
%
% Outputs:
%   SPADdata:   Summed TCSPC histogramsfrom all the imported files
%   inputFiles: Cell of the source data filenames, including the full path
%
% Examples:
%   loadSPADdata('/path/to/file/XYZimage.mat')
%
% Other m-files required: loadCDMdata, analyzeCDM, fitIRF, correctIRFbyCDM,
%                         analyzeCDMmap, analyzeCDMhist, exGauss, exgfit,
%                         resampleHistogramPar, syntheticPhotons (MEX file)
% Subfunctions: none
%
% See also: loadCDMdata

% Jakub Nedbal
% King's College London
% May 2020
% Last Revision: 
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


%% Load the data
% check if the input is a character array or a cell
% both can be treated in a similar way, except the first few lines of code
% The character is a regular expression converted to a cell array.
assert((ischar(SPADfiles) || iscell(SPADfiles)) && ~isempty(SPADfiles), ...
       'SPADcorrection:loadSPADdata:inputInvalid', ...
       ['The input for analysis must be a cell of filenames or ', ...
        'a character array with regular expression of filenames. ', ...
        'The argument supplied is of class "%s".'], ...
       class(SPADfiles));

% Treat the regular expression in the string
if ischar(SPADfiles)
    % Assume that the input is regular expression leading with links to
    % file(s) containing the CDM data.
    % First separate the directory from the filename
    folder = fileparts(SPADfiles);
    % Load names of all files in that folder
    fileStruct = dir(folder);
    % Convert names to a cell
    inputFiles = {fileStruct.name};
    % Add the folder to the names
    inputFiles = cellfun(@(x) [folder, filesep, x], inputFiles, ...
                        'UniformOutput', false);
else
    % Now treat it as a cell
    inputFiles = SPADfiles;
end

% Keep the list of filenames stored for future reference
% correction.files.CDMfiles = filenames;
% Find indices of all filenames matching the input
% fileIndex = find(~cellfun(@isempty, regexp(inputFiles, SPADfiles)));
fileIndex = find(strcmp(inputFiles, SPADfiles));
% assume that this is a link to the file containing the CDM data
% Check that the file exists, that means the index matrix in not empty
assert(~isempty(fileIndex), 'File %s does not exist.', SPADfiles)
% Create a waitbar
cn = 0;
hwaitbar = waitbar(cn, '', 'Name', 'Loading files...');
% Make it a big higher
hwaitbar.Position(4) = hwaitbar.Position(4) + 25;
tic
% Load each file, one by one
for i = fileIndex
    % Check the variables inside the file
    variableInfo = who('-file', inputFiles{i});
    % Check that each file contains a variable called 'XYZimage'
    if any(contains(variableInfo, 'XYZimage'))
        % Update the waitbar
        if ishandle(hwaitbar)
            % Estimated time to completion
            eta = numel(fileIndex) / cn * toc / 86400;
            if isinf(eta); eta = 0; end
            % Update the waitbar
            [~, fn] = fileparts(inputFiles{i});
            waitbar(cn / numel(fileIndex), ...
                    hwaitbar, ...
                    sprintf('%s\n%d of %d\n%s (%s ETA)', ...
                            fn, ...
                            cn, ...
                            numel(fileIndex), ...
                            datestr(toc / 86400, 'HH:MM:SS'), ...
                            datestr(eta, 'HH:MM:SS')));
            cn = cn + 1;
        end
        % Make a comment
        fprintf('Loading SPAD data file %s...\n', inputFiles{i});
        % Load the input file
        load(inputFiles{i}, 'XYZimage')
        % If this is the first instance of of a file with CDM, just
        % keep the data in the results matrix
        if i == fileIndex(1)
            % Make CDM the code density map
            %correction.CDM = XYZimage;
            SPADdata = XYZimage;
        else
            % Add the code density map from the existing file to the
            % stack
            SPADdata = SPADdata + XYZimage;
        end
    else
        % Make a comment for any selected MAT file that does not containt
        % the variable XYZimage
        fprintf(['Warning: File %s does not contain variable ', ...
                 'XYZimage.\n'], inputFiles{i});
    end
end


% Close the waitbar
delete(hwaitbar)