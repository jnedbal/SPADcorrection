% remember the size of the input array
inputSize = size(XYZimage);

% Number of bins in the data
numberBins = int32(inputSize(3));
% Number of pixels in the data
numberPixels = int32(inputSize(1) * inputSize(2));


% Squish linearized bin width into a 2D array for easier indexing
realBinWidth = reshape(correction.binWidth, numberPixels, numberBins)';

% Squish the input histogram into a 2D array for easier indexing
XYZimage = reshape(XYZimage, numberPixels, numberBins)';

% Extract idealized bin width
linBinWidth = reshape(correction.idealBinWidth, numberPixels, numberBins);
linBinWidth = max(linBinWidth, [], 2);

% Find the first calibrated bin for each pixel
[~, firstCalBin] = max(correction.calibratedBins, [], 3);
% Organize into a column vector. Subtract one, because MEX-file indexing
% uses 0, as the least index, not 1.
firstCalBin = int32(firstCalBin(:)) - 1;
% Number of calibrated bins
nrBins = sum(correction.calibratedBins, 3);
nrBins = int32(nrBins(:));
% Since the last bin is not ending at the same point as the
% idealized bin, we need to simulate one extra bin, which is going
% to be discarded at the end
lastCalBin = firstCalBin + nrBins;
