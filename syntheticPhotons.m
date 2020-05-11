% =========================================================================
%  syntheticPhotons.c - generating and binning synthetic photons
%
% This routine resamples histograms of photon arrival times
% It takes the measured bin widths of actual TDC histograms.
% It simulates randomly distributed photons within these histograms.
% It then creates a new histogram with equally sized bins.
%
% Inputs: inputHist, firstCalBin, lastCalBin, realBinWidth, linearBinWidth
% Output: corrHist
%
% inputHist      - raw TDC histogram (MxN uint32 array),
%                  where M is number of bins, N is number of pixels
% firstCalBin    - 1xN int32 vector of first indices of calibrated bins
% lastCalBin     - 1xN int32 vector of last indices of calibrated bins
% realBinWidth   - matrix of calibrated bin widths (MxN double array)
% linearBinWidth - 1xN double vector of linearized bin widths
% 
% corrHist       - MxN uint32 matrix with resampled TDC histograms
%
% The calling syntax is:
%
%		corrHist = arrayProduct(inputHist, ...
%                              firstCalBin, ...
%                              lastCalBin, ...
%                              realBinWidth, ...
%                              linearBinWidth)
%
%
% Copyright 2020 Jakub Nedbal
% BSD license
%	
% =========================================================================
*/