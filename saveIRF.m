function saveIRF(IRFmap, pixel, fname, binWidth)
% Function to save IRF into an ICS file for use in TRI2.
%   IRFmap      is a x, y, t TDC histogram from the SPAD camera
%   pixel       is the pixel coordinate to save
%   fname       is the filename of the ICS file
%   binWidth    is the width of each time bin in picoseconds


% Check if DIPimage exists
if ~exist('writeics', 'file') && ~exist('dipio_imagewriteics', 'file')
    error('DIP Image is not installed')
end

% The data needs to be restructures
IRF = dip_image(IRFmap(pixel(1), pixel(2), :), 'uint16');
% Get the range of the histogram in seconds
range = 1e-9 * size(IRFmap, 3) * binWidth;
% IRF needs to have the correct order of the coordinates. This means time
% along the first axis, and 1x1 pixel along the two other axes.
% save using older dipio_imagewriteics function
if numel(size(IRF)) == 3
    IRF = permute(IRF, [3, 2, 1]);
end
% Save the IRF
exportICS3(IRF, fname, range, binWidth);