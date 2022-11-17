% This sample code demonstrates synchronization of a single channel
% dataset. The RAW file can be downloaded from
% https://doi.org/10.6084/m9.figshare.21563472
% ... and placed in the current working directory.

% Only one channel is specified.
channelPaths = {
    'E8p5_072416_5D_600X30000_0.8VX0.7V_12ms_SP_T20_20-120_600_DCrem.raw'
};

outPath = 'aligned_TDCG';

TDCG_sync(channelPaths, ...
    outPath, ...
    49, ... % The expected period is known.
    [726/600 1.95 875/30000], ... % [X Y Z] scale.
    [600 600], ... % [height width] B-scan dimensions in the RAW file.
    1:30000, ... % Use all 30000 B-scans in the RAW file.
    'rawStackOffsets', 35, ... % The RAW file has a 35 byte header.
    ... % If Imaris is installed, add the ImarisConvert.exe path here.
    'imarisConvertPath', [], ... % Otherwise, no Imaris file will be made.
    'verbosity', 2); % Show and save diagnostic charts.