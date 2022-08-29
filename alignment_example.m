% This is some sample code demonstrating a use case to synchronize a set of
% 24 4D scans programmatically.

% Crops and expected periods are defined individually for each 4D scan,
% then TDCG_sync is run on each. This example does not use Imaris Convert.

clc;
close all;
clear variables;

% These are used as initial guesses for the output period. They are
% automatically refined by the algorithm, so they can be specified
% imprecisely.
expectedPeriods = [
    85 %1
    80 %2
    80 %3
    80 %4
    80 %5
    80 %6
    80 %7
    80 %8
    75 %9
    75 %10
    73 %11
    71 %12
    68 %13
    68 %14
    68 %15
    66 %16
    65 %17
    66 %18
    65 %19
    65 %20
    65 %21
    63 %22
    63 %23
    61 %24
];

% These crops are used only when calculating the TDCG; they are not
% included in the output. Cropping is generally unnecessary, but may
% improve alignment by improving signal to noise ratio of the TDCG. When
% cropping, be careful not to exclude moving structures, even if they are
% outside the heart.
crops = [
    124 170 324 498 %1
    124 170 324 498 %2
    124 170 324 498 %3
    124 170 324 498 %4
    124 170 324 498 %5
    124 170 324 498 %6
    124 170 324 498 %7
    72 37 432 563 %8
    72 37 432 563 %9
    72 37 432 563 %10
    72 37 432 563 %11
    72 37 432 563 %12
    72 37 432 563 %13
    72 37 432 563 %14
    72 37 432 563 %15
    1 1 600 600 %16
    1 1 600 600 %17
    1 1 600 600 %18
    1 1 600 600 %19
    1 1 600 600 %20
    1 1 600 600 %21
    1 1 600 600 %22
    1 1 600 600 %23
    1 1 600 600 %24
];

for i = 1:24
    % The first channel path is the structure, while the remainder can be
    % anything. Typically, the second is Doppler.
    channelPaths = {
        sprintf('F:/E5D_08302016/%02d_structure_20-120.raw', i)
        sprintf('F:/E5D_08302016/%02d_Doppler_40-120.raw', i)
    };

    % The directory to save results to.
    outPath = sprintf('F:/E5D_08302016/aligned/E5D_08302016_align_%02d', i);

    % Scale information shared between the scans. This is only used for
    % Imaris Convert and for certain statistics in the plots; it is not
    % relevant to alignment.
    umPerPixelX = 981/600;   % The X-axis scale for Imaris (microns/px).
    umPerPixelY = 1.95;      % The Y-axis scale for Imaris (microns/px).
    umPerBScan = 1125/20000; % The frame step size in microns/Bscan.

    TDCG_sync_12(channelPaths, ...
        outPath, ...
        expectedPeriods(i), ...
        [umPerPixelX umPerPixelY umPerBScan], ... % [X Y Z] scale.
        [600 600], ... % [height width] B-scan dimensions in the RAW file.
        1:20000, ... % There are 20000 B-scans in the RAW file. Use all.
        'rawStackOffsets', 0, ... % The input RAW files have no headers.
        'verbosity', 2, ... % Show and save diagnostic charts.
        'correlationThreshold', 0.5, ... % 0.4-0.6 is typically good.
        'imarisConvertPath', [], ... % Do not use Imaris Convert.
        'imageCrop', crops(i, :)); % Crop each frame.
end