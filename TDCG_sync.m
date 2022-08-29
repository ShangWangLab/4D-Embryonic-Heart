%% Example:
% % The list of paths to stacks containing associated channel information.
% % Add as many strings to this list as you want.
% % The first channel must always be the structure channel.
% channelPaths = {
%     ['F:/syncEmbryonicHeart2_TDCG/' ...
%     'E5D_24_Structure_600X20000_0.9VX0.9V_20-120_750_debanded.raw']
% };
% 
% % The directory to save output files to.
% outPath = 'F:/syncEmbryonicHeart2_TDCG/';
% 
% % Image parameters.
% rawStackOffsets = 35; % Preamble length before image data begins (units of bytes).
% imageShape = [750 600]; % [height width] (units of pixels).
% imageIndices = 10559 : 19575; % Frame indices to use (1-based).
% 
% expectedHeartPeriod = 60; % Frames per beat (within 10% of the true period).
% umPerPixelX = 978/600; % The X-axis scale for Imaris (microns/px).
% umPerPixelY = 1462/750; % The Y-axis scale for Imaris (microns/px).
% umPerBScan = 0.05625; % The frame step size in microns/Bscan.
% 
% % Specifies the region of the structure to use when calculating the TDCG.
% % Usually, no cropping is necessary, although accuracy may be improved by
% % cropping smaller hearts to highlight the relevant region. This does not
% % effect the size of the aligned output.
% % Standard MATLAB one-based indexing.
% % Format is [row1, column1, row2, column2] in units of pixels.
% imageCrop = [65 40 367 571];
% 
% % Determines the level of detail at while figures will be generated.
% % Level 0: Silent Mode (No logging or figures).
% % Level 1: Regular Mode (Usual logging and important figures).
% % Lever 2: Diagnostic Mode (Show statistics).
% % Level 3: Debug Mode (Show interactive figures).
% verbosity = 2;
% 
% % The path to your installation of ImarisConvert.exe.
% % Set this to [] if you want to disable Imaris file creation.
% imarisConvertPath = 'C:/Program Files/Bitplane/Imaris x64 9.6.0/ImarisConvert.exe';
% 
% % The normalized correlation threshold, above which, a match is considered
% % to be good. Good values are typically between 0.4-0.6.
% correlationThreshold = 0.6;
% 
% % The n = 2k+1 number of time samples to use when averaging a reference frame.
% k = 2;
%
% Finally, run the main function with all the parameters defined:
% TDCG_sync(channelPaths, outPath, expectedHeartPeriod, ...
%     [umPerPixelX umPerPixelY umPerBScan], imageShape, imageIndices
%     'imageCrop', imageCrop, 'rawStackOffsets', rawStackOffsets, ...
%     'verbosity', verbosity, 'imarisConvertPath', imarisConvertPath, ...
%     'correlationThreshold', correlationThreshold, 'k', k ...
% );

% Comment to suppress global variable warnings in the MATLAB editor:
%#ok<*GVMIS>

%% Functions
function TDCG_sync(varargin)

    % You may want to edit this line if using a different version of
    % Imaris.
    defaultImarisConvertPath = ...
        'C:/Program Files/Bitplane/Imaris x64 9.6.0/ImarisConvert.exe';

    %% Parse inputs
    parser = inputParser();
    addRequired(parser, 'channelPaths', @(x) all(cellfun(@ischar, x)));
    addRequired(parser, 'outPath', @(x) exist(fileparts(x), 'dir'));
    addRequired(parser, 'expectedHeartPeriod', @(x) mod(x, 1) == 0);
    addRequired(parser, 'inputResolution', @isnumeric);
    addRequired(parser, 'imageShape', @(x) all(mod(x, 1) == 0));
    addRequired(parser, 'imageIndices', @(x) all(mod(x, 1) == 0));
    addParameter(parser, 'rawStackOffsets', 0, @(x) all(mod(x, 1) == 0 & x >= 0));
    addParameter(parser, 'imageCrop', [], @(x) isempty(x) || numel(x) == 4 && all(mod(x, 1) == 0));
    addParameter(parser, 'verbosity', 2, @(x) mod(x, 1) == 0 && 0 <= x && x <= 3);
    addParameter(parser, 'imarisConvertPath', defaultImarisConvertPath, @(x) exist(x, 'file') || isempty(x));
    addParameter(parser, 'correlationThreshold', 0.6, @(x) -1 <= x && x <= 1);
    addParameter(parser, 'k', 2, @(x) mod(x, 1) == 0 && 1 <= x && x <= 10);
    parse(parser, varargin{:});
    
    global verbosity;
    
    channelPaths = parser.Results.channelPaths;
    outPath = parser.Results.outPath;
    expectedHeartPeriod = parser.Results.expectedHeartPeriod;
    inputResolution = parser.Results.inputResolution;
    imageShape = parser.Results.imageShape;
    imageIndices = parser.Results.imageIndices;
    rawStackOffsets = parser.Results.rawStackOffsets;
    imageCrop = parser.Results.imageCrop;
    verbosity = parser.Results.verbosity;
    imarisConvertPath = parser.Results.imarisConvertPath;
    correlationThreshold = parser.Results.correlationThreshold;
    k = parser.Results.k;

    % Stack offsets can be listed as a single number for all channels, or
    % as a unique offset for each channel.
    if numel(rawStackOffsets) == 1
        rawStackOffsets = repmat(rawStackOffsets, 1, numel(channelPaths));
    elseif numel(rawStackOffsets) ~= numel(channelPaths)
        error('The number of RAW stack offsets must match the number of channels provided.');
    end

    % Paths are assumed to use forward slashes.
    channelPaths = cellfun(@(x) strrep(x, '\', '/'), channelPaths, 'UniformOutput', false);
    outPath = strrep(outPath, '\', '/');
    if ~isempty(imarisConvertPath)
        imarisConvertPath = strrep(imarisConvertPath, '\', '/');
    end
    
    % Ensure no trailing slash on the directory.
    while ~isempty(outPath) && outPath(end) == '/'
        outPath = outPath(1:end-1);
    end
    
    if isempty(imageCrop)
        imageCrop = [1 1 imageShape];
    end

    % Whether or not to save generated figures to a directory.
    saveFigures = verbosity >= 1;

    % Whether or not to log text written to the console.
    enableConsoleLogging = verbosity >= 1;

    % The number of time samples to use when averaging a reference frame.
    nkPhases = 2*k + 1;
    
    umPerPixelX = inputResolution(1);
    umPerPixelY = inputResolution(2);
    umPerBScan = inputResolution(3);

    if saveFigures
        close all hidden;
    end

    if enableConsoleLogging
        diary([outPath '_log.txt']);
        diaryCloser = onCleanup(@() diary('off'));
    end

    if verbosity >= 2
        fprintf('Average with a sliding window of %d sequences.\n', nkPhases);
    end

    %% Load images into a stack and generate the TDCG.
    if verbosity >= 1
        % The number of microns of depth that the scan covers.
        fprintf('Loading a stack of %d images spanning ~%0.3f μm in Z...\n', ...
            numel(imageIndices), numel(imageIndices) * umPerBScan);
    end

    tic;

    % The RAW image stack showing structure is always the first channel path.
    TDCG = loadTDCG(channelPaths{1}, imageShape, imageCrop, imageIndices, rawStackOffsets(1));

    if verbosity >= 1
        fprintf('Loading the structure stack and computing the TDCG took %.3f seconds.\n', toc);
    end

    tic;

    %% Estimate the dominant period to use for the output sequence.
    if verbosity >= 1
        fprintf('Estimating the heart rate...\n');
    end

    % OP is the output period in units of number of B-scans.
    % jPeak is the starting index of the reference sequence in the TDCG.
    [OP, jPeak] = determineOP(TDCG, expectedHeartPeriod);

    %% Generate the initial running average reference sequence.
    if verbosity >= 1
        fprintf('Aligning the TDCG...\n');
    end

    % The valid distance within to check for alignment of a sequence.
    % Add 1 to allow buffer space for the 3-tap LP filter's time-lag.
    searchRadius = ceil(0.1 * OP) + 1;

    if verbosity >= 2
        fprintf('Searching within +/- %d time samples for correlation matches.\n', ...
            searchRadius - 1);
    end

    % The reference sequence to align kPhases to.
    ref = subSpan(TDCG, jPeak, OP);

    kPhases = zeros(nkPhases, 1);
    kPhases(k+1) = jPeak;

    % Find both the previous and next 'k' alignments based on the reference.
    for i = [-k:-1, 1:k]
        j = jPeak + i*OP;

        j1 = j - searchRadius;
        j2 = j + searchRadius + OP - 1;
        lag = findMaxCorrLP(ref, TDCG(j1:j2));
        kPhases(1+k+i) = j1 + lag;
    end

    % Redefine the reference from here on out to be the running average of n
    % aligned sequences.
    nRefs = zeros(OP, nkPhases);
    for i = 1:nkPhases
        nRefs(:, i) = normalize(subSpan(TDCG, kPhases(i), OP));
    end
    avgRef = mean(nRefs, 2);

    %% Allocate the array containing the image indices for each output sequence.
    % Axes are: (heart phase, Z-coordinate).
    % Allocate enough spots to hold the theoretical maximum output size.
    % Crop off unused slots when done populating the array.
    % Provide an extra 200 slots just in case they are needed. It won't take
    % much extra memory, and it provides safety.
    outPhases = zeros(200 + ceil(numel(TDCG)/(OP - searchRadius)), 1);
    outCorrs = zeros(size(outPhases));
    iCenter = 100 + round(numel(outPhases) * jPeak/numel(TDCG));

    % Populate the output with the first found alignments.
    outPhases(iCenter - k : iCenter + k) = kPhases;
    outCorrs(iCenter - k : iCenter + k) = NaN;

    %% Align sequences.
    if verbosity >= 3
        figure(5);
    end

    % TODO: generalize the forward and backward versions of the following two
    % loops into a single function call.

    % Align center-following signal periods while updating the running reference.
    i = iCenter + k;
    j = outPhases(i);
    refSeqs = nRefs; % A circular buffer of reference sequences to average.
    oldestRef = 1; % A pointer to the first reference to remove.
    while j + 2*OP - 1 <= numel(TDCG)
        i = i + 1;
        j = j + OP;

        j1 = j - searchRadius;
        j2 = min(j + searchRadius + OP - 1, numel(TDCG));
        weightedRef = (avgRef + mean(refSeqs, 2))/2;
        [lag, maxCorr] = findMaxCorrLP(weightedRef, TDCG(j1:j2));

        % When a correlation is poor, do not attempt to align or include the
        % sequence as a reference.
        if maxCorr > correlationThreshold
            j = j1 + lag; % Update j to point to the max correlation.

            % Update the running reference sequence.
            refSeqs(:, oldestRef) = normalize(subSpan(TDCG, j, OP));
            oldestRef = iwrap(oldestRef + 1, nkPhases);
        end

        outPhases(i) = j;
        outCorrs(i) = maxCorr;

        if verbosity >= 3
            figure(5);
            subplot(1, 2, 1);
            plot(refSeqs);
            title('Running references');
            subplot(1, 2, 2);
            plot(mean(refSeqs, 2));
            title('Mean sequence');
            sgtitle(sprintf('Forward alignment (i=%d, j=%d)', i, j));
            waitforbuttonpress;
        end
    end
    % Trim empty buffer space off the end.
    outPhases = outPhases(1:i);
    outCorrs = outCorrs(1:i);

    % Align previous signal periods, updating the mean sequence, until less
    % than a period of the signal remains to be aligned.
    i = iCenter - k;
    j = outPhases(i);
    refSeqs = flip(nRefs, 2); % A circular buffer of reference sequences to average.
    oldestRef = 1; % A pointer to the first reference to remove.
    while j - OP >= 1
        i = i - 1;
        j = j - OP;

        j1 = max(j - searchRadius, 1);
        j2 = j + searchRadius + OP - 1;
        weightedRef = (avgRef + mean(refSeqs, 2))/2;
        [lag, maxCorr] = findMaxCorrLP(weightedRef, TDCG(j1:j2));

        if maxCorr > correlationThreshold
            j = j1 + lag; % Update j to point to the max correlation.

            % Update the running reference sequence.
            refSeqs(:, oldestRef) = normalize(subSpan(TDCG, j, OP));
            oldestRef = iwrap(oldestRef + 1, nkPhases);
        end

        outPhases(i) = j;
        outCorrs(i) = maxCorr;

        if verbosity >= 3
            figure(5);
            subplot(1, 2, 1);
            plot(refSeqs);
            title('Running references');
            subplot(1, 2, 2);
            plot(mean(refSeqs, 2));
            title('Mean sequence');
            sgtitle(sprintf('Backward alignment (i=%d, j=%d)', i, j));
            waitforbuttonpress;
        end
    end
    % Trim empty buffer space off the beginning.
    outPhases = outPhases(i:end);
    outCorrs = outCorrs(i:end);

    if verbosity >= 3
        close 5;
    end

    if verbosity >= 1
        fprintf('Alignment took %.3f seconds.\n', toc);
    end

    % Calculate the actual voxel depth, considering how scans may overlap.
    % Note: we count how many *gaps* there are between indices to describe
    % width, not how many indices there are.
    umPerPixelZ = umPerBScan * (outPhases(end) - outPhases(1)) / (numel(outPhases) - 1);
    fprintf('μm/pixel (Z) = %0.5f\n', umPerPixelZ);

    % The median can yield a fraction when two potential outputs are equal.
    % When this happens, choose the integer closest to OP.
    alignedPeriod = roundTowards(median(diff(outPhases)), OP);

    if verbosity >= 1
        fprintf('The median aligned period was %d frames, ', alignedPeriod);
        fprintf('as compared with the initial estimate of %d frames.\n', OP);
    end
    
    % When they differ, the median period is typically more accurate than
    % the initial choice of OP. Having a difference is rare unless you
    % choose a bad initial range to search for OP..
    OP = alignedPeriod;

    if verbosity >= 2
        figure('Name', 'Period statistics');
        subplot(1, 2, 1);
        plot(diff(outPhases));
        xlabel('Output Z frame index');
        ylabel('Frames');
        subplot(1, 2, 2);
        histogram(diff(outPhases));
        xlabel('Frames');
        sgtitle('Alignment step sizes');

        figure;
        plot(outCorrs);
        xlabel('Output Z frame index');
        ylabel('Correlation [-1 +1]');
        title('Best sequence correlations');

        zOut = umPerPixelZ * (0 : numel(outPhases) - 1)';
        zTrue = umPerBScan * (outPhases - outPhases(1));
        zError = zOut - zTrue;
        figure;
        yyaxis left;
        plot(zError);
        ylabel('Deviation (μm)');
        ylimitsLeft = gca().YLim;
        yyaxis right;
        ylim(ylimitsLeft / umPerPixelZ);
        ylabel('Deviation (pixels)');
        xlabel('Output Z frame index');
        title('Distortion in Z-axis');
    end

    %% Write output files.
    tic;
    firstImagePath = outputImages(outPhases, channelPaths, imageShape, ...
        imageIndices, rawStackOffsets, OP, outPath);

    % Save a matrix file listing the indices of each frame.
    saveImageIndices(outPhases, imageIndices, OP, outPath);

    if verbosity >= 1
        fprintf('Writing aligned image data to disk took %.3f seconds.\n', toc);
    end

    % Save all open figures to files by their titles and IDs.
    if saveFigures
        saveAllFigures(outPath);
    end

    tic;

    % Convert the file stack to an Imaris file.
    voxelSize = [umPerPixelX umPerPixelY umPerPixelZ];
    convertStackToImaris(imarisConvertPath, outPath, firstImagePath, voxelSize);

    if verbosity >= 1
        fprintf('Running ImarisConvert took %.3f seconds.\n', toc);
    end
end

function TDCG = loadTDCG(stackPath, imageShape, imageCrop, imageIndices, rawStackOffset)
    global verbosity;
    
    frameSize = prod(imageShape);
    
    %% Load the image stack and compute the TDCG.
    % TDCG is the average squared gradient of each image frame w.r.t. time.
    TDCG = zeros(numel(imageIndices) - 1, 1);
    
    if verbosity >= 1
        fprintf('Loading structure images...\n');
        progressBar = waitbar(0, 'Loading structure images...');
        progressBarCloser = onCleanup(@() close(progressBar));
    end
    
    fid = fopen(stackPath, 'r');
    fileCloser = onCleanup(@() fclose(fid));
    
    % Load the first frame.
    fseek(fid, rawStackOffset + (imageIndices(1)-1)*frameSize, 'bof');
    prevFrame = fread(fid, frameSize, 'uint8=>uint8');
    prevFrame = reshape(prevFrame, [imageShape(2) imageShape(1)])';
    prevFrame = prevFrame(imageCrop(1):imageCrop(3), imageCrop(2):imageCrop(4));
    prevFrame = im2double(prevFrame);
    
    % Load and process all subsequent frames.
    for i = 2:numel(imageIndices)
        if verbosity >= 1 && mod(i, 100) == 0
            progress = i/numel(imageIndices);
            waitbar(progress, progressBar);
        end
        
        fseek(fid, rawStackOffset + (imageIndices(i)-1)*frameSize, 'bof');
        frame = fread(fid, frameSize, 'uint8=>uint8');
        frame = reshape(frame, [imageShape(2) imageShape(1)])';
        frame = frame(imageCrop(1):imageCrop(3), imageCrop(2):imageCrop(4));
        frame = im2double(frame);
        
        TDCG(i-1) = mean2((frame - prevFrame).^2);
        
        prevFrame = frame;
    end
end

function [OP, jPeak] = determineOP(TDCG, expectedHeartPeriod)
    global verbosity;
    
    %% Decide which images to use when estimating the best output period.
    % The magnitude of the time gradient within a period indicates whether that
    % period contains periodic motion.

    % Generate a FIR low pass filter to limit the maximum frequency of the
    % smoothed TDCG to approximately 1/expectedHeartPeriod.
    % Targeting at least 50 dB attenuation of high frequencies.
    LP_taps = firpm(2*expectedHeartPeriod, ...
        [0 0 2/expectedHeartPeriod 1], ...
        [1 1 0 0], ...
        [100 1]);

    % Exclude the edges of the filtered region since they cannot be filtered
    % appropriately.
    TDCG_LPF = conv(TDCG, LP_taps, 'valid');
    TDCG_LPF = padarray(TDCG_LPF, floor(numel(LP_taps)/2), 'replicate', 'both');
    
    % The highpass-filtered TDCG.
    TDCG_HPF = TDCG - TDCG_LPF;
    
    if verbosity >= 2
        figure;
        hold on;
        plot(1:numel(TDCG), TDCG, 'c');
        plot(1:numel(TDCG), TDCG_LPF, 'b');
        xlabel('Input frame index');
        ylabel('TDCG [0-1]');
        legend('TDCG', 'Smoothed TDCG');
        title('Average Time Gradient');
    end

    % The H-score is a metric of determining how likely an index is to
    % contain the center of the heart.
    
    % The H-score has three components:
    
    % 1. The distance from the center of the scan, which is weighted
    % as 1 in the center and falling off to 0 at the edges.
    centerness = 1 - ((0 : numel(TDCG)-1)' / (numel(TDCG)-1) * 2 - 1).^4;
    
    % 2. The smoothed amplitude of the signal, i.e., the average motion.
    minAmp = min(TDCG_LPF);
    maxAmp = max(TDCG_LPF);
    scaledAmp = (TDCG_LPF - minAmp) / (maxAmp - minAmp);
    
    % 3. The smoothed amplitude of the variance of the signal, i.e., how
    % much deviation there is between the fast and slow moving regions.
    LP_taps = firpm(6*expectedHeartPeriod, ...
        [0 0 0.75/expectedHeartPeriod 1], ...
        [1 1 0 0], ...
        [100 1]);
    TDCG_var = TDCG_HPF.^2;
    TDCG_var = conv(TDCG_var, LP_taps, 'valid');
    TDCG_var = padarray(TDCG_var, floor(numel(LP_taps)/2), 'replicate', 'both');
    
    minVar = min(TDCG_var);
    maxVar = max(TDCG_var);
    scaledVar = (TDCG_var - minVar) / (maxVar - minVar);
    
    hScore = centerness .* scaledAmp .* scaledVar;
    
    % The half-way point between min and max will be the threshold.
    periodThresh = (min(hScore) + max(hScore))/2;

    validityCriteria = hScore >= periodThresh;
    validStart = find(validityCriteria, 1, 'first');
    validEnd = find(validityCriteria, 1, 'last');

    if verbosity >= 2
        figure;
        hold on;
        plot(1:numel(TDCG), scaledAmp);
        plot(1:numel(TDCG), scaledVar);
        plot(1:numel(TDCG), centerness);
        xlabel('Input frame index');
        ylabel('Scaled Value [0-1]');
        legend('Scale Amplitude', 'Scale Variance', 'Centerness');
        title('Calculating the H-Score');
        
        figure;
        hold on;
        plot(1:numel(TDCG), hScore, 'k');
        plot([1 numel(TDCG)], [periodThresh periodThresh], 'r');
        plot([validStart validEnd], hScore([validStart validEnd]), 'rx');
        xlabel('Input frame index');
        ylabel('Scaled Value [0-1]');
        legend('H-Score', 'Threshold', 'Bounds');
        title('Likelihood of Heart Center');
    end
    
    %% Determine the repetition period of each sequence.
    % Estimating the period using only high frequencies is to reduce the
    % bias correlation due to changing average signal amplitude.
    valid = TDCG_HPF(validStart:validEnd);
    
    % Search for the ideal period withn 90-110% of the expected period.
    minPeriod = floor(0.9 * expectedHeartPeriod);
    minPeriod = max(minPeriod, 2);
    
    maxPeriod = ceil(1.1 * expectedHeartPeriod);
    maxPeriod = min(maxPeriod, numel(valid)-1);

    if verbosity >= 2
        fprintf('Searching for periods between %d and %d frames.\n', ...
            minPeriod, maxPeriod);
    end
    
    ref = valid(1 : end - maxPeriod);
    
    nTests = maxPeriod - minPeriod + 1;
    correlations = zeros(nTests, 1);
    for i = 1:nTests
        testPeriod = minPeriod + i - 1;
        y = valid(1 + testPeriod : end + testPeriod - maxPeriod);
        correlations(i) = mean(ref .* y);
    end
    
    [~, iMax] = max(correlations);
    
    % - 1: convert 1-based index to 0-base.
    OP = minPeriod + iMax - 1;

    if verbosity >= 1
        fprintf('Determined an output period of %d frames.\n', OP);
    end
    
    if verbosity >= 2
        figure;
        plot(minPeriod : maxPeriod, correlations);
        xlabel('Test Period (# frames)');
        ylabel('Correlation [-1 +1]');
        title('Period Correlations');
    end
    
    % Find the maximum H-Score to ensure a good initial alignment. Within
    % that maximal H-Score, start with the fastest moving peak for a
    % consistent choice of reference phase. Hearts will always align with
    % their fastest moving region at phase = 0.
    [~, jPeak] = max(hScore);
    
    % A peak too close to one end of the scan will result in an indexing
    % error. If the choice of jPeak is too close to either end, ignore the
    % original selection and just pick from the center of the scan.
    if jPeak < 3*OP || jPeak + 3*OP > numel(TDCG)
        if verbosity >= 1
            fprintf(['The heart activity peak identified around frame #%d' ...
                ' was discarded due to proximity to the edge of the scan.\n'], jPeak);
        end
        jPeak = 1 + floor((numel(TDCG) - OP)/2);
    end
    
    [~, jPeakOffset] = max(subSpan(TDCG, jPeak - floor(OP/2), OP));
    jPeak = jPeak + floor(OP/2) + jPeakOffset - 1;
    
    if verbosity >= 1
        fprintf('Peak heart activity identified at frame #%d.\n', jPeak);
    end
end

function firstImagePath = outputImages(outPhases, channelPaths, ...
        imageShape, imageIndices, rawStackOffsets, OP, outPath)
    global verbosity;
    
    frameSize = prod(imageShape);
    
    % Ensure that the output directory exists.
    if ~exist(outPath, 'dir')
        mkdir(outPath);
        
        if verbosity >= 1
            fprintf('Created the stack output directory at "%s".\n', outPath);
        end
    end
    
    % Create a string template with the appropriate number of digits for
    % how many frames there are to write.
    outFilePathTemplate = [ ...
        outPath '/Sync' ...
        '_T%0' int2str(countDigits(OP - 1)) 'd' ...
        '_C%0' int2str(countDigits(numel(channelPaths) - 1)) 'd' ...
        '_Z%0' int2str(countDigits(numel(outPhases) - 1)) 'd.tif' ...
    ];

    firstImagePath = sprintf(outFilePathTemplate, 0, 0, 0);

    if verbosity >= 1
        fprintf('Exporting synchronized image sequence...\n');
        progressBar = waitbar(0, 'Exporting synchronized sequence...');
        progressBarCloser = onCleanup(@() close(progressBar));
    end
    % Channel index, ic (zero-base index).
    for ic = 0:numel(channelPaths)-1
        stackPath = channelPaths{ic + 1};
        fid = fopen(stackPath, 'r');
        fileCloser = onCleanup(@() fclose(fid));
        
        % Output Z index, iz (zero-base index).
        for iz = 0:numel(outPhases)-1
            if verbosity >= 1
                progress = (ic + iz/numel(outPhases))/numel(channelPaths);
                waitbar(progress, progressBar);
            end
            
            % Phase index, ip (zero-base index).
            for ip = 0:OP-1
                % The frame index into the RAW stack (zero-based index).
                iRaw = imageIndices(1) + outPhases(iz + 1) - 2 + ip;
                
                fseek(fid, rawStackOffsets(ic+1) + iRaw*frameSize, 'bof');
                frame = fread(fid, frameSize, 'uint8=>uint8');
                frame = reshape(frame, [imageShape(2) imageShape(1)])';
                imwrite(frame, sprintf(outFilePathTemplate, ip, ic, iz));
            end
        end
    end
end

function saveImageIndices(outPhases, imageIndices, OP, outPath)
    global verbosity;
    
    % Shape: [Z-depth, time]
    imageIndices = outPhases + (imageIndices(1) + (0:OP-1));

    % Save the (z, t)->i image file indices for comparison.
    [outRoot, outName, ~] = fileparts(outPath);
    if ~isempty(outRoot)
        outRoot = [outRoot '/'];
    end
    
    indexFilePath = sprintf('%simageIndices_%s.mat', outRoot, outName);
    
    if verbosity >= 1
        fprintf('Wrote frame index file to "%s".\n', indexFilePath);
    end
    
    save(indexFilePath, 'imageIndices');
end

function saveAllFigures(outPath)
    global verbosity;
    
    [outRoot, outName, ~] = fileparts(outPath);
    if ~isempty(outRoot)
        outRoot = [outRoot '/'];
    end
    
    figureDir = [outRoot outName '_figures'];
    
    % Ensure that the output directory exists.
    if ~exist(figureDir, 'dir')
        mkdir(figureDir);
        
        if verbosity >= 1
            fprintf('Created the figure output directory at "%s".\n', figureDir);
        end
    end
    
    for fig = findall(groot(), 'Type', 'figure')'
        figureName = fig.Name;
        if isempty(figureName)
            figureName = fig.CurrentAxes.Title.String;
        end
        figPath = sprintf('%s/%d_%s.fig', figureDir, fig.Number, figureName);
        savefig(fig, figPath);
        saveas(fig, [figPath(1:end-3) 'png']);
    end
end

function convertStackToImaris(imarisConvertPath, outPath, firstImagePath, voxelSize)
    global verbosity;
    
    [outRoot, outName, ~] = fileparts(outPath);
    if ~isempty(outRoot)
        outRoot = [outRoot '/'];
    end

    if ~isempty(imarisConvertPath)
        if verbosity >= 1
            fprintf('Converting stack to an Imaris file...\n');
        end

        ImarisCommand = strjoin({
            sprintf('"%s"', imarisConvertPath)
            sprintf('--input "%s"', firstImagePath)
            sprintf('--output "%s%s.ims"', outRoot, outName)
            sprintf('--voxelsize "%0.5f %0.5f %0.5f"', voxelSize)
            %'--nthreads 20' % Use 20 threads for image compression (default is 8).
        }, " ");

        if verbosity >= 1
            fprintf('=> %s\n', ImarisCommand);
        end

        % Execute the Imaris Converter command.
        errorCode = system(ImarisCommand);

        if verbosity >= 1
            if errorCode
                fprintf('Imaris file conversion failed with error code %d.\n', errorCode);
            else
                fprintf('Imaris file conversion was successful!\n');
            end
        end
    end
end

% Slide 'ref' along 'space' recording the correlation between the
% overlapping regions. Returns the lag of the maximal correlation.
function [lag, maxCorr] = findMaxCorrLP(ref, space)
    assert(numel(space) >= numel(ref) + 2, ...
        'The search space is not large enough to compare with the reference.');

    width = numel(ref);
    nOffsets = numel(space) - width;
    correlations = zeros(nOffsets, 1);
    ref = normalize(ref); % Normalizing the reference does nothing.
    for i = 1:nOffsets
        y = space(i : i + width - 1);
        % It's important to normalize each exerpt from the space.
        y = normalize(y);
        correlations(i) = mean(ref .* y);
    end
    
    % This applies a 3-tap LP filter to the correlation map.
    % The filter smooths the data to make 'plateaus' convex.
    LP_filter_taps = [0.3 1 0.3] / 1.6;
    corrLPfiltered = conv(correlations, LP_filter_taps, 'valid');
    
    [maxCorr, iMax] = max(corrLPfiltered);
    
    % - 1: convert 1-based index to 0-base.
    % + 1: reverse the time-lag caused by the LP filter.
    lag = iMax - 1 + 1;
end

% Wraps indices, accounting for MATLAB's ridiculous 1-based index scheme.
function j = iwrap(i, width)
    j = mod(i - 1, width) + 1;
end

% Index the first 'width' elements from an 'array' starting at 'index'.
function span = subSpan(array, index, width)
    span = array(index : index + width - 1);
end

% Rounds the floating point number 'x' so that it is the integer closest to
% the integer 'y'.
function x = roundTowards(x, y)
    x = y + fix(x - y);
end

% Count the number of decimal digits in the integer, 'x'.
function n = countDigits(x)
    assert(x >= 0, 'Cannot count digits of a negative number.');
    
    if x == 0
        n = 1;
    else
        n = 1 + floor(log10(x));
    end
end