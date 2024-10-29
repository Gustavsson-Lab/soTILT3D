function [] = f_combineMulticolorOutput()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

scrsz = get(0,'ScreenSize');
numSyncFrames = 25;
useDenoising = 1;
nSample = 1.33;         % index of refraction of sample
nOil = 1.518;           % index of immersion oil
registrationComplete = false;

%% ask user for relevant datafiles
[dataFile, dataPath] = uigetfile({'*.mat'},'Open matlab file from a previous run');
if ~isequal(dataFile,0)
    dataFile = [dataPath dataFile];
    load(dataFile)
end


if ~registrationComplete
    
    %% Ask user for relevant datafiles
    
    dlg_title = 'What channels would you like to transform?';
    prompt = {'Enter channels as 1-letter abbreviations (g,r) without spaces',...
              'Channel to transform into (enter 1-letter abbreviation: g or r)',...
              'Are the channels interleaved (0) or sequential (1)?'...
             };
    def = {'gr','g','1'}; % e.g. if 320x320 of 512x512 in far corner, use 193, 193
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    channelsToCombine = inputdialog{1};
    referenceChannel = inputdialog{2};
    sequentialMeas = str2num(inputdialog{3});
    numChannels = length(channelsToCombine);
    
    color = cell(0);
    for colorIdx = 1:numChannels
        if strcmpi(channelsToCombine(colorIdx),'g')
            color{colorIdx} = 'Green';
        elseif strcmpi(channelsToCombine(colorIdx),'r')
            color{colorIdx} = 'Red';
        else
            color{colorIdx} = 'Black';
        end
    end
    [tformFile tformPath] = uigetfile({'*.mat';'*.*'},'Open 3D_Transform.mat');
    if isequal(tformFile,0)
        error('User cancelled the program');
    end
    load([tformPath tformFile]);

    %% Prepare for data interpolation
    if ~exist('evalCP')
        evalCP = (1:size(matched_cp_reflected,1))';
    end
    x = matched_cp_reflected(evalCP,5);
    y = matched_cp_reflected(evalCP,6);
    z = matched_cp_reflected(evalCP,7);
    F_FRE = TriScatteredInterp(x,y,z,FRE_full(:,1), 'natural');
    F_FRE_X = TriScatteredInterp(x,y,z,FRE_full(:,2), 'natural');
    F_FRE_Y = TriScatteredInterp(x,y,z,FRE_full(:,3), 'natural');
    F_FRE_Z = TriScatteredInterp(x,y,z,FRE_full(:,4), 'natural');
    F_TRE = TriScatteredInterp(x,y,z,TRE_full(:,1), 'natural');
    F_TRE_X = TriScatteredInterp(x,y,z,TRE_full(:,2), 'natural');
    F_TRE_Y = TriScatteredInterp(x,y,z,TRE_full(:,3), 'natural');
    F_TRE_Z = TriScatteredInterp(x,y,z,TRE_full(:,4), 'natural');    
     
    locFiles = {};
    dataSets = [];
    useFids = 1; % use fids unless one or both files does not have a fiducial
    
    for fileNum = 1:numChannels
        
        [locFile, locPath] = uigetfile({'*.mat';'*.*'},['Open data file #' num2str(fileNum) ' ('...
            channelsToCombine(fileNum) ' channel) with filtered molecule fits']);
        if isequal(locFile,0)
            error('User cancelled the program');
        end
        
        locFiles = [locFiles; {[locPath locFile]}];
        % load data
        load([locPath locFile]);
        
        %    Assemble the structure for this dataset
        dataSet.frameNum = frameNum;
        
        dataSet.sigmaX = sigmaX;
        dataSet.sigmaY = sigmaY;
        dataSet.sigmaZ = sigmaZ;
        dataSet.numPhotons = numPhotons;
        dataSet.meanBkgnd = meanBkgnd;
        
        dataSet.fidTrackX = fidTrackX;
        dataSet.fidTrackY = fidTrackY;
        dataSet.fidTrackZ = fidTrackZ;
        dataSet.wlShift = [wlShiftX, wlShiftY];
        dataSet.fidsToUse = fidsToUse;
        dataSet.fidChannel = fidChannel;
        
        % set all not-in-use fiducials to NaNs (could strip out, but want
        % to keep array same size to be consistent)
        clearedFids = ~ismember(1:size(fidTrackX,2),fidsToUse);
        dataSet.fidTrackX(:,clearedFids) = nan;
        dataSet.fidTrackY(:,clearedFids) = nan;
        dataSet.fidTrackZ(:,clearedFids) = nan;
        
        if isempty(fidTrackX)
            dataSet.fidCorrected = 0;
            useFids = 0;
            dataSet.xLoc = xLoc;
            dataSet.yLoc = yLoc;
            dataSet.zLoc = zLoc;
        else
            dataSet.fidCorrected = 1;
            dataSet.xLoc = xLocRaw;
            dataSet.yLoc = yLocRaw;
            dataSet.zLoc = zLocRaw;
        end
        
        if exist('sigma_x_nm','var') 
            dataSet.laserPeakInt = peakIntensity;
            dataSet.laser_x_nm = laser_x_nm;
            dataSet.laser_y_nm = laser_y_nm;
            dataSet.lasersig_x_nm = sigma_x_nm;
            dataSet.lasersig_y_nm = sigma_y_nm;
            clear('sigma_x_nm','sigma_y_nm','laser_x_nm','laser_y_nm','peakIntensity');
        end
        
%         dlg_title = 'Transform dataset';
%         prompt = {'Do you want to transform this dataset?'};
%         def =       { 'Yes' };
%         questiondialog = questdlg(prompt,dlg_title, def);
        dataToTransform = ~strcmpi(channelsToCombine(fileNum),referenceChannel);
        
        % Handle response
        if dataToTransform
                
                % transform SM data
                
                if dataSet.fidCorrected
                    transformedData = transformData([xLocRaw, yLocRaw, zLocRaw],tform);
                else
                    transformedData = transformData([xLoc, yLoc, zLoc],tform);
                end
                dataSet.xLoc_transformed = transformedData(:,1);
                dataSet.yLoc_transformed = transformedData(:,2);
                dataSet.zLoc_transformed = transformedData(:,3);
                transformedDataSet = fileNum;
                
        elseif ~dataToTransform

                dataSet.xLoc_transformed = NaN;
                dataSet.yLoc_transformed = NaN;
                dataSet.zLoc_transformed = NaN;
                dataSet.fidTrackX_transformed =NaN;
                dataSet.fidTrackY_transformed = NaN;
                dataSet.fidTrackZ_transformed = NaN;
                dataSet.fidTrack_interpolated_FRE = NaN;
                dataSet.fidTrack_interpolated_TRE = NaN;
                untransformedDataSet = fileNum;
                %             dataSet.xLoc_transformed = NaN(length(xLoc),1);
                %             dataSet.yLoc_transformed = NaN(length(xLoc),1);
                %             dataSet.zLoc_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrackX_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrackY_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrackZ_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrack_interpolated_FRE = NaN(length(fidTrackX),4);
                %             dataSet.fidTrack_interpolated_TRE = NaN(length(fidTrackX),4);
        end
        
        if dataSet.fidCorrected
            % transform fiducial data
            % Accepts fiducial being different 
            fidToTransform = ~strcmpi(fidChannel,referenceChannel);
            numFids = length(fidsToUse);%size(fidTrackX,2);
            dataSet.fidTrackX_transformed = NaN(size(fidTrackX));
            dataSet.fidTrackY_transformed = NaN(size(fidTrackX));
            dataSet.fidTrackZ_transformed = NaN(size(fidTrackX));
            if fidToTransform
                for fidNum = fidsToUse %1:numFids
                transformedData = transformData([fidTrackX(~isnan(fidTrackX(:,fidNum)),fidNum),...
                    fidTrackY(~isnan(fidTrackY(:,fidNum)),fidNum),...
                    fidTrackZ(~isnan(fidTrackZ(:,fidNum)),fidNum)],tform);
                dataSet.fidTrackX_transformed(~isnan(fidTrackX(:,fidNum)),fidNum) = transformedData(:,1);
                dataSet.fidTrackY_transformed(~isnan(fidTrackY(:,fidNum)),fidNum) = transformedData(:,2);
                dataSet.fidTrackZ_transformed(~isnan(fidTrackZ(:,fidNum)),fidNum) = transformedData(:,3);
                end
            end
        end
        
        
        dataSet.LocFile = locFile;
        dataSet.LocPath = locPath;
        
        dataSet.frameRange = frameRange;
        dataSet.zRange = zRange;
        dataSet.ampRatioLimit = ampRatioLimit;
        dataSet.fitErrorRange = fitErrorRange;
        dataSet.lobeDistBounds = lobeDistBounds;
        dataSet.sigmaBounds = sigmaBounds;
        dataSet.sigmaRatioLimit = sigmaRatioLimit;
        dataSet.numPhotonRange = numPhotonRange;
        dataSet.fidsToUse = fidsToUse;
        dataSet.dataTransformed = dataToTransform;
        dataSet.fidTransformed = fidToTransform;
        
        % Append the structure to the previous structures in the array
        dataSets = [dataSets, dataSet];
        
        
    end
    
    clear LocFile LocPath ampRatioLimit dataSet def dlg_title fileNum
    clear fitErrorRange frameNum frameRange lobeDistBounds meanBkgnd
    clear numPhotonRange numPhotons prompt questiondialog sigmaBounds
    clear sigmaRatioLimit sigmaX sigmaY sigmaZ transformedData
    clear xLoc yLoc zLoc zRange fidTrackX fidTrackY fidTrackZ
    % save('workspace.mat');
    
    
if useFids
    
    %% Identify frames based on the sequence log
    % load data and register sequence log to data frames
    
    [logFile logPath] = uigetfile({'*.dat'},'Open shutter sequence log file');
%     if isequal(logFile,0)
%         error('User cancelled the program');
%     end
    useLog = ~isequal(logFile,0);
    if useLog
        sifLogData =  importdata([logPath logFile]);
        % check: is this complicated proofreading of length really necessary?
        if useFids
            sifLogData = sifLogData(1:size(dataSets(1).fidTrackX,1),:);
        else
            sifLogData = sifLogData(1:max(vertcat(dataSets(:).frameNum)),:);
        end
    end
%     
%     for k = 1:size(sifLogData,1) % kludge fix - should not be needed!!
%         if sifLogData(k,3)==1 && sifLogData(k,2)==1
%             sifLogData(k,1:3) = [0 1 1];
%         elseif sifLogData(k,3)==1
%             sifLogData(k,1:3) = [0 0 1];
%         elseif sifLogData(k,2)==1
%             sifLogData(k,1:3) = [0 1 0];
%         else
%             sifLogData(k,1:3) = [nan nan nan];
%         end
%     end
%     sifLogData(:,4) = [];
%     for k = 1:size(sifLogData,1)
%         if isnan(sifLogData(k,1)) && k < size(sifLogData,1)
%             sifLogData(k,:) = sifLogData(k+1,:);
%         else
%             sifLogData(k,1:3) = [0 0 0];
%         end
%     end

    % switch to logical from integer index
    % define frames to use, either defined by log (interleaved imaging),
    % or order of files (sequential imaging)
    if useLog
        frames = (1:length(sifLogData))';
        greenFrames = sifLogData(:,2) == 1;
        redFrames = sifLogData(:,3) == 1;
            % find shutter transition frames
        startG = find([false; greenFrames(2:end)& ~greenFrames(1:end-1)]);
        endG = find([greenFrames(1:end-1)& ~greenFrames(2:end);false]);

        startR = find([false; redFrames(2:end)& ~redFrames(1:end-1)]);
        endR = find([redFrames(1:end-1)& ~redFrames(2:end);false]);        
    else
        if useFids && ~sequentialMeas
            frames = 1:size(dataSets(1).fidTrackX,1);            
        elseif useFids && sequentialMeas
            frames = 1:(size(dataSets(1).fidTrackX,1)+size(dataSets(2).fidTrackX,1));
        else
            frames = 1:max(vertcat(dataSets(:).frameNum));
        end
    end

    
    % test to make sure frames are called correctly
%     figure; hold on; plot(greenFrames,'g'); plot(redFrames,'r');
%     scatter(startG,ones(size(startG)),'g'); scatter(endG,ones(size(endG)),30,[0 0.5 0]);
%     scatter(startR,ones(size(startR)),'r'); scatter(endR,ones(size(endR)),30,[0.5 0 0]);
%     xlim([100 400]); ylim([0.95 1.05])
    
% old version
%     endFrames_red = redFrames(find(diff(redFrames)== mean(temp(temp>4))));
%     startFrames_green = endFrames_red+1; % this seems very likely to fail if we make any changes to the acquisition parameters....
    
    
    %% Show the decomposed fiducial tracks
    
    fidTracksX = [];
    fidTracksY = [];
    fidTracksZ = [];
    
    % average fiducials from each channel
    if sequentialMeas
        channelFrameOffset = nan(length(dataSets),1);
        for i = 1:length(dataSets)
            channelFrameOffset(i) = length(fidTracksX);
            if dataSets(i).fidTransformed
                fidTracksX = [fidTracksX; nanmean(dataSets(i).fidTrackX_transformed,2)];
                fidTracksY = [fidTracksY; nanmean(dataSets(i).fidTrackY_transformed,2)];
                fidTracksZ = [fidTracksZ; nanmean(dataSets(i).fidTrackZ_transformed,2)];
            else
                fidTracksX = [fidTracksX; nanmean(dataSets(i).fidTrackX,2)];
                fidTracksY = [fidTracksY; nanmean(dataSets(i).fidTrackY,2)];
                fidTracksZ = [fidTracksZ; nanmean(dataSets(i).fidTrackZ,2)];
            end
        end
    else
        channelFrameOffset = zeros(length(dataSets),1);
        for i = 1:length(dataSets)
            if dataSets(i).fidTransformed
%                 transformedDataSet = i;
                fidTracksX = [fidTracksX, nanmean(dataSets(i).fidTrackX_transformed,2)];
                fidTracksY = [fidTracksY, nanmean(dataSets(i).fidTrackY_transformed,2)];
                fidTracksZ = [fidTracksZ, nanmean(dataSets(i).fidTrackZ_transformed,2)];
            else
%                 untransformedDataSet = i;
                fidTracksX = [fidTracksX, nanmean(dataSets(i).fidTrackX,2)];
                fidTracksY = [fidTracksY, nanmean(dataSets(i).fidTrackY,2)];
                fidTracksZ = [fidTracksZ, nanmean(dataSets(i).fidTrackZ,2)];
            end
        end
    end
    
    for i = 1:length(dataSets)
        dataSets(i).frameOffset = channelFrameOffset(i);
        dataSets(i).fidFrames = (1:length(dataSets(i).fidTrackX)) + channelFrameOffset(i);
    end
    
    h_allRegisteredTracks = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2); 
    
    for i = 1:length(dataSets)
        if dataSets(i).fidTransformed
            subplot(3,1,1); hold on;
            plot(dataSets(i).fidFrames,dataSets(i).fidTrackX_transformed, color{i});
            xlabel('frame number');ylabel('X position (nm)');
            
            subplot(3,1,2); hold on;
            plot(dataSets(i).fidFrames,dataSets(i).fidTrackY_transformed, color{i});
            xlabel('frame number');ylabel('Y position (nm)');

            subplot(3,1,3); hold on;
            plot(dataSets(i).fidFrames,dataSets(i).fidTrackZ_transformed, color{i});
            xlabel('frame number');ylabel('Z position (nm)');
        else
            subplot(3,1,1); hold on;
            plot(dataSets(i).fidFrames,dataSets(i).fidTrackX, color{i});
            xlabel('frame number');ylabel('X position (nm)');

            subplot(3,1,2); hold on;
            plot(dataSets(i).fidFrames,dataSets(i).fidTrackY, color{i});
            xlabel('frame number');ylabel('Y position (nm)');

            subplot(3,1,3); hold on;
            plot(dataSets(i).fidFrames,dataSets(i).fidTrackZ, color{i});
            xlabel('frame number');ylabel('Z position (nm)');
        end
    end
    
    % display the mean track also
%     if numFids > 1
%     subplot(3,1,1)
%     plot(nanmean(fidTracksX,2), 'black');
%     xlabel('frame number');ylabel('X position (nm)');
%     hold off
%     subplot(3,1,2)
%     plot(nanmean(fidTracksY,2), 'black');
%     xlabel('frame number');ylabel('Y position (nm)');
%     hold off
%     subplot(3,1,3)
%     plot(nanmean(fidTracksZ,2), 'black');
%     xlabel('frame number');ylabel('Z position (nm)');
%     hold off
%     end

subplot(3,1,1); title('XYZ fiducial drift for each channel');
    %% Pick the appropriate track to use as drift correction
    
    if ~sequentialMeas
        dlg_title = 'Inspect Fiducial Tracks';
        prompt = {'Choose Fiducial Track to denoise and use to correct drift'};
        def =       { 'untransformed'  };
        questiondialog = questdlg(prompt,dlg_title,'untransformed','transformed','correct separately', def);
        % Handle response
        switch questiondialog
            case 'untransformed'
                    chosenFidTrack = untransformedDataSet;
            case 'transformed'
                    chosenFidTrack = transformedDataSet;
            case 'correct separately'
                    chosenFidTrack = nan;
            case 'Cancel'
                error('User cancelled the program');
        end
    else
        chosenFidTrack = 1; % there only exists one concatenated track
    end
    %% Denoise the raw data   
   
    dlg_title = 'Data denoising';
    prompt = {'Denoise data using wavelet (0) or boxcar (1)?',...
              'If using boxcar, what *radius* in frames?',...
             };
    def = {'0','2'}; % e.g. if 320x320 of 512x512 in far corner, use 193, 193
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    useBoxCar = logical(str2num(inputdialog{1}));
    bcRad = str2num(inputdialog{2}); % can we pin to localization precision?
    
    % If jumping between red and green illumination, denoising frames 
    % separately is not a good option - this approach forces each set of
    % tracks to be continuous - use boxcar instead
    
    if sequentialMeas
        if useBoxCar
            for j = 1:size(fidTracksX,2)
                for i = 1:max(dataSets(end).fidFrames)
                    if i <= bcRad||length(frames)-i<=bcRad
                        continue
                    end
                fidTracksX_denoised(i,j)=nanmean(fidTracksX(i-bcRad:i+bcRad,j));
                fidTracksY_denoised(i,j)=nanmean(fidTracksY(i-bcRad:i+bcRad,j));
                fidTracksZ_denoised(i,j)=nanmean(fidTracksZ(i-bcRad:i+bcRad,j));
                end
            end
        else
            [fidTracksX_denoised,fidTracksY_denoised,fidTracksZ_denoised]...
                = f_waveletFidTracks(fidTracksX,fidTracksY,fidTracksZ,1);%0);
        end
    elseif ~sequentialMeas
    
    
    if isnan(chosenFidTrack) % correct channels separately using relevant channel
    
        fidTracksX_denoised_sepFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
        fidTracksY_denoised_sepFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
        fidTracksZ_denoised_sepFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
    
    if useLog
        toSeparate = {'greenFrames','redFrames'};
    else
        toSeparate = {'logical(frames)'''}; % i.e. " logical(frames)' "
    end
    
    for frameSetIdx = 1:length(toSeparate)
        sepFrames = eval(toSeparate{frameSetIdx});
        if useBoxCar
        for j = 1:size(fidTracksX_denoised_sepFrames,2)
            for i = find(sepFrames)'
                if i <= bcRad||length(sepFrames)-i<=bcRad
                    continue
                end
            fidTracksX_denoised_sepFrames(i,j)=nanmean(fidTracksX(i-bcRad:i+bcRad,j));
            fidTracksY_denoised_sepFrames(i,j)=nanmean(fidTracksY(i-bcRad:i+bcRad,j));
            fidTracksX_denoised_sepFrames(i,j)=nanmean(fidTracksZ(i-bcRad:i+bcRad,j));
            end
        end
        else
        [fidTracksX_denoised_sepFrames(sepFrames,frameSetIdx),...
            fidTracksY_denoised_sepFrames(sepFrames,frameSetIdx),...
            fidTracksZ_denoised_sepFrames(sepFrames,frameSetIdx)] = f_waveletFidTracks(...
            fidTracksX(sepFrames,frameSetIdx),...
            fidTracksY(sepFrames,frameSetIdx),...
            fidTracksZ(sepFrames,frameSetIdx),1);%0);
        end
    end
      

        fidTracksX_denoised = fidTracksX_denoised_sepFrames;
        fidTracksY_denoised = fidTracksY_denoised_sepFrames;
        fidTracksZ_denoised = fidTracksZ_denoised_sepFrames;
        
    else % if using a selected fid to do all corrections
        if useBoxCar
            error('box car not implemented for this mode yet')
        else
         [fidTracksX_denoised,fidTracksY_denoised,fidTracksZ_denoised]...
                = f_waveletFidTracks(fidTracksX(:,chosenFidTrack),...
                   fidTracksY(:,chosenFidTrack),fidTracksZ(:,chosenFidTrack),1);
        end
    end
    end %% TODO - check to make sure adding this was correct
    save('workspace.mat')      
    
    
    %% Center the last frames of the chosen denoised fiducial track

%     syncFrames = find(~isnan(fidTracksX_shifted_denoised(:,chosenFidTrack)));
%     syncFrames = syncFrames(end-(numSyncFrames-1):end);
    
        avgDevX = fidTracksX_denoised(:,chosenFidTrack) - nanmean(fidTracksX_denoised(:,chosenFidTrack));
        avgDevY = fidTracksY_denoised(:,chosenFidTrack) - nanmean(fidTracksY_denoised(:,chosenFidTrack));
        avgDevZ = fidTracksZ_denoised(:,chosenFidTrack) - nanmean(fidTracksZ_denoised(:,chosenFidTrack));

        % Show the drift correction to be applied
%     close all; % close extraneous figures

%     color = {'Green', 'Red', 'Blue', 'Cyan', 'Yellow' };
    h_finalCorrections = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);

        subplot(3,1,1)
        plot(avgDevX, color{chosenFidTrack});
        xlabel('frame number');ylabel('X correction (nm)');

        subplot(3,1,2)
        plot(avgDevY, color{chosenFidTrack});
        xlabel('frame number');ylabel('Y correction (nm)');

        subplot(3,1,3)
        plot(avgDevZ, color{chosenFidTrack});
        xlabel('frame number');ylabel('Z correction (nm)');
    
        %% Apply fiduciary corrections

    for i = 1:length(dataSets)
%         dataSets(i).driftCorrX = avgDevX;
%         dataSets(i).driftCorrY = avgDevY;
%         dataSets(i).driftCorrZ = avgDevZ;
%         dataSets(i).devX_denoised = devX_denoised(:,i);
%         dataSets(i).devY_denoised = devY_denoised(:,i);
%         dataSets(i).devZ_denoised = devZ_denoised(:,i);
%         dataSets(i).fidTrackX_denoised = fidTracksX_denoised(:,i);
%         dataSets(i).fidTrackY_denoised = fidTracksY_denoised(:,i);
%         dataSets(i).fidTrackZ_denoised = fidTracksZ_denoised(:,i);

        if dataSets(i).dataTransformed
            dataSets(i).xLoc_driftCorr = dataSets(i).xLoc_transformed - avgDevX(dataSets(i).frameNum+dataSets(i).frameOffset);
            dataSets(i).yLoc_driftCorr = dataSets(i).yLoc_transformed - avgDevY(dataSets(i).frameNum+dataSets(i).frameOffset);
            dataSets(i).zLoc_driftCorr = dataSets(i).zLoc_transformed - avgDevZ(dataSets(i).frameNum+dataSets(i).frameOffset);
        else
            untransformedDataSet=i;
            dataSets(i).xLoc_driftCorr = dataSets(i).xLoc - avgDevX(dataSets(i).frameNum+dataSets(i).frameOffset);
            dataSets(i).yLoc_driftCorr = dataSets(i).yLoc - avgDevY(dataSets(i).frameNum+dataSets(i).frameOffset);
            dataSets(i).zLoc_driftCorr = dataSets(i).zLoc - avgDevZ(dataSets(i).frameNum+dataSets(i).frameOffset);
        end
    end

    % moved lower: this was originally here
%         %% Apply index mismatch corrections
%     % Todo: This is empirical. A better model accouting for index mismatch
%     % needs to developed here.
%     for i = 1:length(dataSets)
%         dataSets(i).zLoc_driftCorr_indexCorr = dataSets(i).zLoc_driftCorr * nSample/nOil;
%     end 

    %% Show the difference between the registered fiducial tracks
%     if ~exist('tformChan') % this currently doesn't work unless transform second dataset
%     
%     figure_h_c = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
%     set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
%         'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
%         'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
%         'DefaultAxesLineWidth',1.2);
%     
%     subplot(4,1,1)
%     plot(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
%     xlabel('frame number');ylabel('X position shift (nm)');
%     avg = nanmean(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
%     stdev = nanstd(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
%     ylim([avg-5*stdev avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,2),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_x']})
%     title({['Fused Fiducial Tracks (raw data)']})
%     
%     subplot(4,1,2)
%     plot(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
%     xlabel('frame number');ylabel('Y position shift (nm)');
%     avg = nanmean(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
%     stdev = nanstd(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
%     ylim([avg-5*stdev avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,3),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_y']})
%     
%     subplot(4,1,3)
%     plot(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
%     xlabel('frame number');ylabel('Z position shift (nm)');
%     avg = nanmean(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
%     stdev = nanstd(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
%     ylim([avg-5*stdev avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,4),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_z']})
%     
%     subplot(4,1,4)
%     euclid_Dist = sqrt(((dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed).^2)+...
%         ((dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed).^2+...
%         ((dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed).^2)));
%     plot(euclid_Dist);
%     xlabel('frame number');ylabel('3D shift (nm)');
%     avg = nanmean(euclid_Dist);
%     stdev = nanstd(euclid_Dist);
%     ylim([0 avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,1),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_3_D']})
%     clear avg
%     end
    %% Show the difference between the denoised and fused fiducial tracks
    % this uses the variables 'fidTracksX_denoised' and so on, which aren't
    % defined...?
%     figure_h_d = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
%     set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
%         'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
%         'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
%         'DefaultAxesLineWidth',1.2);
%     
%     subplot(4,1,1)
%     plot(fidTracksX_denoised(:,1)-fidTracksX_denoised(:,2));
%     xlabel('frame number');ylabel('X position shift (nm)');
%     avg = nanmean(fidTracksX_denoised(:,1)-fidTracksX_denoised(:,2));
%     stdev = nanstd(fidTracksX_denoised(:,1)-fidTracksX_denoised(:,2));
%     ylim([avg-5*stdev avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,2),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_x']})
%     title({['Fused Fiducial Tracks (denoised)']})
%     
%     subplot(4,1,2)
%     plot(fidTracksY_denoised(:,1)-fidTracksY_denoised(:,2));
%     xlabel('frame number');ylabel('Y position shift (nm)');
%     avg = nanmean(fidTracksY_denoised(:,1)-fidTracksY_denoised(:,2));
%     stdev = nanstd(fidTracksY_denoised(:,1)-fidTracksY_denoised(:,2));
%     ylim([avg-5*stdev avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,3),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_y']})
%     
%     subplot(4,1,3)
%     plot(fidTracksZ_denoised(:,1)-fidTracksZ_denoised(:,2));
%     xlabel('frame number');ylabel('Z position shift (nm)');
%     avg = nanmean(fidTracksZ_denoised(:,1)-fidTracksZ_denoised(:,2));
%     stdev = nanstd(fidTracksZ_denoised(:,1)-fidTracksZ_denoised(:,2));
%     ylim([avg-5*stdev avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,4),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_z']})
%     
%     subplot(4,1,4)
%     euclid_Dist = sqrt(((fidTracksX_denoised(:,1)-fidTracksX_denoised(:,2)).^2)+...
%         ((fidTracksY_denoised(:,1)-fidTracksY_denoised(:,2)).^2+...
%         ((fidTracksZ_denoised(:,1)-fidTracksZ_denoised(:,2)).^2)));
%     plot(euclid_Dist);
%     xlabel('frame number');ylabel('3D shift (nm)');
%     avg = nanmean(euclid_Dist);
%     stdev = nanstd(euclid_Dist);
%     ylim([0 avg+7*stdev])
%     hold on
%     plot(dataSets(2).fidTrack_interpolated_TRE(:,1),'red','LineWidth' ,2)
%     hold off
%     legend({['offset = ' num2str(avg),...
%         ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_3_D']})
%     clear avg
    
    %% Show the fiducial tracks in the same coordinate system
    if numFids == 1 && ~sequentialMeas
    figure_h_a = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    subplot(1,2,1)
    for i = 1:length(dataSets)
        if dataSets(i).fidTransformed
            plot(dataSets(i).fidTrackZ_transformed(:,dataSets(i).fidsToUse),...
                dataSets(i).fidTrackX_transformed(:,dataSets(i).fidsToUse),...
                color{i});
        else
            plot(dataSets(i).fidTrackZ(:,dataSets(i).fidsToUse),...
                dataSets(i).fidTrackX(:,dataSets(i).fidsToUse),...
                color{i});
        end
        
        xlabel('z (nm)');ylabel('x (nm)');
        hold on
        
    end
    hold off
    
    subplot(1,2,2)
    for i = 1:length(dataSets)
        if dataSets(i).fidTransformed
            plot(dataSets(i).fidTrackZ_transformed(:,dataSets(i).fidsToUse),...
                dataSets(i).fidTrackY_transformed(:,dataSets(i).fidsToUse),...
                color{i});
        else
            plot(dataSets(i).fidTrackZ(:,dataSets(i).fidsToUse),...
                dataSets(i).fidTrackY(:,dataSets(i).fidsToUse),...
                color{i});
        end
        
        xlabel('z (nm)');ylabel('y (nm)');
        hold on
        
    end
    hold off
    
    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
    end
    %% prompt to save bead registration figures
    [saveFile, savePath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
    savePath = [savePath saveFile filesep];
    mkdir(savePath);
    
%     saveas(figure_h_a,[savePath '3DFidCorrelation.fig']);
%     saveas(figure_h_a,[savePath '3DFidCorrelation.png']);

    saveas(h_allRegisteredTracks,[savePath 'XYZFidTracks.fig']);
    saveas(h_allRegisteredTracks,[savePath 'XYZFidTracks.png']);

    saveas(h_finalCorrections,[savePath 'XYZFidTracks_avg.fig']);
    saveas(h_finalCorrections,[savePath 'XYZFidTracks_avg.png']);

%     if exist('figure_h_c')
%     saveas(figure_h_c,[savePath 'XYZFidMisregistration.fig']);
%     saveas(figure_h_c,[savePath 'XYZFidMisregistration.png']);
%     close(figure_h_c)
%     end
%     saveas(figure_h_d,[savePath 'XYZFidMisregistration_denoised.fig']);
%     saveas(figure_h_d,[savePath 'XYZFidMisregistration_denoised.png']);
%     close(figure_h_d)
    
elseif ~useFids
        disp('You did not use fiducials!')
        % use the 'driftCorr' name to make the code easier, but not that
        % these are NOT fiducial drift corrected - may want to change this
        % to be easier to parse (e.g. use different variable name and an if
        % statement to define which variable name to use for the following)
        for i = 1:length(dataSets);
            
            if dataSets(i).dataTransformed
                transformedDataSet = i;
                dataSets(i).xLoc_driftCorr = dataSets(i).xLoc_transformed;
                dataSets(i).yLoc_driftCorr = dataSets(i).yLoc_transformed;
                dataSets(i).zLoc_driftCorr = dataSets(i).zLoc_transformed;
            else
                dataSets(i).xLoc_driftCorr = dataSets(i).xLoc;
                dataSets(i).yLoc_driftCorr = dataSets(i).yLoc;
                dataSets(i).zLoc_driftCorr = dataSets(i).zLoc;
                untransformedDataSet = i;
            end
        end
        
        
        [saveFile, savePath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
        savePath = [savePath saveFile filesep];
        mkdir(savePath);
end % end fid or not fid section
        %% Apply index mismatch corrections
    % Todo: This is empirical. A better model accouting for index mismatch
    % needs to developed here.
    for i = 1:length(dataSets)
        dataSets(i).zLoc_driftCorr_indexCorr = dataSets(i).zLoc_driftCorr * nSample/nOil;
    end 
    
    %% Clean up
    tform.FRE = FRE;
    tform.TRE = TRE;
    tform.FRE_full = FRE_full;
    tform.TRE_full = TRE_full;
    tform.matched_cp_reflected = matched_cp_reflected;
    tform.matched_cp_transmitted = matched_cp_transmitted	;
    if exist('matched_cp_transmitted_trans') % is this needed?
        tform.matched_cp_transmitted_trans = matched_cp_transmitted_trans;
    end
    registrationComplete = true
    
    clear FRE TRE FRE_full TRE_full matched_cp_reflected matched_cp_transmitted matched_cp_transmitted_trans
    clear avgDevX avgDevY avgDevZ euclid_Dist goodFits h stdev
    clear globalScale kthNeighbor nControlPoints nCores
    clear figure_h_a figure_h_b figure_h_c figure_h_d
    clear devX devX_denoised devY devY_denoised devZ devZ_denoised
    clear fidTracksX fidTracksX_denoised fidTracksY fidTracksY_denoised fidTracksZ fidTracksZ_denoised
    clear x y z xLocPix yLocPix zLoc_IndexCorrected
    close all

    save([savePath 'registeredSMACMData.mat']);
    
end

%% Display the results the final fused SMACM data
useTimeColors = 0;
frameRange = [1, 100000; 1, 100000];
numPhotonRange = [0 100000];

% [whiteLightFile whiteLightPath] = uigetfile({'*.tif';'*.*'},'Open image stack with white light image');
load(locFiles{untransformedDataSet},'whiteLightFile');
% dlg_title = 'Please Input Parameters';
% prompt = {  'Pixel size (in nm)',...
%     'Size of Points in reconstruction',...
%     'White Light Shift X (in nm)',...
%     'White Light Shift Y (in nm)',...
%     };
% def = {    num2str(nmPerPixel), ...
%     '30', ...
%     num2str(dataSets(untransformedDataSet).wlShift(1)), ...
%     num2str(dataSets(untransformedDataSet).wlShift(2)), ...
%     };
% num_lines = 1;
% inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
% 
% nmPerPixel = str2double(inputdialog{1});
% scatterSize = str2double(inputdialog{2});
% wlShiftX = str2double(inputdialog{3});
% wlShiftY = str2double(inputdialog{4});
scatterSize = 30;
wlShiftX = num2str(dataSets(untransformedDataSet).wlShift(1));
wlShiftY = num2str(dataSets(untransformedDataSet).wlShift(2));

% pass = 1;
anotherpass = true;

while anotherpass == true
    close all
    
    
    %% Chose a desired parameter set for reconstruction
    dlg_title = 'Please Input Parameters';
    prompt = {  'Size of points in reconstruction',...
        'Temporal Color Coding',...
        'White light shift X (in nm)',...
        'White light shift Y (in nm)',...
        'First frame (first channel)',...
        'Last frame (first channel)',...
        'First frame (second channel)',...
        'Last frame (second channel)',...
        'Number of photons lower bound',...
        'Number of photons upper bound',...
        };
    def = { ...
        num2str(scatterSize), ...
        num2str(useTimeColors), ...
        num2str(wlShiftX), ...
        num2str(wlShiftY), ...
        num2str(frameRange(1,1)), ...
        num2str(frameRange(1,2)), ...
        num2str(frameRange(2,1)), ...
        num2str(frameRange(2,2)), ...
        num2str(numPhotonRange(1)), ...
        num2str(numPhotonRange(2)), ...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    scatterSize = str2double(inputdialog{1});
    useTimeColors = str2double(inputdialog{2});
    wlShiftX = str2double(inputdialog{3});
    wlShiftY = str2double(inputdialog{4});
    frameRange = [str2double(inputdialog{5}) str2double(inputdialog{6});...
                  str2double(inputdialog{7}) str2double(inputdialog{8})];
    numPhotonRange = [str2double(inputdialog{9}) str2double(inputdialog{10})];
    
     %% Plot the white light image if specified
        if ~any((whiteLightFile==0))
            if exist('whiteLightPath') && ~any(strfind(whiteLightFile,whiteLightPath)) % if whiteLightFile not already complete path
                whiteLightFile = [whiteLightPath whiteLightFile];
            end
            whiteLightInfo = imfinfo(whiteLightFile);
            whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
%             im = tiffread2(whiteLightFile, 1, length(whiteLightInfo));
            % average white light images together to get better SNR
            for a = 1:length(whiteLightInfo)
                whiteLight = whiteLight + double(imread(whiteLightFile,a, ...
                    'Info', whiteLightInfo));
%                 whiteLight = whiteLight + double(im(a).data);
            end
%             clear im
            % resize white light to the size of the ROI of the single molecule fits
            if ~exist('tformChan') || strcmp(tformChan,'r') || strcmp(tformChan,'R')
                ROI_initial = [1, 1, min(350,size(whiteLight,1)-1), min(350,size(whiteLight,2)-1)];
            elseif strcmp(tformChan,'y') || strcmp(tformChan,'Y')
                ROI_initial = [200, 200, size(whiteLight,1)-200, size(whiteLight,2)-200];
            end
            whiteLight = whiteLight(ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1,ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1);
            % rescale white light image to vary from 0 to 1
            %         whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
            borderStart = 50; borderEnd = 230;
            wlCenter = whiteLight(borderStart:borderEnd,...
                borderStart:borderEnd);
            whiteLight = (whiteLight-min(wlCenter(:)))/...
                (max(wlCenter(:))-min(wlCenter(:)));
            whiteLight(whiteLight(:)>1) = 1; whiteLight(whiteLight(:)<0) = 0;
            
            
            [xWL yWL] = meshgrid((ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1) * nmPerPixel + wlShiftX, ...
                (ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1) * nmPerPixel + wlShiftY);
        end
    
    %% ask user what region to plot in superresolution image
    
%     if pass ~= 1
%         if whiteLightFile ~= 0
%             whiteLightInfo = imfinfo([whiteLightPath whiteLightFile]);
%             whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
%             % average white light images together to get better SNR
%             for a = 1:length(whiteLightInfo)
%                 whiteLight = whiteLight + double(imread([whiteLightPath whiteLightFile],a, ...
%                     'Info', whiteLightInfo));
%             end
%             % resize white light to the size of the ROI of the single molecule fits
%             ROI_initial = [1, 1, 270, 270];
%             whiteLight = whiteLight(ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1,ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1);
%             % rescale white light image to vary from 0 to 1
%             whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
%             [xWL yWL] = meshgrid((ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1) * nmPerPixel + wlShiftX, ...
%                 (ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1) * nmPerPixel + wlShiftY);
%         end
%     end
    
    figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    
    if whiteLightFile ~= 0
        xRange = xWL(1,:);
        yRange = yWL(:,1);
        % pick region that contains background
        imagesc(xRange,yRange,whiteLight);axis image;colormap gray;
    else
        xRange = [min(vertcat(dataSets.xLoc_driftCorr)) max(vertcat(dataSets.xLoc_driftCorr))];
        yRange = [min(vertcat(dataSets.yLoc_driftCorr)) max(vertcat(dataSets.yLoc_driftCorr))];
        [xBl, yBl] = meshgrid(round(xRange(1)):100:round(xRange(2)),...
                              round(yRange(1)):100:round(yRange(2)));
        imagesc(yBl(:,1),xBl(1,:),zeros(size(xBl)),[-1 0]); axis image; colormap gray;
    end
    
    
    hold on;
    
    for i = 1:length(dataSets)
        
        scatter(dataSets(i).xLoc_driftCorr,dataSets(i).yLoc_driftCorr,1,'filled', color{i});
        xlim([min(dataSets(i).xLoc_driftCorr) max(dataSets(i).xLoc_driftCorr)]);
        ylim([min(dataSets(i).yLoc_driftCorr) max(dataSets(i).yLoc_driftCorr)]);
        xlabel('x (nm)');ylabel('y (nm)');
        axis ij;
        hold on
        
    end
    
    [ROI, xi, yi] = roipoly;
    plot(xi, yi, 'Color','black', 'LineWidth',2)
    hold off
    
    %% filter out localizations outside of ROI
    
    croppedDataSets = [];
    interpolated_FREs = [];
    interpolated_TREs = [];
    plotRange = [];
    f = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','k','renderer','opengl', 'Toolbar', 'figure');
    if whiteLightFile~=0
        %imagesc(xRange,yRange,whiteLight);axis image;colormap gray;hold on;
        [x,y,z] = meshgrid(xRange,yRange,[-2000 2000]);
        xslice = []; yslice = []; zslice = -600;
        h=slice(x,y,z,repmat(whiteLight,[1 1 2]),xslice,yslice,zslice,'nearest');
        set(h,'EdgeColor','none','FaceAlpha',0.75);
        colormap gray; 
    end
    
    hold on; grid on;
    
    for i = 1:length(dataSets)
        
        xLoc = dataSets(i).xLoc_driftCorr;
        yLoc = dataSets(i).yLoc_driftCorr;
        zLoc = dataSets(i).zLoc_driftCorr;
        zLoc_indexCorr = dataSets(i).zLoc_driftCorr_indexCorr;
        frameNum = dataSets(i).frameNum;
        
        validPoints = inpolygon(xLoc,yLoc,xi, yi);
        validPoints = validPoints & frameNum >= frameRange(i,1) & frameNum <= frameRange(i,2);
        xLoc = xLoc(validPoints);
        yLoc = yLoc(validPoints);
        zLoc = zLoc(validPoints);
        zLoc_indexCorr = zLoc_indexCorr(validPoints);
        
        %% Assemble the croppedDataSet structure
        
        croppedDataSet.frameNum =  dataSets(i).frameNum(validPoints);
        croppedDataSet.xLoc = xLoc;
        croppedDataSet.yLoc = yLoc;
        croppedDataSet.zLoc = zLoc;
        croppedDataSet.zLoc_indexCorr = zLoc_indexCorr;
        croppedDataSet.sigmaX = dataSets(i).sigmaX(validPoints);
        croppedDataSet.sigmaY = dataSets(i).sigmaY(validPoints);
        croppedDataSet.sigmaZ = dataSets(i).sigmaZ(validPoints);
        croppedDataSet.numPhotons = dataSets(i).numPhotons(validPoints);
        croppedDataSet.meanBkgnd = dataSets(i).meanBkgnd(validPoints);
        if dataSets(i).dataTransformed == 1
            if ~exist('F_FRE')&&exist('tform')
                F_FRE = tform.interpolationObjects.F_FRE;
                F_TRE = tform.interpolationObjects.F_TRE;
            end
            interpolated_FRE = F_FRE(xLoc,yLoc,zLoc);
            interpolated_TRE = F_TRE(xLoc,yLoc,zLoc);
            croppedDataSet.interpolated_FRE = interpolated_FRE;
            croppedDataSet.interpolated_TRE = interpolated_TRE;
            interpolated_FREs = [interpolated_FREs; interpolated_FRE];
            interpolated_TREs = [interpolated_TREs; interpolated_TRE];
        else
            croppedDataSet.interpolated_FRE = nan;
            croppedDataSet.interpolated_TRE = nan;
        end
        
        %% Display the results in 3D
        % plot 3D scatterplot of localizations with white light
        scatter3(xLoc,yLoc,zLoc,scatterSize,'filled',color{i});
        axis vis3d equal;
        
        plotRange = [plotRange; [min(xLoc) max(xLoc) min(yLoc) max(yLoc) min(zLoc) max(zLoc)]];
        plotRange = max(plotRange,[],1);
        
        xlim([plotRange(1) plotRange(2)]);
        ylim([plotRange(3) plotRange(4)]);
        zlim([plotRange(5) plotRange(6)]);
        xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
        title({[num2str(length(interpolated_FREs)) ' transformed localizations'];...
            ['Mean FRE = ' num2str(nanmean(interpolated_FREs)) ' nm'];...
            ['Mean TRE = ' num2str(nanmean(interpolated_TREs)) ' nm'];...
            'no index correction applied'},...
            'color','w');
        set(gca,'color','k');
        set(gca,'xcolor','w');set(gca,'ycolor','w');set(gca,'zcolor','w');
        
        %        %% Display the results in 3D with index correction
        %         f_2 = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','k','renderer','opengl', 'Toolbar', 'figure');
        %         if whiteLightFile~=0
        %             %imagesc(xRange,yRange,whiteLight);axis image;colormap gray;hold on;
        %             [x,y,z] = meshgrid(xRange,yRange,[-2000 2000]);
        %             xslice = []; yslice = []; zslice = -600;
        %             h=slice(x,y,z,repmat(whiteLight,[1 1 2]),xslice,yslice,zslice,'nearest');
        %             set(h,'EdgeColor','none','FaceAlpha',0.75);
        %             colormap gray; hold on;
        %         end
        %
        %         % plot 3D scatterplot of localizations with white light
        %         scatter3(xLoc,yLoc,zLoc_indexCorr,scatterSize,'filled',color{i});
        %         axis vis3d equal;
        %
        %         xlim([plotRange(1) plotRange(2)]);
        %         ylim([plotRange(3) plotRange(4)]);
        %         zlim([plotRange(5) plotRange(6)]);
        %         xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
        %         title({[num2str(length(interpolated_FREs)) ' transformed localizations'];...
        %             ['Mean FRE = ' num2str(mean(interpolated_FREs)) ' nm'];...
        %             ['Mean TRE = ' num2str(mean(interpolated_TREs)) ' nm'];...
        %             'index correction applied'},...
        %             'color','w');
        %         set(gca,'color','k');
        %         set(gca,'xcolor','w');set(gca,'ycolor','w');set(gca,'zcolor','w');
        
        croppedDataSets = [croppedDataSets, croppedDataSet];
        
    end
    hold off
    
    %% Construct a questdlg with three options
    
    % f = figure;
    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    % disp('This will print immediately');
    uiwait(gcf);
    % disp('This will print after you click Continue');
    %     close(f);
    
    dlg_title = 'Replot';
    prompt = {'Would you like to replot with a different parameter set?'};
    def =       { 'Yes'  };
    questiondialog = questdlg(prompt,dlg_title, def);
    % Handle response
    switch questiondialog
        case 'Yes'
%             pass = pass + 1;
        case 'No'
            anotherpass = false;
        case 'Cancel'
            error('User cancelled the program');
    end
    
end

%% clean up
if exist('F_FRE_X') % this is a kludge to allow reopening old outputs from this function
tform.interpolationObjects.F_FRE = F_FRE;
tform.interpolationObjects.F_FRE_X = F_FRE_X;
tform.interpolationObjects.F_FRE_Y = F_FRE_Y;
tform.interpolationObjects.F_FRE_Z = F_FRE_Z;
tform.interpolationObjects.F_TRE = F_TRE;
tform.interpolationObjects.F_TRE_X = F_TRE_X;
tform.interpolationObjects.F_TRE_Y = F_TRE_Y;
tform.interpolationObjects.F_TRE_Z = F_TRE_Z;
end
clear F_FRE F_FRE_X F_FRE_Y F_FRE_Z F_TRE F_TRE_X F_TRE_Y F_TRE_Z
clear a anotherpass bead croppedDataSet def dlg_title f h i
clear interpolated_FRE interpolated_FREs interpolated_TRE interpolated_TREs
clear num_lines pass prompt questiondialog
clear xLoc xslice yLoc yslice zLoc zLoc_indexCorr zslice
clear x xRange xWL xi y yRange yWL yi z

%% prompt to save data and figures
[newFile, newPath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
if ~isequal(newFile,0)
% newPath = [newPath filesep];
mkdir(newPath);
save([newPath filesep newFile '_multicolorSMACM.mat']);
else
    newPath = savePath;
    newFile = 'output';
    save([savePath filesep newFile 'MulticolorSMACM.mat']);
end

saveas(gcf,[newPath newFile 'MulticolorSMACM_3D.fig']);
close
saveas(gcf,[newPath newFile 'MulticolorSMACM_2D.fig']);
close all

end