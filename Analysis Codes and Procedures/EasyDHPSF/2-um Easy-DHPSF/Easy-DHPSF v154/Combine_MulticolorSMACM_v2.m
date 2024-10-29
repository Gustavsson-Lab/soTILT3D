function [] = Combine_MulticolorSMACM()
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
    
    [tformFile tformPath] = uigetfile({'*.mat';'*.*'},'Open 3D_Transform.mat');
    if isequal(tformFile,0)
        error('User cancelled the program');
    end
    load([tformPath tformFile]);

    %% Prepare for data interpolation
    x = matched_cp_reflected(:,5);
    y = matched_cp_reflected(:,6);
    z = matched_cp_reflected(:,7);
    F_FRE = TriScatteredInterp(x,y,z,FRE_full(:,1), 'natural');
    F_FRE_X = TriScatteredInterp(x,y,z,FRE_full(:,2), 'natural');
    F_FRE_Y = TriScatteredInterp(x,y,z,FRE_full(:,3), 'natural');
    F_FRE_Z = TriScatteredInterp(x,y,z,FRE_full(:,4), 'natural');
    F_TRE = TriScatteredInterp(x,y,z,TRE_full(:,1), 'natural');
    F_TRE_X = TriScatteredInterp(x,y,z,TRE_full(:,2), 'natural');
    F_TRE_Y = TriScatteredInterp(x,y,z,TRE_full(:,3), 'natural');
    F_TRE_Z = TriScatteredInterp(x,y,z,TRE_full(:,4), 'natural');
    
    [whiteLightFile whiteLightPath] = uigetfile({'*.tif';'*.*'},'Open image stack with white light image');
    
    [LocFile LocPath] = uigetfile({'*.mat';'*.*'},'Open data file file #1 with filtered molecule fits');
    if isequal(LocFile,0)
        error('User cancelled the program');
    end
    load([LocPath LocFile]);
    
    fileNum = 1;
    LocFiles = {};
    dataSets = [];
    
    while ~isequal(LocFile,0)
        
        LocFiles = [LocFiles; {[LocPath LocFile]}];
        % load data
        load([LocPath LocFile]);
        
        %    Assemble the structure for this dataset
        dataSet.frameNum = frameNum;
        dataSet.xLoc = xLoc;
        dataSet.yLoc = yLoc;
        dataSet.zLoc = zLoc;
        dataSet.sigmaX = sigmaX;
        dataSet.sigmaY = sigmaY;
        dataSet.sigmaZ = sigmaZ;
        dataSet.numPhotons = numPhotons;
        dataSet.meanBkgnd = meanBkgnd;
        dataSet.fidTrackX = fidTrackX;
        dataSet.fidTrackY = fidTrackY;
        dataSet.fidTrackZ = fidTrackZ;
        dataSet.wlShift = [wlShiftX, wlShiftY];
        
        dlg_title = 'Transform dataset';
        prompt = {'Do you want to transform this dataset?'};
        def =       { 'Yes' };
        questiondialog = questdlg(prompt,dlg_title, def);
        % Handle response
        switch questiondialog
            case 'Yes'
                
                % transform SM data
                dataSet.transformedDataset = true;
                transformedData = transformData([xLoc, yLoc, zLoc],tform);
                dataSet.xLoc_transformed = transformedData(:,1);
                dataSet.yLoc_transformed = transformedData(:,2);
                dataSet.zLoc_transformed = transformedData(:,3);
                
                % transform fiducial data
                % ToDo:  Code in the possibility for multiple fiducial in the
                % field-of-view.  Currently the code only handles a single one.
                dataSet.fidTrackX_transformed = NaN(length(fidTrackX),1);
                dataSet.fidTrackY_transformed = NaN(length(fidTrackX),1);
                dataSet.fidTrackZ_transformed = NaN(length(fidTrackX),1);
                transformedData = transformData([fidTrackX(~isnan(fidTrackX)),...
                    fidTrackY(~isnan(fidTrackY)),...
                    fidTrackZ(~isnan(fidTrackZ))],tform);
                dataSet.fidTrackX_transformed(~isnan(fidTrackX)) = transformedData(:,1);
                dataSet.fidTrackY_transformed(~isnan(fidTrackY)) = transformedData(:,2);
                dataSet.fidTrackZ_transformed(~isnan(fidTrackZ)) = transformedData(:,3);
                
                dataSet.fidTrack_interpolated_FRE = NaN(length(fidTrackX),4);
                dataSet.fidTrack_interpolated_TRE = NaN(length(fidTrackX),4);
                dataSet.fidTrack_interpolated_FRE(~isnan(fidTrackX),:) = [...
                    F_FRE(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                    F_FRE_X(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                    F_FRE_Y(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                    F_FRE_Z(transformedData(:,1),transformedData(:,2),transformedData(:,3)) ];
                dataSet.fidTrack_interpolated_TRE(~isnan(fidTrackX),:) = [...
                    F_TRE(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                    F_TRE_X(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                    F_TRE_Y(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                    F_TRE_Z(transformedData(:,1),transformedData(:,2),transformedData(:,3)) ];
                
            case 'No'
                
                dataSet.transformedDataset = false;
                dataSet.xLoc_transformed = NaN;
                dataSet.yLoc_transformed = NaN;
                dataSet.zLoc_transformed = NaN;
                dataSet.fidTrackX_transformed =NaN;
                dataSet.fidTrackY_transformed = NaN;
                dataSet.fidTrackZ_transformed = NaN;
                dataSet.fidTrack_interpolated_FRE = NaN;
                dataSet.fidTrack_interpolated_TRE = NaN;
                
                %             dataSet.xLoc_transformed = NaN(length(xLoc),1);
                %             dataSet.yLoc_transformed = NaN(length(xLoc),1);
                %             dataSet.zLoc_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrackX_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrackY_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrackZ_transformed = NaN(length(xLoc),1);
                %             dataSet.fidTrack_interpolated_FRE = NaN(length(fidTrackX),4);
                %             dataSet.fidTrack_interpolated_TRE = NaN(length(fidTrackX),4);
                
            case 'Cancel'
                error('User cancelled the program');
        end
        
        dataSet.LocFile = LocFile;
        dataSet.LocPath = LocPath;
        
        dataSet.frameRange = frameRange;
        dataSet.zRange = zRange;
        dataSet.ampRatioLimit = ampRatioLimit;
        dataSet.fitErrorRange = fitErrorRange;
        dataSet.lobeDistBounds = lobeDistBounds;
        dataSet.sigmaBounds = sigmaBounds;
        dataSet.sigmaRatioLimit = sigmaRatioLimit;
        dataSet.numPhotonRange = numPhotonRange;
        
        
        % Append the structure to the previous structures in the array
        dataSets = [dataSets, dataSet];
        
        fileNum = fileNum+1;
        [LocFile LocPath] = uigetfile({'*.mat';'*.*'},...
            ['Open data file #' num2str(fileNum) ' with filtered molecule fits']);
        
    end
    
    clear LocFile LocPath ampRatioLimit dataSet def dlg_title fileNum
    clear fitErrorRange frameNum frameRange lobeDistBounds meanBkgnd
    clear numPhotonRange numPhotons prompt questiondialog sigmaBounds
    clear sigmaRatioLimit sigmaX sigmaY sigmaZ transformedData
    clear xLoc yLoc zLoc zRange fidTrackX fidTrackY fidTrackZ
    
    %% Identify frames based on the sequence log
    % load data and register sequence log to data frames
    
    [logFile logPath] = uigetfile({'*.dat'},'Open shutter sequence log file');
    if isequal(logFile,0)
        error('User cancelled the program');
    end

    sifLogData =  importdata([logPath logFile]);
    frames_green = find(sifLogData(:,2) == 1);
    frames_red = find(sifLogData(:,3) == 1);
    selectedFrames = unique(sort([frames_green; frames_red]));
    % find shutter transition frames
    %     endFrames_green = frames_green(find(diff(frames_green)==36));
    temp = diff(frames_red);
    endFrames_red = frames_red(find(diff(frames_red)== mean(temp(temp>4))));
    startFrames_green = endFrames_red+1;
    clear temp

    %% Show the decomposed fiducial tracks
    
    fidTracksX = [];
    fidTracksY = [];
    fidTracksZ = [];
    
    for i = 1:length(dataSets)
        if dataSets(i).transformedDataset
            transformedDataSet = i;
            fidTracksX = [fidTracksX, dataSets(i).fidTrackX_transformed];
            fidTracksY = [fidTracksY, dataSets(i).fidTrackY_transformed];
            fidTracksZ = [fidTracksZ, dataSets(i).fidTrackZ_transformed];
        else
            untransformedDataSet = i;
            fidTracksX = [fidTracksX, dataSets(i).fidTrackX];
            fidTracksY = [fidTracksY, dataSets(i).fidTrackY];
            fidTracksZ = [fidTracksZ, dataSets(i).fidTrackZ];
        end
    end
    
    color = {'Green', 'Red', 'Blue', 'Cyan', 'Yellow' };
    figure_h_a = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    
    for i = 1:length(dataSets)
        if dataSets(i).transformedDataset
            subplot(3,1,1)
            plot(dataSets(i).fidTrackX_transformed, color{i});
            xlabel('frame number');ylabel('X position (nm)');
            title('Raw red and green channel fiducial tracks after registration')
            hold on
            subplot(3,1,2)
            plot(dataSets(i).fidTrackY_transformed, color{i});
            xlabel('frame number');ylabel('Y position (nm)');
            hold on
            subplot(3,1,3)
            plot(dataSets(i).fidTrackZ_transformed, color{i});
            xlabel('frame number');ylabel('Z position (nm)');
            hold on
        else
            subplot(3,1,1)
            plot(dataSets(i).fidTrackX, color{i});
            xlabel('frame number');ylabel('X position (nm)');
            hold on
            subplot(3,1,2)
            plot(dataSets(i).fidTrackY, color{i});
            xlabel('frame number');ylabel('Y position (nm)');
            hold on
            subplot(3,1,3)
            plot(dataSets(i).fidTrackZ, color{i});
            xlabel('frame number');ylabel('Z position (nm)');
            hold on
        end
        hold on
    end
    hold off
%     % display the mean track also
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

    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
        
    %% Pick the appropriate track to use as drift correction
    
    dlg_title = 'Inspect Fiducial Tracks';
    prompt = {'Choose Fiducial Track to denoise and use to correct drift'};
    questiondialog = questdlg(prompt,dlg_title, 'untransformed', 'transformed', 'Cancel', 'untransformed'  );
    % Handle response
    switch questiondialog
        case 'untransformed'
                chosenFidTrack = untransformedDataSet;
        case 'transformed'
                chosenFidTrack = transformedDataSet;
        case 'Cancel'
            error('User cancelled the program');
    end
    
    
    %% Denoise the raw data   
   
    dlg_title = 'Denoise the tracks';
    prompt = {'Try to get the smoothest line possible that still represents the drift'};
    questiondialog = questdlg(prompt,dlg_title, 'Continue', 'Cancel', 'Continue'  );
    % Handle response
    switch questiondialog
        case 'Continue'
        case 'Cancel'
            error('User cancelled the program');
    end
    
    % Denoise the red illuminated frames 
    fidTracksX_denoised_redFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
    fidTracksY_denoised_redFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
    fidTracksZ_denoised_redFrames = nan(size(fidTracksX,1),size(fidTracksX,2));

    [fidTracksX_denoised_redFrames(frames_red,transformedDataSet),...
        fidTracksY_denoised_redFrames(frames_red,transformedDataSet),...
        fidTracksZ_denoised_redFrames(frames_red,transformedDataSet)] = f_waveletFidTracks(...
        fidTracksX(frames_red,transformedDataSet),...
        fidTracksY(frames_red,transformedDataSet),...
        fidTracksZ(frames_red,transformedDataSet),0);
    [fidTracksX_denoised_redFrames(frames_red,untransformedDataSet),...
        fidTracksY_denoised_redFrames(frames_red,untransformedDataSet),...
        fidTracksZ_denoised_redFrames(frames_red,untransformedDataSet)] = f_waveletFidTracks(...
        fidTracksX(frames_red,untransformedDataSet),...
        fidTracksY(frames_red,untransformedDataSet),...
        fidTracksZ(frames_red,untransformedDataSet),0);
    
    % Denoise the green illuminated frames 
    fidTracksX_denoised_greenFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
    fidTracksY_denoised_greenFrames = nan(size(fidTracksX,1),size(fidTracksX,2));
    fidTracksZ_denoised_greenFrames = nan(size(fidTracksX,1),size(fidTracksX,2));

    [fidTracksX_denoised_greenFrames(frames_green,transformedDataSet),...
        fidTracksY_denoised_greenFrames(frames_green,transformedDataSet),...
        fidTracksZ_denoised_greenFrames(frames_green,transformedDataSet)] = f_waveletFidTracks(...
        fidTracksX(frames_green,transformedDataSet),...
        fidTracksY(frames_green,transformedDataSet),...
        fidTracksZ(frames_green,transformedDataSet),0);
    [fidTracksX_denoised_greenFrames(frames_green,untransformedDataSet),...
        fidTracksY_denoised_greenFrames(frames_green,untransformedDataSet),...
        fidTracksZ_denoised_greenFrames(frames_green,untransformedDataSet)] = f_waveletFidTracks(...
        fidTracksX(frames_green,untransformedDataSet),...
        fidTracksY(frames_green,untransformedDataSet),...
        fidTracksZ(frames_green,untransformedDataSet),0);

    % Show the piecewise denoised tracks of the chosen channel
    color = {'Green', 'Red', 'Blue', 'Cyan', 'Yellow' };
    figure_h_b = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    subplot(3,1,1)
    plot(fidTracksX_denoised_greenFrames(:,chosenFidTrack), color{1});
    hold on
    plot(fidTracksX_denoised_redFrames(:,chosenFidTrack), color{2});
    hold off
    xlabel('frame number');ylabel('X position (nm)');
    title('Denoised red and green illuminated chosen fiducial track')
    
    subplot(3,1,2)
    plot(fidTracksY_denoised_greenFrames(:,chosenFidTrack), color{1});
    hold on
    plot(fidTracksY_denoised_redFrames(:,chosenFidTrack), color{2});
    hold off
    xlabel('frame number');ylabel('Y position (nm)');

    subplot(3,1,3)
    plot(fidTracksZ_denoised_greenFrames(:,chosenFidTrack), color{1});
    hold on
    plot(fidTracksZ_denoised_redFrames(:,chosenFidTrack), color{2});
    hold off
    xlabel('frame number');ylabel('Z position (nm)');
    
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
    
    %% Compute the shifts at the transition Frames
    shiftX =  fidTracksX_denoised_greenFrames(startFrames_green,:) - ...
        fidTracksX_denoised_redFrames(endFrames_red,:);
    shiftY =  fidTracksY_denoised_greenFrames(startFrames_green,:) - ...
        fidTracksY_denoised_redFrames(endFrames_red,:);
    shiftZ =  fidTracksZ_denoised_greenFrames(startFrames_green,:) - ...
        fidTracksZ_denoised_redFrames(endFrames_red,:);
    
    shifts_mean = [nanmean(shiftX,1);nanmean(shiftY,1);nanmean(shiftZ,1)]
    shifts_std = [nanstd(shiftX,1);nanstd(shiftY,1);nanstd(shiftZ,1)]
    
    lb = shifts_mean-1.5*shifts_std;
    ub = shifts_mean+1.5*shifts_std;
    
    shifts_mean_filt = [[mean(shiftX(shiftX(:,1)>=lb(1,1) & shiftX(:,1)<=ub(1,1),1)),...
        mean(shiftX(shiftX(:,2)>=lb(1,2) & shiftX(:,2)<=ub(1,2),2))];...
        [mean(shiftY(shiftY(:,1)>=lb(2,1) & shiftY(:,1)<=ub(2,1),1)),...
        mean(shiftY(shiftY(:,2)>=lb(2,2) & shiftY(:,2)<=ub(2,2),2))];...
        [mean(shiftZ(shiftZ(:,1)>=lb(3,1) & shiftZ(:,1)<=ub(3,1),1)),...
        mean(shiftZ(shiftZ(:,2)>=lb(3,2) & shiftZ(:,2)<=ub(3,2),2))]];
    
    shifts_std_filt = [[std(shiftX(shiftX(:,1)>=lb(1,1) & shiftX(:,1)<=ub(1,1),1)),...
        std(shiftX(shiftX(:,2)>=lb(1,2) & shiftX(:,2)<=ub(1,2),2))];...
        [std(shiftY(shiftY(:,1)>=lb(2,1) & shiftY(:,1)<=ub(2,1),1)),...
        std(shiftY(shiftY(:,2)>=lb(2,2) & shiftY(:,2)<=ub(2,2),2))];...
        [std(shiftZ(shiftZ(:,1)>=lb(3,1) & shiftZ(:,1)<=ub(3,1),1)),...
        std(shiftZ(shiftZ(:,2)>=lb(3,2) & shiftZ(:,2)<=ub(3,2),2))]];
    
    figure_h_c = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    subplot(3,1,1)
    hist(fidTracksX_denoised_greenFrames(startFrames_green,:) - fidTracksX_denoised_redFrames(endFrames_red,:),-150:5:150)
    xlim([-160 160])
    legend({'green channel';'red channel'})
    title('Shift distributions due to changing illumination')
    subplot(3,1,2)
    hist(fidTracksY_denoised_greenFrames(startFrames_green,:) - fidTracksY_denoised_redFrames(endFrames_red,:),-150:5:150)
    xlim([-160 160])
    legend({'green channel';'red channel'})
    subplot(3,1,3)
    hist(fidTracksZ_denoised_greenFrames(startFrames_green,:) - fidTracksZ_denoised_redFrames(endFrames_red,:),-150:5:150)
    xlim([-160 160])
    legend({'green channel';'red channel'})
    
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
    
    %% shift the green illuminated portion of the fiducial tracks to fit the red illuminated portion 

    fidTracksX_shifted = fidTracksX;
    fidTracksY_shifted = fidTracksY;
    fidTracksZ_shifted = fidTracksZ;

    fidTracksX_shifted(frames_green,1) = fidTracksX(frames_green,1) - shifts_mean_filt(1,1);
    fidTracksX_shifted(frames_green,2) = fidTracksX(frames_green,2) - shifts_mean_filt(1,2);
    fidTracksY_shifted(frames_green,1) = fidTracksY(frames_green,1) - shifts_mean_filt(2,1);
    fidTracksY_shifted(frames_green,2) = fidTracksY(frames_green,2) - shifts_mean_filt(2,2);
    fidTracksZ_shifted(frames_green,1) = fidTracksZ(frames_green,1) - shifts_mean_filt(3,1);
    fidTracksZ_shifted(frames_green,2) = fidTracksZ(frames_green,2) - shifts_mean_filt(3,2);
    
    % Show the illumination induced shifted fiducial tracks to make sure the offsets are taken care of 
    color = {'Green', 'Red', 'Blue', 'Cyan', 'Yellow' };
    figure_h_d = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    for i = 1:length(dataSets)
        subplot(3,1,1)
        plot(fidTracksX_shifted(:,i), color{i});
        xlabel('frame number');ylabel('X position (nm)');
        title('Corrected fiducial tracks after illumination induced shift correction')
        hold on
        subplot(3,1,2)
        plot(fidTracksY_shifted(:,i), color{i});
        xlabel('frame number');ylabel('Y position (nm)');
        hold on
        subplot(3,1,3)
        plot(fidTracksZ_shifted(:,i), color{i});
        xlabel('frame number');ylabel('Z position (nm)');
        hold on
    end
    hold off
    
    %% denoise the corrected tracks

        dlg_title = 'Denoise the tracks';
    prompt = {'Try to get the smoothest line possible that still represents the drift'};
    questiondialog = questdlg(prompt,dlg_title, 'Continue', 'Cancel', 'Continue'  );
    % Handle response
    switch questiondialog
        case 'Continue'
        case 'Cancel'
            error('User cancelled the program');
    end
    
    [fidTracksX_shifted_denoised(:,transformedDataSet),...
        fidTracksY_shifted_denoised(:,transformedDataSet),...
        fidTracksZ_shifted_denoised(:,transformedDataSet)] = f_waveletFidTracks(...
        fidTracksX_shifted(:,transformedDataSet),...
        fidTracksY_shifted(:,transformedDataSet),...
        fidTracksZ_shifted(:,transformedDataSet),0);
    [fidTracksX_shifted_denoised(:,untransformedDataSet),...
        fidTracksY_shifted_denoised(:,untransformedDataSet),...
        fidTracksZ_shifted_denoised(:,untransformedDataSet)] = f_waveletFidTracks(...
        fidTracksX_shifted(:,untransformedDataSet),...
        fidTracksY_shifted(:,untransformedDataSet),...
        fidTracksZ_shifted(:,untransformedDataSet),0);
    
        % Show the chromatic shifted fiducial tracks and the denoising result 
    color = {'Green', 'Red', 'Blue', 'Cyan', 'Yellow' };
    figure_h_e = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    for i = 1:length(dataSets)
        subplot(3,1,1)
        plot(fidTracksX_shifted(:,i), color{i});
        xlabel('frame number');ylabel('X position (nm)');
        title('Corrected fiducial tracks and denoising result')
        hold on
        subplot(3,1,2)
        plot(fidTracksY_shifted(:,i), color{i});
        xlabel('frame number');ylabel('Y position (nm)');
        hold on
        subplot(3,1,3)
        plot(fidTracksZ_shifted(:,i), color{i});
        xlabel('frame number');ylabel('Z position (nm)');
        hold on
    end
        for i = 1:length(dataSets)
        subplot(3,1,1)
        plot(fidTracksX_shifted_denoised(:,i), color{i});
        xlabel('frame number');ylabel('X position (nm)');
        hold on
        subplot(3,1,2)
        plot(fidTracksY_shifted_denoised(:,i), color{i});
        xlabel('frame number');ylabel('Y position (nm)');
        hold on
        subplot(3,1,3)
        plot(fidTracksZ_shifted_denoised(:,i), color{i});
        xlabel('frame number');ylabel('Z position (nm)');
        hold on
    end
    hold off
    
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);

    %% Center the last frames of the denoised fiducial track

    syncFrames = find(~isnan(fidTracksX(:,1)));
    syncFrames = syncFrames(end-(numSyncFrames-1):end);

    avgDevX(:,1) = fidTracksX_shifted_denoised(:,1) - mean(fidTracksX_shifted_denoised(syncFrames,1),1);
    avgDevY(:,1) = fidTracksY_shifted_denoised(:,1) - mean(fidTracksY_shifted_denoised(syncFrames,1),1);
    avgDevZ(:,1) = fidTracksZ_shifted_denoised(:,1) - mean(fidTracksZ_shifted_denoised(syncFrames,1),1);
    
    avgDevX(:,2) = fidTracksX_shifted_denoised(:,2) - mean(fidTracksX_shifted_denoised(syncFrames,2),1);
    avgDevY(:,2) = fidTracksY_shifted_denoised(:,2) - mean(fidTracksY_shifted_denoised(syncFrames,2),1);
    avgDevZ(:,2) = fidTracksZ_shifted_denoised(:,2) - mean(fidTracksZ_shifted_denoised(syncFrames,2),1);
    
        % Show the drift correction to be applied
    color = {'Green', 'Red', 'Blue', 'Cyan', 'Yellow' };
    figure_h_f = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);

        subplot(3,1,1)
        plot(avgDevX(:,chosenFidTrack), color{chosenFidTrack});
        xlabel('frame number');ylabel('X position (nm)');
        title('Applied drift correction')

        subplot(3,1,2)
        plot(avgDevY(:,chosenFidTrack), color{chosenFidTrack});
        xlabel('frame number');ylabel('Y position (nm)');

        subplot(3,1,3)
        plot(avgDevZ(:,chosenFidTrack), color{chosenFidTrack});
        xlabel('frame number');ylabel('Z position (nm)');
        
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);

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

        if dataSets(i).transformedDataset
            dataSets(i).xLoc_driftCorr = dataSets(i).xLoc_transformed - avgDevX(dataSets(i).frameNum,chosenFidTrack);
            dataSets(i).yLoc_driftCorr = dataSets(i).yLoc_transformed - avgDevY(dataSets(i).frameNum,chosenFidTrack);
            dataSets(i).zLoc_driftCorr = dataSets(i).zLoc_transformed - avgDevZ(dataSets(i).frameNum,chosenFidTrack);
        else
            dataSets(i).xLoc_driftCorr = dataSets(i).xLoc - avgDevX(dataSets(i).frameNum,chosenFidTrack);
            dataSets(i).yLoc_driftCorr = dataSets(i).yLoc - avgDevY(dataSets(i).frameNum,chosenFidTrack);
            dataSets(i).zLoc_driftCorr = dataSets(i).zLoc - avgDevZ(dataSets(i).frameNum,chosenFidTrack);
        end
    end
    
        %% Apply index mismatch corrections
    % Todo: This is empirical. A better model accouting for index mismatch
    % needs to developed here.
    for i = 1:length(dataSets)
        dataSets(i).zLoc_driftCorr_indexCorr = dataSets(i).zLoc_driftCorr * nSample/nOil;
    end
    

    %% Show the difference between the registered fiducial tracks
    figure_h_g = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
    set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
        'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
        'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
        'DefaultAxesLineWidth',1.2);
    
    subplot(4,1,1)
    plot(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
    xlabel('frame number');ylabel('X position shift (nm)');
    avg = nanmean(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
    stdev = nanstd(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
    ylim([avg-5*stdev avg+7*stdev])
    hold on
    plot(dataSets(2).fidTrack_interpolated_TRE(:,2),'red','LineWidth' ,2)
    hold off
    legend({['offset = ' num2str(avg),...
        ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_x']})
    title({['Fused Fiducial Tracks (raw data)']})
    
    subplot(4,1,2)
    plot(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
    xlabel('frame number');ylabel('Y position shift (nm)');
    avg = nanmean(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
    stdev = nanstd(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
    ylim([avg-5*stdev avg+7*stdev])
    hold on
    plot(dataSets(2).fidTrack_interpolated_TRE(:,3),'red','LineWidth' ,2)
    hold off
    legend({['offset = ' num2str(avg),...
        ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_y']})
    
    subplot(4,1,3)
    plot(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
    xlabel('frame number');ylabel('Z position shift (nm)');
    avg = nanmean(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
    stdev = nanstd(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
    ylim([avg-5*stdev avg+7*stdev])
    hold on
    plot(dataSets(2).fidTrack_interpolated_TRE(:,4),'red','LineWidth' ,2)
    hold off
    legend({['offset = ' num2str(avg),...
        ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_z']})
    
    subplot(4,1,4)
    euclid_Dist = sqrt(((dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed).^2)+...
        ((dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed).^2+...
        ((dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed).^2)));
    plot(euclid_Dist);
    xlabel('frame number');ylabel('3D shift (nm)');
    avg = nanmean(euclid_Dist);
    stdev = nanstd(euclid_Dist);
    ylim([0 avg+7*stdev])
    hold on
    plot(dataSets(2).fidTrack_interpolated_TRE(:,1),'red','LineWidth' ,2)
    hold off
    legend({['offset = ' num2str(avg),...
        ' +/- ' num2str(stdev) ' nm']; ['interpolated TRE_3_D']})
    clear avg
    
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
    
 %%   
%     %% Show the difference between the denoised and fused fiducial tracks
%     figure_h_h = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
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
%     
%     %% Show the fiducial tracks in the same coordinate system
%     figure_h_i = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
%     set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
%         'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
%         'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
%         'DefaultAxesLineWidth',1.2);
%     
%     for i = 1:length(dataSets)
%         if dataSets(i).transformedDataset
%             scatter3(dataSets(i).fidTrackX_transformed,...
%                 dataSets(i).fidTrackY_transformed,...
%                 dataSets(i).fidTrackZ_transformed,...
%                 5,'filled', color{i});
%         else
%             scatter3(dataSets(i).fidTrackX,...
%                 dataSets(i).fidTrackY,...
%                 dataSets(i).fidTrackZ,...
%                 5,'filled', color{i});
%         end
%         
%         xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
%         axis vis3d equal;
%         hold on
%         
%     end
%     hold off
%     
%     h = uicontrol('Position',[20 20 200 40],'String','Continue',...
%         'Callback','uiresume(gcbf)');
%     uiwait(gcf);
%     
    %% prompt to save bead registration figures
    [saveFile, savePath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
    savePath = [savePath saveFile '/'];
    mkdir(savePath);
    
    
    
    
    saveas(figure_h_a,[savePath 'XYZFidTracks_raw.fig']);
    saveas(figure_h_a,[savePath 'XYZFidTracks_raw.png']);
    close(figure_h_a)
    saveas(figure_h_b,[savePath 'XYZFidTracks_piecewiseDenoised.fig']);
    saveas(figure_h_b,[savePath 'XYZFidTracks_piecewiseDenoised.png']);
    close(figure_h_b)
    saveas(figure_h_c,[savePath 'ChromaticShiftDist.fig']);
    saveas(figure_h_c,[savePath 'ChromaticShiftDist.png']);
    close(figure_h_c)
    saveas(figure_h_d,[savePath 'XYZFidTracks_corrected.fig']);
    saveas(figure_h_d,[savePath 'XYZFidTracks_corrected.png']);
    close(figure_h_d)
    saveas(figure_h_e,[savePath 'XYZFidTracks_correctedDenoised.fig']);
    saveas(figure_h_e,[savePath 'XYZFidTracks_correctedDenoised.png']);
    close(figure_h_e)
    saveas(figure_h_f,[savePath 'AppliedDriftCorrection.fig']);
    saveas(figure_h_f,[savePath 'AppliedDriftCorrection.png']);
    close(figure_h_f)
    saveas(figure_h_g,[savePath 'XYZFidMisregistration.fig']);
    saveas(figure_h_g,[savePath 'XYZFidMisregistration.png']);
    close(figure_h_g)
%     saveas(figure_h_h,[savePath 'XYZFidMisregistration_denoised.fig']);
%     saveas(figure_h_h,[savePath 'XYZFidMisregistration_denoised.png']);
%     close(figure_h_h)
%     saveas(figure_h_i,[savePath '3DFidCorrelation.fig']);
%     saveas(figure_h_i,[savePath '3DFidCorrelation.png']);
%     close(figure_h_i)  

    
    %% Clean up
    tform.FRE = FRE;
    tform.TRE = TRE;
    tform.FRE_full = FRE_full;
    tform.TRE_full = TRE_full;
    tform.matched_cp_reflected = matched_cp_reflected;
    tform.matched_cp_transmitted = matched_cp_transmitted	;
    tform.matched_cp_transmitted_trans = matched_cp_transmitted_trans;
    
    % Need to clean up further
    fidCorrData.avgDevX = avgDevX;
    
    
    
    registrationComplete = true
    
    clear FRE TRE FRE_full TRE_full matched_cp_reflected matched_cp_transmitted matched_cp_transmitted_trans
    clear avgDevX avgDevY avgDevZ euclid_Dist goodFits h stdev
    clear globalScale kthNeighbor nControlPoints nCores
    clear figure_h_a figure_h_b figure_h_c figure_h_d
    clear devX devX_denoised devY devY_denoised devZ devZ_denoised
    clear fidTracksX fidTracksX_denoised fidTracksY fidTracksY_denoised fidTracksZ fidTracksZ_denoised
    clear x y z xLocPix yLocPix zLoc_IndexCorrected
    close all

    save([savePath 'Combined_MulticolorSMACMData.mat']);
    
end

%% Display the results the final fused SMACM data
useTimeColors = 0;
frameRange = [1 100000];
numPhotonRange = [0 100000];

dlg_title = 'Please Input Parameters';
prompt = {  'Pixel size (in nm)',...
    'Size of Points in reconstruction',...
    'White Light Shift X (in nm)',...
    'White Light Shift Y (in nm)',...
    };
def = {    '125.78', ...
    '30', ...
    num2str(dataSets(untransformedDataSet).wlShift(1)), ...
    num2str(dataSets(untransformedDataSet).wlShift(2)), ...
    };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

nmPerPixel = str2double(inputdialog{1});
scatterSize = str2double(inputdialog{2});
wlShiftX = str2double(inputdialog{3});
wlShiftY = str2double(inputdialog{4});

pass = 1;
anotherpass = true;

while anotherpass == true
    close all
    
    %% Plot the white light image if specified
    if whiteLightFile ~= 0
        whiteLightInfo = imfinfo([whiteLightPath whiteLightFile]);
        whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
        % average white light images together to get better SNR
        for a = 1:length(whiteLightInfo)
            whiteLight = whiteLight + double(imread([whiteLightPath whiteLightFile], ...
                'Info', whiteLightInfo));
        end
        % resize white light to the size of the ROI of the single molecule fits
        ROI_initial = [1, 1, 270, 270];
        whiteLight = whiteLight(ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1,ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1);
        % rescale white light image to vary from 0 to 1
        whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
        [xWL yWL] = meshgrid((ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1) * nmPerPixel + wlShiftX, ...
            (ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1) * nmPerPixel + wlShiftY);
    end
    
    
    %% Chose a desired parameter set for reconstruction
    dlg_title = 'Please Input Parameters';
    prompt = {  'Size of points in reconstruction',...
        'Temporal Color Coding',...
        'White light shift X (in nm)',...
        'White light shift Y (in nm)',...
        'First frame',...
        'Last frame',...
        'Number of photons lower bound',...
        'Number of photons upper bound',...
        };
    def = { ...
        num2str(scatterSize), ...
        num2str(useTimeColors), ...
        num2str(wlShiftX), ...
        num2str(wlShiftY), ...
        num2str(frameRange(1)), ...
        num2str(frameRange(2)), ...
        num2str(numPhotonRange(1)), ...
        num2str(numPhotonRange(2)), ...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    scatterSize = str2double(inputdialog{1});
    useTimeColors = str2double(inputdialog{2});
    wlShiftX = str2double(inputdialog{3});
    wlShiftY = str2double(inputdialog{4});
    frameRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
    numPhotonRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];
    
    
    %% ask user what region to plot in superresolution image
    
    if whiteLightFile ~= 0
        xRange = xWL(1,:);
        yRange = yWL(:,1);
        % pick region that contains background
        figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(xRange,yRange,whiteLight);axis image;colormap gray;
        hold on;
    end
    
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
        colormap gray; hold on;
    end
    
    for i = 1:length(dataSets)
        
        xLoc = dataSets(i).xLoc_driftCorr;
        yLoc = dataSets(i).yLoc_driftCorr;
        zLoc = dataSets(i).zLoc_driftCorr;
        zLoc_indexCorr = dataSets(i).zLoc_driftCorr;
        
        validPoints = inpolygon(xLoc,yLoc,xi, yi);
        
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
        if dataSets(i).transformedDataset == 1
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
            ['Mean FRE = ' num2str(mean(interpolated_FREs)) ' nm'];...
            ['Mean TRE = ' num2str(mean(interpolated_TREs)) ' nm'];...
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
            pass = pass + 1;
        case 'No'
            anotherpass = false;
        case 'Cancel'
            error('User cancelled the program');
    end
    
end

%% clean up

tform.interpolationObjects.F_FRE = F_FRE;
tform.interpolationObjects.F_FRE_X = F_FRE_X;
tform.interpolationObjects.F_FRE_Y = F_FRE_Y;
tform.interpolationObjects.F_FRE_Z = F_FRE_Z;
tform.interpolationObjects.F_TRE = F_TRE;
tform.interpolationObjects.F_TRE_X = F_TRE_X;
tform.interpolationObjects.F_TRE_Y = F_TRE_Y;
tform.interpolationObjects.F_TRE_Z = F_TRE_Z;

clear F_FRE F_FRE_X F_FRE_Y F_FRE_Z F_TRE F_TRE_X F_TRE_Y F_TRE_Z
clear a anotherpass bead croppedDataSet def dlg_title f h i
clear interpolated_FRE interpolated_FREs interpolated_TRE interpolated_TREs
clear num_lines pass prompt questiondialog
clear xLoc xslice yLoc yslice zLoc zLoc_indexCorr zslice
clear x xRange xWL xi y yRange yWL yi z

%% prompt to save data and figures
[saveFile, savePath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
savePath = [savePath saveFile '/'];
mkdir(savePath);

if ~isequal(saveFile,0)
    save([savePath saveFile '_multicolorSMACM.mat']);
end

saveas(gcf,[savePath saveFile(1:length(saveFile)-4) '_multicolorSMACM_3D.fig']);
close
saveas(gcf,[savePath saveFile(1:length(saveFile)-4) '_multicolorSMACM_2D.fig']);
close all

end

