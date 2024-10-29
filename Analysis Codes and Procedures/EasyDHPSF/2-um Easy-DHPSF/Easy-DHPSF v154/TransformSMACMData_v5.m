%%%% 2018-11-07 - for use with 2 color DH data and DH fiducial in the same channel
% corrects drift, applies index mismatch correction, and calculates
% sigmaXYZ
%%%% 2019-01-31 register data, register fiducial, then apply relative
%%%% fiducial correction to registered data

clear all
close all

scrsz = get(0,'ScreenSize');
nSample = 1.33;         % index of refraction of sample
nOil = 1.518;           % index of immersion oil

numChannels = 2;

numSyncFrames = 25;
useDenoising = 1;
moleFiles = {};
fidFiles = {};

% use local cals questdlg
dlg_title = 'Use local cals?';
prompt = { 'Do you want to try to use spatially dependent calibrations?' };
def =       { 'No'  };
questiondialog = questdlg(prompt,dlg_title, def);
% Handle response
switch questiondialog
    case 'Yes'
        spatialCorr=true;
    case 'No'
        spatialCorr=false;
    case 'Cancel'
        error('User cancelled the program');
end

for ii=1:numChannels
    
    if ii == 1
        [easyFile, easyPath] = uigetfile({'*.mat';'*.*'}, 'Open green easy DHPSF file');
    elseif ii == 2
        [easyFile, easyPath] = uigetfile({'*.mat';'*.*'}, 'Open red easy DHPSF file');
    end
    
    load([easyPath easyFile]);
    
    fitFilePrefix = s.fitFilePrefix;
    fidFilePrefix = s.fidFilePrefix;
    channel = s.channel;
    calFile = [s.calFilePrefix 'calibration.mat'];
    useFidCorrections = 1;
    
    %% Get the raw fiduciary data from the easy DHPSf file:
    for fileNum=1:length(fidFilePrefix)
        %         fidFiles = [fidFiles; {[fidPath fidFile]}];
        % load data
        load([fidFilePrefix{fileNum} 'raw fits.mat'],'PSFfits','numFrames','numMoles');
        
        fidFiles = [fidFiles; {[fidFilePrefix{fileNum}]}]
        
        if fileNum == 1
            tempPSFfits = PSFfits(:,1:23);
            numFramesInFiles = numFrames;
        else
            numFramesInFiles = [numFramesInFiles numFrames];
            PSFfits(:,1) = PSFfits(:,1) + sum(numFramesInFiles(1:fileNum-1));
            tempPSFfits = [tempPSFfits; PSFfits(:,1:23)];
        end
        
    end
    PSFfits = tempPSFfits;
    numFrames = sum(numFramesInFiles);
    clear tempPSFfits;
    
    %% compute movement of fiduciaries
    
    fidTrackX = NaN(size(PSFfits,1),numMoles);
    fidTrackY = NaN(size(PSFfits,1),numMoles);
    fidTrackZ = NaN(size(PSFfits,1),numMoles);
    
    devX = zeros(numFrames,numMoles);
    devY = zeros(numFrames,numMoles);
    devZ = zeros(numFrames,numMoles);
    goodFitFlag = zeros(numFrames,numMoles);
    numPhotons = zeros(numFrames,numMoles);
    avgDevX = zeros(numFrames,1);
    avgDevY = zeros(numFrames,1);
    avgDevZ = zeros(numFrames,1);
    numValidFits = zeros(numFrames,1);
    
    syncFrames = zeros(1,numSyncFrames);
    lastGoodFrame = numFrames;
    for a = 1:numSyncFrames
        while sum(PSFfits(PSFfits(PSFfits(:,1)==lastGoodFrame,13)>0,2)) ~= numMoles/2*(1+numMoles)
            lastGoodFrame = lastGoodFrame - 1;
        end
        syncFrames(a)=lastGoodFrame;
        lastGoodFrame = lastGoodFrame - 1;
    end
    
    for molecule = 1:numMoles
        %% extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        
        
        %% Correct the moleculeFitParam for the index mismatch
        
        % compute deviation with respect to bead location averaged over last
        % numSyncFrames frames of the movie
        devX(:,molecule) = moleculeFitParam(:,21) -nanmean(moleculeFitParam(goodFit,21));
        devY(:,molecule) = moleculeFitParam(:,22) -nanmean(moleculeFitParam(goodFit,22));
        devZ(:,molecule) = moleculeFitParam(:,23) -nanmean(moleculeFitParam(goodFit,23));
        % if particle was fit successfully, add its movement to the average
        avgDevX = avgDevX + goodFit.*devX(:,molecule);
        avgDevY = avgDevY + goodFit.*devY(:,molecule);
        avgDevZ = avgDevZ + goodFit.*devZ(:,molecule);
        numValidFits = numValidFits + goodFit;
    end
    tempAvgDevX = avgDevX./numValidFits;
    tempAvgDevY = avgDevY./numValidFits;
    tempAvgDevZ = avgDevZ./numValidFits;
    
    
    for molecule = 1:numMoles
        % extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        % only use good fits as defined by fitting function
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        
        % raw positions of the fiducial tracks
        fidTrackX(goodFit,molecule) = moleculeFitParam(goodFit,21);
        fidTrackY(goodFit,molecule) = moleculeFitParam(goodFit,22);
        fidTrackZ(goodFit,molecule) = moleculeFitParam(goodFit,23);
        %     numPhotons(:,molecule) = moleculeFitParam(:,17)
        
    end
    
    for cc = 1:fileNum
        %         [moleFile, molePath] = uigetfile({'*.mat';'*.*'},...
        %             ['Open data file #' num2str(cc) ' with PSF localizations']);
        %         load([fitFilePrefix{cc} 'molecule fits.mat'],'totalPSFfits','numFrames');
        load([fitFilePrefix{cc} 'molecule fits.mat']);
        
        moleFiles = [moleFiles; {[fitFilePrefix{cc}]}];
        % load data
        %         load([molePath moleFile]);
        
        if cc == 1
            tempPSFfits = totalPSFfits(:,1:27);
        else
            totalPSFfits(:,1) = totalPSFfits(:,1) + sum(numFramesInFiles(1:cc-1));
            tempPSFfits = [tempPSFfits; totalPSFfits(:,1:27)];
        end
    end
    totalPSFfits = tempPSFfits;
    numFrames = sum(numFramesInFiles);
    
    
    if spatialCorr
        totalPSFfits = f_makeLocalCals(totalPSFfits,calFile,'SMACM',0,1);
    end
    
    
    clear tempPSFfits;
    
    % De-noise the fiduciary tracks
    avgDevX = tempAvgDevX;
    avgDevY = tempAvgDevY;
    avgDevZ = tempAvgDevZ;
    
    [avgDevX_denoised,avgDevY_denoised,avgDevZ_denoised] = f_waveletFidTracks(avgDevX,avgDevY,avgDevZ,0);
    
    %%%%
    if ii ==1
        avgDevX_denoised = avgDevX_denoised - avgDevX_denoised(1,1);
        avgDevY_denoised = avgDevY_denoised - avgDevY_denoised(1,1);
        avgDevZ_denoised = avgDevZ_denoised - avgDevZ_denoised(1,1);
%     elseif ii == 2
%         avgDevX_denoised = avgDevX_denoised - avgDevX_denoised(end,1);
%         avgDevY_denoised = avgDevY_denoised - avgDevY_denoised(end,1);
%         avgDevZ_denoised = avgDevZ_denoised - avgDevZ_denoised(end,1);
    end
    %%%%
        
    
    % apply fiduciary corrections
    if useDenoising == 1
        xFidCorrected = totalPSFfits(:,25) - avgDevX_denoised(totalPSFfits(:,1));
        yFidCorrected = totalPSFfits(:,26) - avgDevY_denoised(totalPSFfits(:,1));
        zFidCorrected = totalPSFfits(:,27) - avgDevZ_denoised(totalPSFfits(:,1));
    else
        xFidCorrected = totalPSFfits(:,25) - avgDevX(totalPSFfits(:,1));
        yFidCorrected = totalPSFfits(:,26) - avgDevY(totalPSFfits(:,1));
        zFidCorrected = totalPSFfits(:,27) - avgDevZ(totalPSFfits(:,1));
    end
    totalPSFfits = [totalPSFfits xFidCorrected yFidCorrected zFidCorrected];
    clear tempAvgDevX tempAvgDevY tempAvgDevZ;
    
    if ii==1 % green
        totalPSFfits_green = totalPSFfits;
        avgDevX_denoised_green = avgDevX_denoised;
        avgDevY_denoised_green = avgDevY_denoised;
        avgDevZ_denoised_green = avgDevZ_denoised;
    elseif ii == 2 % red
        totalPSFfits_red = totalPSFfits;
        avgDevX_denoised_red = avgDevX_denoised;
        avgDevY_denoised_red = avgDevY_denoised;
        avgDevZ_denoised_red = avgDevZ_denoised;
    end
    
end
    
%% load in the transform 
[tformFile, tformPath] = uigetfile({'*.mat';'*.*'},'Open 3D_Transform.mat');
if isequal(tformFile,0)
    error('User cancelled the program');
end

% tformFile = '3D_Transform_lwquadratic.mat';
% tformPath = 'C:\UserFiles\UsersNOTBackedUp\Temp\2014-04-09_L1-PAM-B-eYFP\FittingCode-Sept3rd-2014\';

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


% Columns: 
frameCol = 1;
matchConCol = 6; 
amp1Col = 7; 
amp2Col = 8; 
xLoc1Col = 9; 
yLoc1Col = 10; 
xLoc2Col = 11; 
yLoc2Col = 12; 
sigma1Col = 13; 
sigma2Col = 14; 
bkgndMeanCol = 15; 
fitErrorCol = 16; 
goodFitFlagCol = 17; 
xCenterCol = 18; 
yCenterCol = 19; 
AngleCol = 20; 
PhotonCol = 21;
interlobeCol = 22; 
ampRatioCol = 23; 
sigmaRatioCol = 24; 
xCol = 25;
yCol = 26;
zCol = 27;
xFidCol = 28;
yFidCol = 29; 
zFidCol = 30; 

fileNum = 1;
LocFiles = {};
dataSets = [];

%% load in the GREEN .mat file
% Add it to the total data structure 
% clear everything 
% [LocFile, LocPath] = uigetfile({'*.mat';'*.*'},'Load in the GREEN data');
% if isequal(LocFile,0)
%     error('User cancelled the program');
% end
% load([LocPath LocFile]);
% 
% LocFiles = [LocFiles; {[LocPath LocFile]}];

clear totalPSFfits
totalPSFfits = totalPSFfits_green;
clear avgDevX_denoised avgDevY_denoised avgDevZ_denoised
avgDevX_denoised = avgDevX_denoised_green;
avgDevY_denoised = avgDevY_denoised_green;
avgDevZ_denoised = avgDevZ_denoised_green;

% crop dataset
dlg_title = 'Crop dataset';
prompt = {'Do you want to crop this dataset?'};
def =       { 'Yes' };
questiondialog = questdlg(prompt,dlg_title, def);

switch questiondialog
    case 'Yes'
        figure;
        scatter(totalPSFfits(:,xCol),totalPSFfits(:,yCol),'.');
        axis image;hold on;
        title('click an ROI around the part you want to crop')
        
        [xi,yi] = ginput;
        xi = [xi; xi(1)];
        yi = [yi; yi(1)];
        plot(xi, yi, 'k','LineWidth',2);
        
        hold off
        
        validPoints = inpolygon(totalPSFfits(:,xCol),totalPSFfits(:,yCol),xi, yi);
        
        temp = totalPSFfits(validPoints,:); clear totalPSFfits;
        totalPSFfits = temp; clear temp;
    case 'No'
        % dont crop
end

goodidxnotnan = ~isnan(totalPSFfits(:,zCol));
tempPSFfits = totalPSFfits(goodidxnotnan,:);
clear totalPSFfits
totalPSFfits = tempPSFfits;
clear tempPSFfits

dataSet.frame = totalPSFfits(:,frameCol); 
dataSet.matchConf = totalPSFfits(:,matchConCol); 
dataSet.amp1 = totalPSFfits(:,amp1Col); 
dataSet.amp2 = totalPSFfits(:,amp2Col); 
dataSet.xLoc1 = totalPSFfits(:,xLoc1Col); 
dataSet.yLoc1 = totalPSFfits(:,yLoc1Col); 
dataSet.xLoc2 = totalPSFfits(:,xLoc2Col); 
dataSet.yLoc2 = totalPSFfits(:,yLoc2Col); 
dataSet.sigma1 = totalPSFfits(:,sigma1Col); 
dataSet.sigma2 = totalPSFfits(:,sigma2Col); 
dataSet.bkgndMean = totalPSFfits(:,bkgndMeanCol); 
dataSet.fitError = totalPSFfits(:,fitErrorCol); 
dataSet.goodFitFlag = totalPSFfits(:,goodFitFlagCol); 
dataSet.xCenter = totalPSFfits(:,xCenterCol); 
dataSet.yCenter = totalPSFfits(:,yCenterCol); 
dataSet.angle = totalPSFfits(:,AngleCol); 
dataSet.photons = totalPSFfits(:,PhotonCol); 
dataSet.interlobe = totalPSFfits(:,interlobeCol); 
dataSet.ampRatio = totalPSFfits(:,ampRatioCol); 
dataSet.sigmaRatio = totalPSFfits(:,sigmaRatioCol); 
dataSet.x = totalPSFfits(:,xCol); 
dataSet.y = totalPSFfits(:,yCol); 
dataSet.z = totalPSFfits(:,zCol); 
dataSet.LocFile = easyFile; 
dataSet.LocPath = easyPath; 

if size(totalPSFfits,2)>27
    dataSet.xFid = totalPSFfits(:,xFidCol); 
    dataSet.yFid = totalPSFfits(:,yFidCol); 
    dataSet.zFid = totalPSFfits(:,zFidCol); 
    dataSet.fidTrackX = avgDevX_denoised; 
    dataSet.fidTrackY = avgDevY_denoised; 
    dataSet.fidTrackZ = avgDevZ_denoised; 
else
    dataSet.xFid = NaN; 
    dataSet.yFid = NaN; 
    dataSet.zFid = NaN; 
    dataSet.fidTrackX = NaN; 
    dataSet.fidTrackY = NaN; 
    dataSet.fidTrackZ = NaN;     
end

% transform green data: 
% plot(totalPSFfits(:,xFidCol),totalPSFfits(:,yFidCol),'b.'); 
% 
% transformedData = transformData([totalPSFfits(:,xCol),totalPSFfits(:,yCol),totalPSFfits(:,zCol)],tform); 
% % transformedData = transformData([totalPSFfits(:,xFidCol),totalPSFfits(:,yFidCol),totalPSFfits(:,zFidCol)],tform); 
dataSet.xTrans = NaN; 
dataSet.yTrans = NaN; 
dataSet.zTrans = NaN; 

% Append the structure to the previous structures in the array
dataSets = [dataSets, dataSet];

clear totalPSFfits avgDevX_denoised avgDevY_denoised avgDevZ_denoised 

%% load in the RED .mat file
% Add it to the total data structure 
% clear everything 

% [LocFile, LocPath] = uigetfile({'*.mat';'*.*'},'Load in the RED data');
% if isequal(LocFile,0)
%     error('User cancelled the program');
% end
% load([LocPath LocFile]);
% 
% LocFiles = [LocFiles; {[LocPath LocFile]}];

clear totalPSFfits
totalPSFfits = totalPSFfits_red;
clear avgDevX_denoised avgDevY_denoised avgDevZ_denoised
avgDevX_denoised = avgDevX_denoised_red;
avgDevY_denoised = avgDevY_denoised_red;
avgDevZ_denoised = avgDevZ_denoised_red;

% crop dataset
dlg_title = 'Crop dataset';
prompt = {'Do you want to crop this dataset?'};
def =       { 'Yes' };
questiondialog = questdlg(prompt,dlg_title, def);

switch questiondialog
    case 'Yes'
        figure;
        scatter(totalPSFfits(:,xCol),totalPSFfits(:,yCol),'.');
        axis image;hold on;
        title('click an ROI around the part you want to crop')
        
        [xi,yi] = ginput;
        xi = [xi; xi(1)];
        yi = [yi; yi(1)];
        plot(xi, yi, 'k','LineWidth',2);
        
        hold off
        
        validPoints = inpolygon(totalPSFfits(:,xCol),totalPSFfits(:,yCol),xi, yi);
        
        temp = totalPSFfits(validPoints,:); clear totalPSFfits;
        totalPSFfits = temp; clear temp;
    case 'No'
        % dont crop
end

goodidxnotnan = ~isnan(totalPSFfits(:,zCol));
tempPSFfits = totalPSFfits(goodidxnotnan,:);
clear totalPSFfits
totalPSFfits = tempPSFfits;
clear tempPSFfits

dataSet.frame = totalPSFfits(:,frameCol);           % 1
dataSet.matchConf = totalPSFfits(:,matchConCol);    % 2 
dataSet.amp1 = totalPSFfits(:,amp1Col);             % 3
dataSet.amp2 = totalPSFfits(:,amp2Col);             % 4
dataSet.xLoc1 = totalPSFfits(:,xLoc1Col);           % 5
dataSet.yLoc1 = totalPSFfits(:,yLoc1Col);           % 6
dataSet.xLoc2 = totalPSFfits(:,xLoc2Col);           % 7
dataSet.yLoc2 = totalPSFfits(:,yLoc2Col);           % 8
dataSet.sigma1 = totalPSFfits(:,sigma1Col);         % 9
dataSet.sigma2 = totalPSFfits(:,sigma2Col);         % 10 
dataSet.bkgndMean = totalPSFfits(:,bkgndMeanCol);   % 11
dataSet.fitError = totalPSFfits(:,fitErrorCol);     % 12
dataSet.goodFitFlag = totalPSFfits(:,goodFitFlagCol); % 13
dataSet.xCenter = totalPSFfits(:,xCenterCol);         % 14
dataSet.yCenter = totalPSFfits(:,yCenterCol);         % 15
dataSet.angle = totalPSFfits(:,AngleCol);             % 16
dataSet.photons = totalPSFfits(:,PhotonCol);          % 17
dataSet.interlobe = totalPSFfits(:,interlobeCol);     % 18
dataSet.ampRatio = totalPSFfits(:,ampRatioCol);       % 19
dataSet.sigmaRatio = totalPSFfits(:,sigmaRatioCol);   % 20
dataSet.x = totalPSFfits(:,xCol);                     % 21
dataSet.y = totalPSFfits(:,yCol);                     % 22
dataSet.z = totalPSFfits(:,zCol);                     % 23
dataSet.LocFile = easyFile; 
dataSet.LocPath = easyPath; 

if size(totalPSFfits,2)>27
    dataSet.xFid = totalPSFfits(:,xFidCol); 
    dataSet.yFid = totalPSFfits(:,yFidCol); 
    dataSet.zFid = totalPSFfits(:,zFidCol); 
    dataSet.fidTrackX = avgDevX_denoised; 
    dataSet.fidTrackY = avgDevY_denoised; 
    dataSet.fidTrackZ = avgDevZ_denoised; 
else
    dataSet.xFid = NaN; 
    dataSet.yFid = NaN; 
    dataSet.zFid = NaN; 
    dataSet.fidTrackX = NaN; 
    dataSet.fidTrackY = NaN; 
    dataSet.fidTrackZ = NaN;     
end


% transform red data into green channel: 
transformedData = transformData([totalPSFfits(:,xCol),totalPSFfits(:,yCol),totalPSFfits(:,zCol)],tform); 
% transformedData = transformData([totalPSFfits(:,xFidCol),totalPSFfits(:,yFidCol),totalPSFfits(:,zFidCol)],tform); 
dataSet.xTrans = transformedData(:,1); 
dataSet.yTrans = transformedData(:,2); 
dataSet.zTrans = transformedData(:,3); 

% transform red fiducial into green channel:
transformedFid = transformData([avgDevX_denoised, avgDevY_denoised, avgDevZ_denoised],tform);

% zero red fid track to last red frame
dataSet.fidTrackX = nan; dataSet.fidTrackX = transformedFid(:,1);
dataSet.fidTrackX = dataSet.fidTrackX - dataSet.fidTrackX(end,1);
dataSet.fidTrackY = nan; dataSet.fidTrackY = transformedFid(:,2);
dataSet.fidTrackY = dataSet.fidTrackY - dataSet.fidTrackY(end,1);
dataSet.fidTrackZ = nan; dataSet.fidTrackZ = transformedFid(:,3);
dataSet.fidTrackZ = dataSet.fidTrackZ - dataSet.fidTrackZ(end,1);


% Append the structure to the previous structures in the array
dataSets = [dataSets, dataSet];

%% clear all except the dataSets

clearvars -except dataSets LocFiles

% RED TRANSFORMED GREEN NOT TRANSFORMED
% figure
% plot3(dataSets(1,2).xTrans,dataSets(1,2).yTrans,dataSets(1,2).zTrans,'r.')
% hold on
% plot3(dataSets(1,1).x,dataSets(1,1).y,dataSets(1,1).z,'g.')
% axis image 

% plot(dataSets(1,2).xTrans,dataSets(1,2).yTrans,'r.')
% hold on
% plot(dataSets(1,1).x,dataSets(1,1).y,'g.')
% axis image 


% correct fiducial in green channel
dataSets(1).xFid = dataSets(1).x - dataSets(1).fidTrackX(dataSets(1).frame); 
dataSets(1).yFid = dataSets(1).y - dataSets(1).fidTrackY(dataSets(1).frame); 
dataSets(1).zFid = dataSets(1).z - dataSets(1).fidTrackZ(dataSets(1).frame); 

% correct fiducial in red channel
dataSets(2).xFid = dataSets(2).xTrans - dataSets(2).fidTrackX(dataSets(2).frame); 
dataSets(2).yFid = dataSets(2).yTrans - dataSets(2).fidTrackY(dataSets(2).frame); 
dataSets(2).zFid = dataSets(2).zTrans - dataSets(2).fidTrackZ(dataSets(2).frame); 

figure; 
plot(dataSets(1,2).xFid,dataSets(1,2).yFid,'r.')
hold on
plot(dataSets(1,1).xFid,dataSets(1,1).yFid,'g.')
axis image 

% plot3(dataSets(1,2).xFid,dataSets(1,2).yFid,dataSets(1,2).zFid,'r.')
% hold on
% plot3(dataSets(1,1).xFid,dataSets(1,1).yFid,dataSets(1,1).zFid,'g.')
% axis image 

nSample = 1.33;         % index of refraction of sample
nOil = 1.518;           % index of immersion oil

%% apply index mismatch correction
for i = 1:length(dataSets)
        dataSets(i).zLoc_driftCorr_indexCorr = dataSets(i).zFid * nSample/nOil;
end 

%% calculate sigmaXYZ
% calculate localization precision
amplitude =  [361035.867260138,22.2956414971275;...   %   [A1x  A2x]
    348907.934759022,28.3183226442783;...   %   [A1y  A2y]
    840446.405407229,23.3314294806927];      %   [A1z  A2z] 8A back values

for i=1:length(dataSets)
    
% Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings

dataSets(i).sigmaX = sqrt(amplitude(1,1) .* (1./dataSets(i).photons) + ...
    amplitude(1,1)*4*amplitude(1,2) .* dataSets(i).bkgndMean./(dataSets(i).photons).^2 + ...
    amplitude(1,1) .* (1./dataSets(i).photons) .* sqrt((2*amplitude(1,2)*(dataSets(i).bkgndMean./dataSets(i).photons))./(1+(4*amplitude(1,2)*(dataSets(i).bkgndMean./dataSets(i).photons)))));

dataSets(i).sigmaY = sqrt(amplitude(2,1) .* (1./dataSets(i).photons) + ...
    amplitude(2,1)*4*amplitude(2,2) .* dataSets(i).bkgndMean./(dataSets(i).photons).^2 + ...
    amplitude(2,1) .* (1./dataSets(i).photons) .* sqrt((2*amplitude(2,2)*(dataSets(i).bkgndMean./dataSets(i).photons))./(1+(4*amplitude(2,2)*(dataSets(i).bkgndMean./dataSets(i).photons)))));

dataSets(i).sigmaZ = sqrt(amplitude(3,1) .* (1./dataSets(i).photons) + ...
    amplitude(3,1)*4*amplitude(3,2) .* dataSets(i).bkgndMean./(dataSets(i).photons).^2 + ...
    amplitude(3,1) .* (1./dataSets(i).photons) .* sqrt((2*amplitude(3,2)*(dataSets(i).bkgndMean./dataSets(i).photons))./(1+(4*amplitude(3,2)*(dataSets(i).bkgndMean./dataSets(i).photons)))));

end

%% save the data 
[saveFilePrefix, savePath] = uiputfile({'*.*'},'Enter a prefix for the combined data sets');
if isequal(saveFilePrefix,0)
    error('User cancelled the program');
end
saveFile = [savePath saveFilePrefix];


save([saveFile '.mat'])