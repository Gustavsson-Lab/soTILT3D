% Copyright (c)2018, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names
% of its contributors may be used to endorse or promote products derived
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [outputFilePrefix,logFile,logPath] = ...
    f_fitSMs(dataFile,dataPath,calFile,calBeadIdx,templateFrames,peakThreshold,...
    darkFile,logFile,logPath,boxRadius,channel, sigmaBounds,gaussianFilterSigma,minDistBetweenSMs,...
    lobeDistBounds,conversionGain,nmPerPixel,EMGain,templateLocs,threshFile,ROI,nhaData,angledNHA)
% f_fitSMs is a module in easy_dhpsf that finds the positions of likely
% DH-PSF profiles by matching to a series of templates generated in
% f_calDHPSF and prepared in f_calSMidentification. These are then more
% precisely localized using a double gaussian fit and corrected for drift
% using the results from f_trackFiducials.

if exist('threshFile')
    load(threshFile,'usePolyROI','medianFilter','windowSize','medianBlurSigma','dispRaw','interpVal')
    if usePolyROI
        load(threshFile,'FOVmask','FOVmask1');
    end
end
if nhaData
    peakThreshold=peakThreshold*10000; %cancels out the call to divide in easy_dhpsf
end

printOutputFrames = 0;

if printOutputFrames == 1 % this will save all process images (correlation image, raw data, and reconstruction)
    mkdir('output images');
end

% initialize parameters
scrsz = get(0,'ScreenSize');
conversionFactor = conversionGain/EMGain;

ampRatioLimit = 0.5;
if nhaData
    ampRatioLimit = 0.7;
end
sigmaRatioLimit = 0.4;

% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');


%% preprocessing for files

outputFilePrefix = cell(1,length(dataFile));

% can comment the code below in if it is desirable to select specific files
% see f_calSMidentification for how to do this

% if length(dataFile) > 1
%     dlg_title = 'Select Files';
%     prompt = {  'Choose files for thresholding' };
%     def = {num2str(1:length(dataFile))};
%     num_lines = 1;
%     inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
%     if isempty(inputdialog)
%         error('User cancelled the program')
%     end
%     selectedFiles = str2num(inputdialog{1});
% else
%     selectedFiles = 1:length(dataFile);
% end

selectedFiles = 1:length(dataFile);

% allows user to select subset of frames, at the beginning, for all files
%dlg_title = 'Select Frames (default: use all frames)';
dlg_title = 'Estimate laser profile?';
num_lines = 1;
def = {};
prompt = {};

%%%%%%%%%%%%%%%%
% def = {'Yes'};
% prompt = {'Do you want to estimate the laser profile?'};
%%%%%%%%%%%%%%% this comment and the one below can be uncommented to allow
%%%%%%%%%%%%%%% selection of specific frames: however, this is *slow*
% populates 'fileInfo' and 'numFrames' for all files, and generates the
% fields for the frame selection dlg
fileInfoAll=cell(length(selectedFiles),1);
numFramesAll = zeros(length(selectedFiles),1);
def = cell(length(selectedFiles)+2,1);
prompt = cell(length(selectedFiles)+2,1);

for i = 1:length(selectedFiles)
    fileInfoAll{i} = imfinfo([dataPath dataFile{selectedFiles(i)}]);
    numFramesAll(i) = length(fileInfoAll{i});
    def{i} = ['[1:' num2str(numFramesAll(i)) ']'];
    prompt{i} = ['Choose frames for ' dataFile{selectedFiles(i)}];
end
    def{end-1} = 'No';
    prompt{end-1} = 'Do you want to estimate the laser profile?';
    def{end} = 'No';
    prompt{end} = 'Do you want to try to find true lobe sigmas?';

%%% edit so that dialog box fits in screen
if length(selectedFiles) <= 10
    % small number of files, display all at once
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    allinputdialog = inputdialog;
else
    allinputdialog = {};
    for i =1:ceil(length(prompt)/10)
        startFile = (i-1)*10 + 1; 
        if i == ceil(length(prompt)/10)
            endFile = length(prompt);
        else
            endFile = startFile + 9;
        end
        inputdialog = inputdlg(prompt(startFile:endFile), dlg_title, num_lines, def(startFile:endFile));
        allinputdialog = [allinputdialog; inputdialog];
    end
end    
    
% inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

%%%%%%%%%%%%%%%** this must be uncommented if selecting frames **
for i = 1:length(selectedFiles)
    framesAll{i} = str2num(allinputdialog{i});
end
findLaserInt = ~strcmp(allinputdialog{end-1},'No');
findTrueSigma = ~strcmp(allinputdialog{end},'No');

if findTrueSigma % this uses constrained sigmas for the fit of location, etc., and then separately finds the sigmas
    sigSearchBounds = [0.9 4];
    sigOptions = optimset('FunValCheck','on','Diagnostics','off','Jacobian','off', 'Display', 'off');
end
% This information should be passed from the earlier execution of
% this code in f_calSMidentification, but we get another chance to
% select it if it wasn't specified
if isequal(logFile,0)
    [logFile, logPath] = uigetfile({'*.dat';'*.*'},...
        'Open sequence log file(s) corresponding to image stack(s) (optional: hit cancel to skip)',...
        'MultiSelect', 'on');
    if isequal(logPath,0)
        logFile = 'not specified';
    end
    if ischar(logFile)
        logFile = cellstr(logFile);
    end
end
% sets absolute frame number for relating data to sequence log
absFrameNum = 1;

dlg_title = 'Filtering control points?';
prompt = {  'If control points: how many stationary frames for each position? (Cancel if not CP)',...
    };
def = {    '20', ...
    };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(inputdialog)
    maxNumMeasurement = str2double(inputdialog{1});
    filterCPs = true;
    peakThreshold = peakThreshold*0.25; %%% edit to make threshold smaller - useful for beads during CP scan
else
    filterCPs = false;
end


%% begin fitting loop over files
for stack = selectedFiles % = 1:length(dataFile)
    
    %%% this is an efficient way to set these variables if frames were
    %%% selected above
        fileIdx = find(selectedFiles == stack);
        fileInfo = fileInfoAll{fileIdx};
        frames = framesAll{fileIdx};
        numFrames = numFramesAll(fileIdx);
    
    %%% these variables must be set this way if frames were not selected
%     fileInfo = imfinfo([dataPath dataFile{stack}]);
%     numFrames = length(fileInfo);
%     frames = 1:numFrames;
    
    imgHeight = fileInfo(1).Height;
    imgWidth = fileInfo(1).Width;
    %% create output log filenames
    
    % saves in labeled directory if a channel is selected
    if channel == '0'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'molecule fits  ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
    else
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep channel(1) ' molecule fits  ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
    end
    mkdir(outputFilePrefix{stack});
    
    if stack == 1 %%% begin init
        %% miscellaneous bookkeeping on the darkfile / template / data
        % Compute darkAvg counts
        
        if ~isequal(darkFile,0)
            % Computes average of darkAvg frames for background subtraction
            darkFileInfo = imfinfo(darkFile);
            numDarkFrames = length(darkFileInfo);
            darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
%             im = tiffread2(darkFile,1,numDarkFrames);
            for frame = 1:numDarkFrames
%                 darkAvg = darkAvg + double(im(frame).data);
                  darkAvg = darkAvg + double(imread(darkFile,frame,'Info',darkFileInfo));
            end
            darkAvg = darkAvg/numDarkFrames;
            if ~isequal(size(darkAvg),[imgHeight imgWidth])
                warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
                darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
            end
        else
            darkAvg = 0;
        end
        clear darkFileInfo
        
        % Compute average image
        
        %         avgImg = zeros(imgHeight,imgWidth);
        %         for a = 1:200
        %             avgImg = avgImg + double(imread([dataPath dataFile{stack}],a,'Info',fileInfo)) - darkAvg;
        %             %    avgImg = avgImg + double(imread([dataPath dataFile],a)) - darkAvg;
        %             %Deleted 'Info' -AC 6/21
        %         end
        %         avgImg = avgImg/200;
        
        %define variables related to the templates
%         if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
%             templateInfo = imfinfo(templateFile);
%             if templateInfo(1).Height ~= templateInfo(1).Width
%                 error('Template is not square');
%             end
%             templateSize = templateInfo(1).Height;
%         else
            load(threshFile,'template');
            templateSize = size(template,2);
%         end
        numTemplates = length(templateFrames);
        templateColors = jet(numTemplates);
        
        % make sure ROI is an even number of pixels: should also be done in
        % f_calSMidentification
        if mod(ROI(3),2)==1
            ROI(3) = ROI(3)-1;
        end
        if mod(ROI(4),2)==1
            ROI(4) = ROI(4)-1;
        end
        cropWidth = ROI(3);
        cropHeight = ROI(4);
        
        %% trace out FoV if computing laser profile
        
        if findLaserInt
            temp = inputdlg({'What was the laser power at the objective? (in mW)'},...
                'Input laser power',...
                1,...
                {'0'});
            powerAtObjective = str2double(temp{1})/1000;
        end
        
        if findLaserInt == 1;
            avgImg = zeros(imgHeight,imgWidth);
            avgImgFrames = min(200,length(frames));
%             im = tiffread2([dataPath dataFile{stack}], 1, avgImgFrames);
            for a = 1:avgImgFrames
%                 avgImg = avgImg + double(im(a).data) - darkAvg;
                avgImg = avgImg + double(imread([dataPath dataFile{stack}],frames(a),'Info',fileInfo)) - darkAvg;
            end
            avgImg = avgImg/avgImgFrames;
%             clear im
            
            if ~exist('FOVmask')
            hFOVmaskFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            imagesc(avgImg(ROI(2):ROI(2)+ROI(4)-1, ...
                ROI(1):ROI(1)+ROI(3)-1),[0 300]);
            axis image;colorbar;colormap hot;
            title('Trace out the field of view boundaries by clicking at vertices');
            FOVmask = roipoly;
            close(hFOVmaskFig);
            end
        end
        
        %% prepare template for template matching
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
%             
%             if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
%                 templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
%                     [(cropHeight-size(template,2))/2 ...
%                     (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
%             else
                templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
%             end
            
            % multiplying by conjugate of template in FT domain is equivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
        
        
        
    end %%% end init
    
    
    %% Identify frames to analyze based on the sequence log
    % load data and register sequence log to data frames
    
    if ~isequal(logPath,0)
        
        if strcmp(logFile{stack}(end-2:end),'dat')
            logType = 1; % frames are in sifLog format
            
            if length(logFile) == length(dataFile) % if one sif file/data file
                sifLogData =  importdata([logPath logFile{stack}]);
                sifLogData = sifLogData(1:numFrames,:);
            else % i.e., if only one sif file for all data files
                sifLogData =  importdata([logPath logFile{1}]);
                sifLogData = sifLogData(absFrameNum:absFrameNum+numFrames-1,:);
                absFrameNum = absFrameNum + numFrames;
            end
            if filterCPs   %%% Option 1: Select only stationary steps (for CP)
                % works by selecting periods between the '-1' markers that signify motion
                % if the length of the period != step dwell time (i.e., = translation period), skip
                %
                motionFrames = find(sifLogData(:,1)==-1);  % This finds -1 entries in the shutters correspoding to moving frames
                validPeriod = diff(motionFrames)==maxNumMeasurement+1;
                validFrames = [];
                validMotionFrames = motionFrames(validPeriod);
                beginScan = vertcat(false, validPeriod(2:end-1) & ~validPeriod(1:end-2) & validPeriod(3:end), false);
                beginScan(2) = validPeriod(2) && validPeriod(3);
                beginScanFrames = motionFrames(beginScan)+3;
                for i = 1:length(validMotionFrames)
                    initFrame = validMotionFrames(i);
                    temp = (initFrame+3:initFrame+maxNumMeasurement)';
                    validFrames = [validFrames; temp];
                end
                frames = intersect(validFrames,frames);
            else   %%% Option 2: Select only frames that are illuminated (SMACM 2-color)
                if channel == 'g' %&& ~filterCPs   % use an intersect in case frames are limited by user
                    frames = intersect(find(sifLogData(:,2) == 1),frames);
                elseif channel == 'r' %&& ~filterCPs
                    frames = intersect(find(sifLogData(:,3) == 1),frames);
                end
            end
            
        elseif strcmp(logFile{stack}(end-2:end), 'mat')
            % mat file output for control point scan using PI stage
            logType = 3;
            
            if filterCPs
                % load in relevant vars from mat file
                load([logPath logFile{stack}],'fVec_cp','nFrames_per_step','numCPScans');
                
                fVec_cp = fVec_cp .* -1; % flip sign
                sifLogData = zeros(length(fVec_cp).*nFrames_per_step,2);
                sifLogData(:,1) = 1:length(fVec_cp).*nFrames_per_step; %col1 is frameNum
                
                zVecRep = [];
                for ii=1:length(fVec_cp)
                    zVecRep = [zVecRep; repmat(fVec_cp(ii),[nFrames_per_step,1])];
                end
                sifLogData(:,2) = zVecRep-mean(fVec_cp); %fVec_cp is in microns
                sifLogData(:,2) = sifLogData(:,2).*1000; % convert to nm
                
                %%%%%%%%%%%%%
                beginScanFrames = [];
                beginScanCurr = 1;
                for ii= 1:numCPScans
                    beginScanFrames = [beginScanFrames; beginScanCurr];
                    beginScanCurr = beginScanCurr + length(sifLogData(:,1))./numCPScans;
                end
                
                validFrames = sifLogData(:,1);
                frames = intersect(sifLogData(:,1),frames);
                %%%%%%%%%%%%
                clear beginScanCurr zVecRep
                
            end
        end % end check what type of extension
    
    end % end check to make sure user loaded sif log
    
    %% if NHA data, choose points to fit
    if filterCPs && nhaData
        holeLocsAll = cell(length(beginScanFrames),1);
        dataSeed = [];
        holeLocsSeed = nan(4,2);
        tryCorrelation = 1;
        for idx = 1:length(beginScanFrames)
            chosenFrame = beginScanFrames(idx)+3;%%%%
%             im = tiffread2([dataPath dataFile{stack}],chosenFrame);
%             data = double(im.data) - darkAvg;
%             clear im
            data = double(imread([dataPath dataFile{stack}],chosenFrame,'Info',fileInfo))-darkAvg;
            data = data(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1); 
                % define array of points in NHA from user for this region
            if isequal(idx,1) || ~tryCorrelation
                dataSeed = data;
                hPicks=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
                imagesc(data); axis image;colorbar;colormap hot; 
                title('Select four fiducials that form a large box. Press mousewheel to eliminate specific holes.');
                hold on;
                % User will click on beads in image and then text will be imposed on the
                % image corresponding to the bead number.
                % moleLocs is a n by 2 array -- x y values for each bead
                holeLocs = [];
                erase = ones(0,2);
                n = 0;
                % Loop, collecting bead locations and drawing text to mark them
                but = 1;
                while but == 1
                    [xi,yi,but] = ginput(1);
                    if isempty(xi)
                        break
                    end
                    if but == 2;
                        erase = [erase; round([xi yi])];
                        text(xi,yi,'*','color','red','fontsize',13,'fontweight','bold');
                        but = 1;
                        continue
                    end
                    n = n+1;
                    text(xi,yi,num2str(n),'color','white','fontsize',13,'fontweight','bold');
                    holeLocs(n,:) = round([xi yi]);
                end
                hold off;
%                 angledNHA = 1; % assuming this is true, can be looped in from f_calSMidentificaiton
        %         edgeVerts = bwboundaries(FOVmask);
        %         edgeVerts = [edgeVerts{1}(:,2),edgeVerts{1}(:,1)];
                load(calFile,'absLocs');
                chullCals = convhull(absLocs(:,1),absLocs(:,2));
                acceptPoints = [absLocs(chullCals,1),absLocs(chullCals,2)]/nmPerPixel;
                holeLocsSeed = holeLocs;
                close(hPicks)
            else
                dataCorr = xcorr2(data,dataSeed);
                [~,maxIdx] = max(abs(dataCorr(:)));
                [yPk, xPk] = ind2sub(size(dataCorr),maxIdx(1));
                corrOffset = [xPk - size(data,2), yPk - size(data,1)];
                holeLocs = bsxfun(@plus,holeLocsSeed,corrOffset);
            end
            holeLocs = f_fill_NHA(holeLocs,nmPerPixel,2500,angledNHA);
            holeLocs = holeLocs(inpolygon(holeLocs(:,1),holeLocs(:,2),acceptPoints(:,1)-ROI(1)+1,acceptPoints(:,2)-ROI(2)+1),:);
            if erase
                for eraseIdx = 1:size(erase,1)
                    [~,toErase] = min((holeLocs(:,1)-erase(eraseIdx,1)).^2+(holeLocs(:,2)-erase(eraseIdx,2)).^2);
                    holeLocs(toErase,:) = [];
                end
            end
            
            hPicks=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            imagesc(data);axis image;colorbar;colormap hot; 
            for n = 1:length(holeLocs)
                text(holeLocs(n,1),holeLocs(n,2),num2str(n),'color','white','fontsize',13,'fontweight','bold');
            end
            saveas(hPicks,[outputFilePrefix{stack} 'Picked points frame ' num2str(chosenFrame) '.png']);
            close(hPicks);
            holeLocsAll{idx} = holeLocs;
        end
    end
    %% initialize figure and variables and start fitting
    
    
    hSMFits=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    %totalPSFfits=[frame#, loc# w/in frame, xLocation yLocation matchingTemplateNumber matchConfidence (6),
    %amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2 bkgndMean(15)
    %totalFitError goodFit xCenter yCenter angle numPhotons
    %interlobeDistance, amplitude ratio, sigma ratio(24)
    %corrected X, corrected Y, corrected Z (27)]    
    totalPSFfits = zeros(20000, 6+18+3);
    numPSFfits = 0;
    startTime = tic;
    logFlag = 0;
    cropWidth = ROI(3);
    cropHeight = ROI(4);
    bkgndImgTotal = zeros(length(ROI(2):ROI(2)+ROI(4)-1),...
        length(ROI(1):ROI(1)+ROI(3)-1));
    numbkgndImg = 0;
    
    if size(frames,1) > 1 % make sure it will work in the for loop (need 1xn)
        frames = frames';
    end
    
    % set up median filter
    dataWindow = nan([2*floor(windowSize/2)+1,size(bkgndImgTotal)]);
    if medianBlurSigma~=0
        medBlurFilt = fspecial('gaussian',100,medianBlurSigma);
    else
        medBlurFilt=1;
    end
    lastFrame = nan;
    
    % set up mean BG for laser spot calculation
    meanBG = nan(length(frames),1);
    meanSignal = nan(length(frames),1);
    xavg = ROI(3)-9:ROI(3);
    yavg = ROI(4)-9:ROI(4);
    if ~medianFilter&&~nhaData
        disp('finding background as average of all frames, may take a while')
        avgImg = zeros([length(yavg),length(xavg)]);
        for a = frames
%             im = tiffread2([dataPath dataFile{stack}], a);
%             newData = double(im.data) - darkAvg; clear im
            newData = double(imread([dataPath dataFile{stack}],a,'Info',fileInfo)) - darkAvg;
            avgImg = avgImg + newData(yavg,xavg);
        end
        avgImg = avgImg/length(frames);
    end
    
    for c=frames
        if isempty(frames) % if [] has been selected for this file
            break
        end
        
        % import data from image and filter appropriately
        currIdx = find(frames==c); % current index within 'frames'
%         im = tiffread2([dataPath dataFile{stack}],c);
%         data = double(im.data) - darkAvg; clear im
        data = double(imread([dataPath dataFile{stack}],c,'Info',fileInfo))-darkAvg;
        data = data(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);        % crop data to ROI
        if nhaData 
            data=data.*(FOVmask1);
        end
        if usePolyROI
            dataRing=data(FOVmask&~FOVmask1);
            data(~FOVmask)=median(dataRing);
            % the below can smooth the data near the edge of the FOV -
            % would be useful if using wavelet bg subtraction, e.g.
%             for i = 1:10
%             dataBlur=imfilter(data,gaussFilt,'replicate');
%             data(~FOVmask)=dataBlur(~FOVmask);
%             end
%             clear dataBlur
        end
        % subtract the background and continue
        if medianFilter
            
            if isnan(lastFrame) || currIdx==interpVal*round(currIdx/interpVal);
                [bkgndImgMed, dataWindow] = f_medianFilter([dataPath dataFile{stack}], fileInfo, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                lastFrame = c;
            end
            
            bkgndImg = bkgndImgMed;
            
            
            if medianBlurSigma~=0
                bkgndImg(~FOVmask)=median(reshape(bkgndImg(~FOVmask1&FOVmask),[],1));
                bkgndImg = imfilter(bkgndImg,medBlurFilt,'replicate');
            end
        elseif nhaData
            % use prctile of pixel intensity as BG estimate (not that
            % important in this context / hard to estimate given NHA fill
            % fraction)
            bgAmpEst = prctile(data(FOVmask),20);
            bkgndImg = bgAmpEst*ones(size(data));%median(data(:)).*ones(size(data));
            clear bgAmpEst
        else
             bkgndImg = ones(size(data))*mean(avgImg(:));
%             if c==frames(1)
%                 disp('set to no background subtraction!');
%             end
%             bkgndImg = zeros(size(data)); 
%             bkgndImg = f_waveletBackground(data);
        end
        
        if usePolyROI && medianFilter % otherwise the displayed image is wrong
            bkgndImg(~FOVmask) = data(~FOVmask);
        end
        
        bkgndImgTotal = bkgndImgTotal + bkgndImg;
        numbkgndImg = numbkgndImg +1;
        
        data = data - bkgndImg;
        
        % begin template matching in earnest
        
    
                       
        dataFT = fft2(data,cropHeight,cropWidth);
        maxPeakImg = zeros(cropHeight,cropWidth);
        % matrix PSFLocs stores information about double helices that were
        % found via template matching
        % rows are different matches
        % [xLocation yLocation matchingTemplateNumber matchConfidence];
        PSFLocs = zeros(100,4);
        numPSFLocs = 0;
        for b=1:numTemplates
            % try no prefiltering
            %H = 1;
            % try phase correlation
            %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
            % try weighted phase correlation (emphasizing low frequency
            % components
            % dataFT had zero amplitude 7/20/2016 when analyzing simulated
            % data and using a constant background => div-by-0 error
            H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
            H(isnan(H)) = 0;
            if nhaData
                H = gaussianFilter;
            end
            % normalize H so it doesn't add any energy to template match
            %H = H / sqrt(sum(abs(H(:)).^2));
            
            peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));
%             peakImg = abs(peakImg);%%%%%%%%%%%
            if exist('FOVmask') % prevents finding matches outside FOV (due to e.g. wrapping around)
                peakImg=peakImg.*FOVmask;
            end
            % normalize response of peakImg by dividing by number of pixels in
            % data
            %peakImg = peakImg / (cropHeight*cropWidth);
            maxPeakImg = max(maxPeakImg, peakImg);
            
            %threshold = mean(peakImg(:))+peakThreshold*std(peakImg(:));
            peakImg(peakImg < peakThreshold(stack,b)) = peakThreshold(stack,b);
            
            if isreal(peakImg) && sum(sum(isnan(peakImg)))==0
                temp = find(imregionalmax(peakImg));
            else
                % Write log file
                fileID = fopen('peakImg log.txt','a');
                fprintf(fileID,[datestr(now) '\r\n' dataFile{stack}]);
                fprintf(fileID,'\r\nROI: [%d %d %d %d]',ROI);
                fprintf(fileID,'\r\nFrame: %d\r\npeakImg matrix:\r\n',a);
                dlmwrite('peakImg log.txt',peakImg,'-append','delimiter','\t','newline','pc')
                fprintf(fileID,'\r\n********************NEXT********************\r\n\r\n');
                fclose('all');
                
                logFlag = logFlag + 1;
                peakImg(isnan(peakImg)) = 0; %inserted to deal with NaNs -AC 6/22
                temp = find(imregionalmax(real(peakImg)));
            end
            
            
            % make sure threshold didn't eliminate all peaks and create
            % lots of matches
            if length(temp) < cropHeight*cropWidth/2;
                [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                    [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                numPSFLocs = numPSFLocs+length(temp);
            end
        end
        clear H dataFT peakImg
        
        %% filter out extraneous matches due to very strong signals
        
        if numPSFLocs > 0
            % sort location matrix in decending order of confidence
            temp = sortrows(PSFLocs(1:numPSFLocs,:),-4);
            % copy most confident match to list of locations
            PSFLocs(1,:) = temp(1,:);
            numPSFLocs = 1;
            for b=2:size(temp,1)
                % make sure that this candidate location is a minimum distance away
                % from all other candidate locations
                if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
                    % add it to list of locations
                    numPSFLocs = numPSFLocs + 1;
                    PSFLocs(numPSFLocs,:) = temp(b,:);
                end
            end
        end
        
        % alternate way of fitting, ensure get correct points/angle
    if ~filterCPs || ~nhaData
            totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,1:6) = ...
            [c*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
    else
        
%         if ismember(c,beginScanFrames) %% CLIP
%             % define array of points in NHA from user for this region
%         hLocs=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%         imagesc(data); axis image;colorbar;colormap hot; 
%         title('Select four fiducials that form a large box. Press mousewheel to eliminate specific holes.');
%         hold on;
%         % User will click on beads in image and then text will be imposed on the
%         % image corresponding to the bead number.
%         % moleLocs is a n by 2 array -- x y values for each bead
%         holeLocs = [];
%         erase = ones(0,2);
%         n = 0;
%         % Loop, collecting bead locations and drawing text to mark them
%         but = 1;
%         while but == 1
%             [xi,yi,but] = ginput(1);
%             if isempty(xi)
%                 break
%             end
%             if but == 2;
%                 erase = [erase; round([xi yi])];
%                 text(xi,yi,'*','color','red','fontsize',13,'fontweight','bold');
%                 but = 1;
%                 continue
%             end
%             n = n+1;
%             text(xi,yi,num2str(n),'color','white','fontsize',13,'fontweight','bold');
%             holeLocs(n,:) = round([xi yi]);
%         end
%         hold off;
%         angledNHA = 1;
% %         edgeVerts = bwboundaries(FOVmask);
% %         edgeVerts = [edgeVerts{1}(:,2),edgeVerts{1}(:,1)];
%         load(calFile,'absLocs');
%         chullCals = convhull(absLocs(:,1),absLocs(:,2));
%         acceptPoints = [absLocs(chullCals,1),absLocs(chullCals,2)]/nmPerPixel;
%         holeLocs = f_fill_NHA(holeLocs,nmPerPixel,2500,angledNHA);
%         holeLocs = holeLocs(inpolygon(holeLocs(:,1),holeLocs(:,2),acceptPoints(:,1),acceptPoints(:,2)),:);
%         if erase
%             for eraseIdx = 1:size(erase,1)
%                 [~,toErase] = min((holeLocs(:,1)-erase(eraseIdx,1)).^2+(holeLocs(:,2)-erase(eraseIdx,2)).^2);
%                 holeLocs(toErase,:) = [];
%             end
%         end
%         close(hLocs)
%         hLocs=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%             imagesc(data(ROI(2):ROI(2)+ROI(4)-1, ...
%             ROI(1):ROI(1)+ROI(3)-1));axis image;colorbar;colormap hot; 
%         for n = 1:length(holeLocs)
%             text(holeLocs(n,1),holeLocs(n,2),num2str(n),'color','white','fontsize',13,'fontweight','bold');
%         end
% 
%         numPoints = size(holeLocs,1);
% %         totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,1:6) = ...
% %             [c*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
%         
%         end %% CLIP
            startFrameIdx = find(c>=beginScanFrames,1,'last'); %%%%
            holeLocs = holeLocsAll{startFrameIdx};
            numPoints = size(holeLocs,1);
            % use previous array of points found
            % [xLocation yLocation matchingTemplateNumber matchConfidence];
            phaseCorrLocs = PSFLocs(1:numPSFLocs,:);
            numPSFLocs = numPoints;
            PSFLocs = [holeLocs, repmat([mode(phaseCorrLocs(:,3)),max(phaseCorrLocs(:,4))],numPoints,1)];
            totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,1:6) = ...
            [c*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];

    end % end alternative NHA fitting method
        %% do fitting to extract exact locations of DH-PSFs
        
        % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2 bkgndMean(9)
        %  totalFitError goodFit xCenter yCenter angle numPhotons interlobeDistance amplitude ratio sigma ratio]
        PSFfits = zeros(numPSFLocs, 18);
        % create reconstructed DH-PSF image from fitted data
        reconstructImg = zeros(cropHeight, cropWidth);
        
        for b=1:numPSFLocs
            %% prepare parameters for fine fitting of PSFs
            
            % create indices to use for fitting
            [xIdx, yIdx] = meshgrid(PSFLocs(b,1)-boxRadius:PSFLocs(b,1)+boxRadius, ...
                PSFLocs(b,2)-boxRadius:PSFLocs(b,2)+boxRadius);
            % make sure indices are inside ROI
            if min(xIdx(:)) < 1
                xIdx = xIdx + (1-min(xIdx(:)));
            end
            if max(xIdx(:)) > cropWidth
                xIdx = xIdx - (max(xIdx(:))-cropWidth);
            end
            if min(yIdx(:)) < 1
                yIdx = yIdx + (1-min(yIdx(:)));
            end
            if max(yIdx(:)) > cropHeight
                yIdx = yIdx - (max(yIdx(:))-cropHeight);
            end
            
            % compute initial parameters from the location of two spots in
            % the templates
            fitParam(3) = PSFLocs(b,1) + templateLocs(PSFLocs(b,3),1)-(templateSize/2+0.5);
            fitParam(4) = PSFLocs(b,2) + templateLocs(PSFLocs(b,3),2)-(templateSize/2+0.5);
            fitParam(5) = PSFLocs(b,1) + templateLocs(PSFLocs(b,3),3)-(templateSize/2+0.5);
            fitParam(6) = PSFLocs(b,2) + templateLocs(PSFLocs(b,3),4)-(templateSize/2+0.5);
            % make sure initial guess lies within the ROI: if not, move on
            if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                    || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                    || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                    || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                PSFfits(b,11) = -1000;
                continue;
            end
            %  fitParam(1) = data(round(fitParam(4)),round(fitParam(3)))-bkgndMean;
            %  fitParam(2) = data(round(fitParam(6)),round(fitParam(5)))-bkgndMean;
            fitParam(1) = data(round(fitParam(4)),round(fitParam(3)));
            fitParam(2) = data(round(fitParam(6)),round(fitParam(5)));
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
            lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
                sigmaBounds(1) sigmaBounds(1)];
            %             upperBound = [max(max(data(yIdx(:,1),xIdx(1,:))))-bkgndMean ...
            %                 max(max(data(yIdx(:,1),xIdx(1,:))))-bkgndMean ...
            %                 max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
            %                 sigmaBounds(2) sigmaBounds(2)];
            upperBound = [max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
                sigmaBounds(2) sigmaBounds(2)];
            
            %% Fit with lsqnonlin
            
            %             [fitParam,temp,residual,exitflag] = lsqnonlin(@(x) ...
            %                 doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),bkgndMean,xIdx,yIdx),...
            %                 fitParam,lowerBound,upperBound,options);
            [fitParam,~,residual,exitflag] = lsqnonlin(@(x) ...
                f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                fitParam,lowerBound,upperBound,options);
            
            if findTrueSigma
                lowerBoundSig = [sigSearchBounds(1) sigSearchBounds(1)];
                upperBoundSig = [sigSearchBounds(2) sigSearchBounds(2)];
                [fitParamSigma,~,~,~] = lsqnonlin(@(x) ...
                f_doubleGaussianVector([fitParam(1:6) x],data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                fitParam(7:8),lowerBoundSig,upperBoundSig,sigOptions);
                fitParam(7:8) = fitParamSigma;
            end
           
            
            fittedBkgndMean = mean2(bkgndImg(yIdx(:,1),xIdx(1,:)))*conversionFactor;
            
            PSFfits(b,1:11) = [fitParam fittedBkgndMean sum(abs(residual)) exitflag];
            
            %% compute derived paramters from fine fitting output
            
            % shift coordinates relative to entire dataset (not just ROI)
            PSFfits(b,3) = PSFfits(b,3) + ROI(1)-1;
            PSFfits(b,4) = PSFfits(b,4) + ROI(2)-1;
            PSFfits(b,5) = PSFfits(b,5) + ROI(1)-1;
            PSFfits(b,6) = PSFfits(b,6) + ROI(2)-1;
            % Calculate midpoint between two Gaussian spots
            % shift coordinates relative to entire dataset (not just ROI) and
            % convert from pixels to nm
            PSFfits(b,12) = ((fitParam(3)+fitParam(5))/2 + ROI(1)-1)*nmPerPixel;
            PSFfits(b,13) = ((fitParam(4)+fitParam(6))/2 + ROI(2)-1)*nmPerPixel;
            
            % Below is the calculation of the angle of the two lobes.
            % Remember that two vertical lobes is focal plane because camera
            % outputs data that is rotated. Therefore, we want y2>y1 for all
            % angle calculations (so that -90<=angle<=90, and we use swap
            % the use of x and y for the atan2 calculation.
            x1 = fitParam(3);
            x2 = fitParam(5);
            y1 = fitParam(4);
            y2 = fitParam(6);
            % swap if y1>y2
            if (y1 > y2)
                tx = x1; ty = y1;
                x1 = x2; y1 = y2;
                x2 = tx; y2 = ty;
                clear tx ty;
            end
            %Finds the angle
            PSFfits(b,14) = atan2(-(x2-x1),y2-y1) * 180/pi;
            clear x1 x2 y1 y2;
            
            %Below is a way to count the photons. It integrates the box and
            %subtracts the boxarea*offset from the fit. It is inherently flawed
            %if there happen to be bright pixels inside of the fitting region.
            totalCounts = sum(sum(data(yIdx(:,1),xIdx(1,:))));
            %  PSFfits(b,15) = (totalCounts-(2*boxRadius+1)^2*bkgndMean)*conversionFactor;  % bkgndMean = 0
            PSFfits(b,15) = totalCounts*conversionFactor;  % bkgndMean = 0
            
            %The interlobe distance
            lobeDist = sqrt((fitParam(3)-fitParam(5)).^2 + ...
                (fitParam(4)-fitParam(6)).^2);
            PSFfits(b,16) = lobeDist;
            
            %Amplitude Ratio
            ampRatio = abs(fitParam(1) - fitParam(2))/sum(fitParam(1:2));
            PSFfits(b,17) = ampRatio;
            
            % Gaussian width Ratio
            sigmaRatio = abs(fitParam(7) - fitParam(8))/sum(fitParam(7:8));
            PSFfits(b,18) = sigmaRatio;
            
            %% Now evaluate the goodness of the fits
            
            % Conditions for fits (play with these):
            % (1) Amplitude of both lobes > 0
            % (2) All locations x1,y1, x2,y2 lie inside area of small box
            % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
            % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
            % (5) Make sure amplitudes are within 100% of one another
            % (6) Make sure totalFitError/(total number of photons) < 1.05
            
            if exitflag > 0
                % absolute amplitude > 0?
                if fitParam(1)<0 || fitParam(2)<0
                    PSFfits(b,11) = -1001;
                end
                % peaks inside box?
                if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                        || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                        || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                        || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                    PSFfits(b,11) = -1002;
                end
                % absolute sigma size for either lobe within bounds?
                if findTrueSigma && (...
                        fitParam(7)-0.001<=sigSearchBounds(1) || fitParam(8)-0.001<=sigSearchBounds(1) ...
                        || fitParam(7)+0.001>=sigSearchBounds(2) || fitParam(8)+0.001>=sigSearchBounds(2) )
                    PSFfits(b,11) = -1003;
                elseif ~findTrueSigma && (...
                        fitParam(7)<=sigmaBounds(1) || fitParam(8)<=sigmaBounds(1) ...
                        || fitParam(7)>=sigmaBounds(2) || fitParam(8)>=sigmaBounds(2) )
                    PSFfits(b,11) = -1003;
                end
                % sigma ratio of lobes less than limit?
                if sigmaRatio > sigmaRatioLimit;
                    PSFfits(b,11) = -1004;
                end
                % interlobe distance within bounds?
                if lobeDist < lobeDistBounds(1) || lobeDist > lobeDistBounds(2)
                    PSFfits(b,11) = -1005;
                end
                % amplitude ratio of lobes less than limit?
                if ampRatio > ampRatioLimit;
                    PSFfits(b,11) = -1006;
                end
                % normalized error within limit?
                if PSFfits(b,10)*conversionFactor/PSFfits(b,15) > 4.0  || ...
                        PSFfits(b,10)*conversionFactor/PSFfits(b,15) < 0.0
                    PSFfits(b,11) = -1007;
                end
                
            end
            
            % if the fit is good, add it to the reconstructed image
            if PSFfits(b,11) > 0
                [xIdx, yIdx] = meshgrid(1:cropWidth,1:cropHeight);
                reconstructImg = reconstructImg + ...
                    fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                    +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
            end
        end
        
        totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,6+1:6+18) = PSFfits;
        numPSFfits = numPSFfits+numPSFLocs;
        
        %%  plot results of template matching and fitting
        
        if dispRaw
            dataToPlot = data+bkgndImg;
            reconstToPlot = reconstructImg + bkgndImg;
            dataTitle = ['Frame ' num2str(c) ': raw data - darkAvg counts'];
        else
            dataToPlot = data;
            reconstToPlot = reconstructImg;
            dataTitle = ['Frame ' num2str(c) ': raw data - darkAvg & bkgnd'];
        end
        set(0,'CurrentFigure',hSMFits);
        subplot('Position',[0.025 0.025 .85/3 .95]);
        imagesc(maxPeakImg,[0 3*min(peakThreshold(stack,:))]);axis image;
        title({'Peaks correspond to likely template matches' ...
            [num2str(numPSFLocs) ' matches found']});
        
        subplot('Position',[0.075+.85/3 0.025 .85/3 .95]);
        imagesc(dataToPlot);axis image;colormap hot;
        hold on;
        for b=1:numPSFLocs
            %plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
            %    'MarkerSize', 15*PSFLocs(b,4)/peakThreshold(b), ...
            %    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
            plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                'MarkerSize', 15*PSFLocs(b,4)/peakThreshold(stack,PSFLocs(b,3)), ...
                'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
            %Changed peakThreshold(b) to peakThreshold(PSFLocs(b,3)) because b
            %does not seem to logically correspond to the correct template, whereas the third
            %column of PSFLocs is defined as corresponding to a specific
            %template. Furthermore, whenever numPSFLocs > length(peakThreshold)
            %there is an error. -AC 6/22
        end
        hold off;
        title({dataTitle ...
            ['ROI [xmin ymin width height] = ' mat2str(ROI)]});
        
        subplot('Position',[0.125+2*.85/3 0.025 .85/3 .95]);
        %         imagesc(reconstructImg+bkgndMean,[min(data(:)) max(data(:))]);axis image;
        imagesc(reconstToPlot,[min(dataToPlot(:)) max(dataToPlot(:))]);axis image;
        title({'Image reconstructed from fitted matches' ...
            [num2str(sum(PSFfits(:,11)>0)) ' successful fits']});
        
        drawnow;
        if printOutputFrames == 1
            set(gcf,'PaperPositionMode','auto');
            saveas(hSMFits, ['output images' filesep 'frame ' num2str(c) '.tif']);
        end
        
        meanSignal(currIdx) = mean(data(FOVmask));
        meanBG(currIdx) = mean(bkgndImg(FOVmask));
    end
    elapsedTime = toc(startTime);
    totalPSFfits = totalPSFfits(1:numPSFfits,:);
    clear data bkgnd residual fileInfo maxPeakImg reconstructImg xIdx yIdx temp;
    %     clear data bkgnd residual fileInfo maxPeakImg reconstructImg templateFT xIdx yIdx temp;
    
    if logFlag ~= 0
        sprintf('There were %d frames in which the peakImg matrix contained complex numbers or NaNs. See log file "peakImg log.txt" for more details',logFlag)
        logFlag = 0;
    end
    close(hSMFits)
    fps = length(frames)/elapsedTime
    moleculesPerSec = numPSFfits/elapsedTime
    
    %movie2avi(M,'output_v1.avi','fps',10,'Compression','FFDS');
    
    %% Measure the Gaussian Laser Intensity Distribution
    if findLaserInt == 1;
        % This takes the average bkgndImg found using f_waveletBackground
        % and tries to extract a Gaussian laser profile from it, assuming the
        % background intensity is uniformly proportional to laser intensity
        bkgndImg_avg = bkgndImgTotal/numbkgndImg;
        
        [laser_x_nm, laser_y_nm ,sigma_x_nm, sigma_y_nm, theta, peakIntensity, waist]...
            = f_findGaussianLaserProfile...
            (bkgndImg_avg, FOVmask, nmPerPixel, powerAtObjective, ROI);
        
    end
    
    %% Translate angle into corrected x,y positions and z position
    % set values as nan if beyond range of calibration
    load(calFile);
    goodFit_forward = logical(squeeze(goodFit_f(1,calBeadIdx,:)));
    totalPSFfits(:,25) = totalPSFfits(:,18) ...
        - interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),...
        squeeze(meanX(1,calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline',nan);
    totalPSFfits(:,26) = totalPSFfits(:,19) ...
        - interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),...
        squeeze(meanY(1,calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline',nan);
    totalPSFfits(:,27) = interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),...
        squeeze(z(1,calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline',nan);
    
    %% output data to external file
    
    totalPSFfits_columns = {'frame number', 'localization number in frame',...
        'x position (pixel)','y position (pixel)','template matched','template match strength',...
        'lobe1 amplitude','lobe2 amplitude','lobe1 x pos (pix)','lobe1 y pos (pix)',...
        'lobe2 x pos (pix)', 'lobe2 y pos (pix)', 'lobe1 sigma (pix)', 'lobe2 sigma (pix)',...
        'mean background (photons/pix)','fit error','fit flag','x position (nm)',...
        'y position (nm)','DH-PSF angle (deg)','total deteted photons (- bg)','interlobe distance (pix)',...
        'lobe amplitude ratio','lobe sigma ratio','calibration-corrected x (nm)',...
        'calibration-corrected y (nm)','z from calibration (nm)'};
    
    save([outputFilePrefix{stack} 'molecule fits.mat']);
    
    
end
end
