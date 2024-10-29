% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
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

function [templateFrames, ROI, dataFile, dataPath, darkFile, logFile,...
            logPath, EMGain, templateLocs, outputFilePrefix,nhaData] = f_calSMidentification(calFile,calBeadIdx,...
            templateFile, boxRadius,channel,sigmaBounds,gaussianFilterSigma,minDistBetweenSMs,horizontalDHPSF)
% f_calSMidentification is a module in easy_dhpsf that prepares the
% templates from f_calDHPSF and uses them to generate a series of template
% matches. These are then used to judge an appropriate threshold for 
% f_fitSMs. This module also sets the file and other parameters used for 
% f_fitSMs, as well as some parameters for f_trackFiducials.
% if this is true, select rect ROI for channel and poly for excluding fids,
% edge of FoV outside iris, any other impinging 'non-image' feature
usePolyROI = true;

% initializing variables for ROI selection
FOVmask = [];

% Instrument Specific Parameters

dlg_title = 'Set EM Gain';
prompt = {  'EM Gain (1 if no gain):' };
def = {'300'};

num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

if isempty(inputdialog)
    error('User cancelled the program')
end

EMGain = str2double(inputdialog{1});
if EMGain < 1 || isnan(EMGain)
    warning('EMGain should be >= 1. Setting to 1...');
    EMGain = 1;
end

frameNum = 1;
scrsz = get(0,'ScreenSize');
% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');

%% ask user for relevant datafiles

[dataFile, dataPath] = uigetfile({'*.tif';'*.*'},...
    'Open SMACM image stack(s) for data processing',...
    'MultiSelect', 'on');
if isequal(dataFile,0)
    error('User cancelled the program');
end

if ischar(dataFile)
    dataFile = cellstr(dataFile);
end

% allows user to limit the number of files in case there are many, many
% large files (i.e., very long acquisitions where the sample does not
% change very much)
if length(dataFile) > 1
    dlg_title = 'Select Files';
    prompt = {  'Choose files for thresholding' };
    def = {num2str(1:length(dataFile))};
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(inputdialog)
        error('User cancelled the program')
    end
    selectedFiles = str2num(inputdialog{1});
else
    selectedFiles = 1:length(dataFile);
end

% allows user to select subset of frames, at the beginning, for all files
dlg_title = 'Select Frames (default: use all frames)';
num_lines = 1;
def = {};
prompt = {};
% populates 'fileInfo' and 'numFrames' for all files, and generates the
% fields for the frame selection dlg
for i = 1:length(selectedFiles)
    fileInfoAll{i} = imfinfo([dataPath dataFile{selectedFiles(i)}]);
    numFramesAll(i) = length(fileInfoAll{i});
    def{i} = ['[1:' num2str(numFramesAll(i)) ']'];
    prompt{i} = ['Choose frames for ' dataFile{selectedFiles(i)}];
end
   
%%%% edit so that dialog box fits in screen
if length(selectedFiles) <= 10
    % small number of files, display all at once
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    allinputdialog = inputdialog;
else
    allinputdialog = {};
    for i =1:ceil(length(selectedFiles)/10)
        startFile = (i-1)*10 + 1; 
        if i == ceil(length(selectedFiles)/10)
            endFile = length(selectedFiles);
        else
            endFile = startFile + 9;
        end
        inputdialog = inputdlg(prompt(startFile:endFile), dlg_title, num_lines, def(startFile:endFile));
        allinputdialog = [allinputdialog; inputdialog];
    end
end

for i = 1:length(selectedFiles)
    framesAll{i} = str2num(allinputdialog{i});
end
    
for stack = selectedFiles
    
    fileIdx = find(selectedFiles == stack);
    fileInfo = fileInfoAll{fileIdx};
    frames = framesAll{fileIdx};
    numFrames = numFramesAll(fileIdx);
    
%     % lets user specify a smaller number of frames in case of
%     % large files with many SMs, since templates will be sampled quickly.
%     % For simplicity, only activated when only when there is only one file
%     % (makes batch processing work quickly without user input)
%     if length(selectedFiles) == 1
%         dlg_title = 'Select Frames';
%         prompt = {  'Choose frames for thresholding' };
%         def = {[ '[1:' num2str(numFrames) ']' ]};
%         num_lines = 1;
%         inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
% 
%         if isempty(inputdialog)
%             error('User cancelled the program')
%         end
%         frames = str2num(inputdialog{1});
%     else
%         frames = 1:numFrames;
%     end
    imgHeight = fileInfo(1).Height;
    imgWidth = fileInfo(1).Width;
    
    if stack == selectedFiles(1)
        
%         [templateFile, templatePath] = uigetfile({'*.tif;*.mat';'*.*'},'Open image stack or MATLAB *.mat file of a DH-PSF template');
%         if isequal(templateFile,0)
%             error('User cancelled the program');
%         end
        
        
%         if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
%             
%             templateInfo = imfinfo(templateFile);
%             if templateInfo(1).Height ~= templateInfo(1).Width
%                 error('Template is not square');
%             end
%             templateSize = templateInfo(1).Height;
%         else
            load(templateFile);
%             if nhaData
                  % an alternate way to do this is to just manually edit the template.mat files to cut out the neighboring PSFs  
%                 template=template(:,5:22,5:22); % crops out lobes in corners with high-frequency psfs (why not fork into main code?)
%             end
            templateSize = size(template,2);
%         end
        clear templateFrames 
        load(calFile);
        
        if ~exist('templateFrames'); % if already calculated
        goodFit_forward = logical(squeeze(goodFit_f(1,calBeadIdx,:)));
        % compute templates to use based upon fitted DG angles
        maxAngle = 10*floor(min([max(meanAngles(:)) 80])/10);
        minAngle = 10*ceil(max([min(meanAngles(:)) -80])/10); % round to 10 degrees
        templateFrames = interp1(squeeze(meanAngles(1,calBeadIdx,goodFit_forward)),...
            1:length(meanAngles(1,calBeadIdx,goodFit_forward)),linspace(minAngle,maxAngle,6),'nearest'); %90:-30:-60
        end
        % if some angles are out of range, use the ends of the template
        % stack. first determine whether frames are increasing or
        % decreasing so the the 'ends' are chosen correctly.
        
        % use 2nd frame since 1st seems to give wrong angle (at least for
        % NHAs)
        if any(isnan(templateFrames))
            if mean(diff(templateFrames)) >= 0
                endFrames = [2 sum(goodFit_forward)];
            elseif mean(diff(templateFrames),'omitnan') < 0
                endFrames = [sum(goodFit_forward) 2];
            else
                endFrames = nan;
            end
            if isnan(templateFrames(1))
                templateFrames(1) = endFrames(1);
            end
            if isnan(templateFrames(end))
                templateFrames(end) = endFrames(2);
            end
        end

        temp = inputdlg({['Input sets of frames corresponding to each template or the index of templates to use (ex. ' ...
            mat2str(templateFrames) ')']},...
            'Input template numbers',1, ...
            {mat2str(templateFrames)}); ...     % {'[8:4:32]'});
        templateFrames = str2num(temp{1});
        
%         logPath = 0;
%         logFile = 0;
        [logFile, logPath] = uigetfile({'*.dat';'*.*'},...
            'Open sequence log file(s) corresponding to image stack(s) (optional: hit cancel to skip)',...
            'MultiSelect', 'on');
        if ischar(logFile)
            logFile = cellstr(logFile);
        end

        [darkFile, darkPath] = uigetfile({'*.tif';'*.*'},'Open image stack with dark counts (same parameters as SMACM data)');
%         if isequal(darkFile,0)
%             error('User cancelled the program');
%         end
        
%inputdlg(prompt,title,nl,def,options)
         temp = inputdlg({'Use (M)edian filtering, (W)avelet filtering, treat as (N)HA?, or n(O) subtraction?',...
                          'If median filtering, what sigma size in pix (0 for no smoothing)',...
                          'If median filtering, what total window size in frames?',...
                          'If median filtering, number of frames to "interpolate"',...
                          'When fitting data, display raw (not background-subtracted) data?'},...
                'Input filtering options',1, ...
                {'M','15','101','10','1'});
        filterType = temp{1};
        if strcmp(filterType,'m')
            filterType = 'M';
        elseif strcmp(filterType,'w')
            filterType = 'W';
        elseif strcmp(filterType,'n')
            filterType = 'N';
        end
        
        medianBlurSigma = str2num(temp{2});
        windowSize = str2num(temp{3});
        interpVal = str2num(temp{4}); % number of frames to interpolate in median bg estimation. 1 uses each frame.
        dispRaw = logical(str2num(temp{5}));
        
        if strcmp(filterType,'N')
            nhaData = true;
            medianFilter = false;
            waveFilter = false;
        elseif strcmp(filterType,'M')
            medianFilter = true;
            nhaData = false;
            waveFilter = false;
        elseif strcmp(filterType,'W')
            medianFilter = false;
            nhaData = false;
            waveFilter = true;
        else
            medianFilter = false;
            nhaData = false;
            waveFilter = false;            
        end
    end
    
    %% create output log filenames
    % saves in labeled directory if a channel is selected
    if channel == '0'
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep 'threshold ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
    else
        outputFilePrefix{stack} = [dataPath dataFile{stack}(1:length(dataFile{stack})-4) filesep channel(1) ' threshold ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
    end
    
    mkdir(outputFilePrefix{stack});
    
    if stack== selectedFiles(1)
        %% Compute darkAvg counts
        
        if ~isequal(darkFile,0)
            darkFile = [darkPath darkFile];
            % Computes average of dark frames for background subtraction
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
                warning('Dark count image and data image stack are not the same size. Resizing darkAvg count image...');
                darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
            end
        else
            darkAvg = 0;
        end
        clear darkFileInfo;
        
 
        
        %% opens and processes templates of DH-PSF for template matching
        
%         if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
%             
%             numTemplates = size(templateFrames,1);
%             templateColors = jet(numTemplates);
%             template = zeros(numTemplates,templateSize,templateSize);
%             templateLocs = zeros(numTemplates,5);
%             fitParam = zeros(1,8);
%             [xIdx, yIdx] = meshgrid(1:templateSize,1:templateSize);
%             hTemplate=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%             set(hTemplate,'Visible','off');
%             for a=1:numTemplates
%                 for b=templateFrames(a,1):templateFrames(a,2)
%                     template(a,:,:) = squeeze(template(a,:,:)) + ...
%                         double(imread([templatePath templateFile],b,'Info',templateInfo));
%                 end
%                 % make minimum count level in template equal to 0
%                 template(a,:,:) = template(a,:,:) - min(min(template(a,:,:)));
%                 % normalize energy contained (sum of all counts) in the template
%                 template(a,:,:) = template(a,:,:) / sum(sum(template(a,:,:)));
%                 % finally, make mean of template equal to 0
%                 template(a,:,:) = template(a,:,:) - mean(mean(template(a,:,:)));
%                 
%                 % find two largest peaks in template
%                 [tempY, tempX] = ind2sub([templateSize templateSize],find(imregionalmax(template(a,:,:))));
%                 temp = sortrows([tempX tempY template(sub2ind([numTemplates templateSize templateSize],a*ones(length(tempX),1),tempY,tempX))],-3);
%                 
%                 % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
%                 fitParam(3) = temp(1,1);
%                 fitParam(4) = temp(1,2);
%                 fitParam(5) = temp(2,1);
%                 fitParam(6) = temp(2,2);
%                 fitParam(1) = temp(1,3);
%                 fitParam(2) = temp(2,3);
%                 fitParam(7) = 1.8;
%                 fitParam(8) = 1.8;
%                 lowerBound = [0 0 1 1 1 1 sigmaBounds(1) sigmaBounds(1)];
%                 upperBound = [max(max(template(a,:,:))) max(max(template(a,:,:))) ...
%                     templateSize templateSize templateSize templateSize ...
%                     sigmaBounds(2) sigmaBounds(2)];
%                 
%                 % Fit with lsqnonlin
%                 fitParam = lsqnonlin(@(x) ...
%                     f_doubleGaussianVector(x,squeeze(template(a,:,:)),0,xIdx,yIdx),...
%                     fitParam,lowerBound,upperBound,options);
%                 
%                 templateLocs(a,1:2) = fitParam(3:4);
%                 templateLocs(a,3:4) = fitParam(5:6);
%                 % calculate rough angle between peaks
%                 templateLocs(a,5) = 180/pi*atan2(templateLocs(a,2)-templateLocs(a,4), ...
%                     templateLocs(a,1)-templateLocs(a,3));
%                 
%                 subplot(1,numTemplates,a);imagesc(squeeze(template(a,:,:)));
%                 axis image;colormap hot;colorbar;
%                 hold on;
%                 plot(templateLocs(a,1),templateLocs(a,2),'.','MarkerEdgeColor', templateColors(a,:));
%                 plot(templateLocs(a,3),templateLocs(a,4),'.','MarkerEdgeColor', templateColors(a,:));
%                 title({['Template ' mat2str(templateFrames(a,:))] ...
%                     ['Angle = ' num2str(templateLocs(a,5)) ' deg']});
%             end
%             imwrite(frame2im(getframe(hTemplate)),[outputFilePrefix{stack} 'templates.tif']);
%             clear templateInfo tempX tempY temp xIdx yIdx;
% 
%         else
            templateFrames = templateFrames';
            numTemplates = size(templateFrames,1);
            templateColors = jet(numTemplates);
            templateLocs = zeros(numTemplates,5);
            fitParam = zeros(1,8);
            [xIdx, yIdx] = meshgrid(1:templateSize,1:templateSize);
            
            
            hTemplate=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            for a=1:numTemplates
                
                % make minimum count level in template equal to 0
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    - min(min(template(templateFrames(a),:,:)));
                % normalize energy contained (sum of all counts) in the template
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    / sum(sum(template(templateFrames(a),:,:)));
                % finally, make mean of template equal to 0
                if ~nhaData
                    template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                        - mean(mean(template(templateFrames(a),:,:)));
                end
                
                % find two largest peaks in template
                hFilt = fspecial('gaussian', 3, 0.5);
                blurTempl = imfilter(squeeze(template(templateFrames(a),:,:)),hFilt);
                [tempY, tempX] = ind2sub([templateSize templateSize],find(imregionalmax(blurTempl)));
                temp = sortrows([tempX tempY template(sub2ind(size(template),...
                    templateFrames(a)*ones(length(tempX),1),tempY,tempX))],-3);
                
                % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
                fitParam(3) = temp(1,1);
                fitParam(4) = temp(1,2);
                fitParam(5) = temp(2,1);
                fitParam(6) = temp(2,2);
                fitParam(1) = temp(1,3);
                fitParam(2) = temp(2,3);
                fitParam(7) = mean(sigmaBounds);
                fitParam(8) = mean(sigmaBounds);
                lowerBound = [0 0 1 1 1 1 sigmaBounds(1) sigmaBounds(1)];
                upperBound = [max(max(template(templateFrames(a),:,:))) max(max(template(templateFrames(a),:,:))) ...
                    templateSize templateSize templateSize templateSize ...
                    sigmaBounds(2) sigmaBounds(2)];
                
                % Fit with lsqnonlin
                fitParam = lsqnonlin(@(x) ...
                    f_doubleGaussianVector(x,squeeze(template(templateFrames(a),:,:)),0,xIdx,yIdx),...
                    fitParam,lowerBound,upperBound,options);
                
                
                x1 = fitParam(3);
                x2 = fitParam(5);
                y1 = fitParam(4);
                y2 = fitParam(6);
                if (~horizontalDHPSF && (y1 > y2) )  ||  (horizontalDHPSF && x1 > x2)
                    tx = x1; ty = y1;
                    x1 = x2; y1 = y2;
                    x2 = tx; y2 = ty;
                    clear tx ty;
                end
                %Finds the angle

                templateLocs(a,1:2) = fitParam(3:4);
                templateLocs(a,3:4) = fitParam(5:6);
                
                % calculate rough angle between peaks
                
                if ~horizontalDHPSF
                    templateLocs(a,5) = atan2(-(x2-x1),y2-y1) * 180/pi;
                elseif horizontalDHPSF
                    templateLocs(a,5) = atan2(-(y2-y1),x2-x1) * 180/pi;
                end
                
%                 templateLocs(a,5) = 180/pi*atan2(templateLocs(a,2)-templateLocs(a,4), ...
%                     templateLocs(a,1)-templateLocs(a,3));
                
                subplot(1,numTemplates,a);
                imagesc(squeeze(template(templateFrames(a),:,:)));
                axis image;colormap hot;colorbar;
                hold on;
                plot(templateLocs(a,1),templateLocs(a,2),'.','MarkerEdgeColor', templateColors(a,:));
                plot(templateLocs(a,3),templateLocs(a,4),'.','MarkerEdgeColor', templateColors(a,:));
                title({['Frames ' mat2str(templateFrames(a,:))] ...
                    ['Angle = ' num2str(templateLocs(a,5)) ' deg']});
                %         drawnow
                %         pause(1)
            end
            imwrite(frame2im(getframe(hTemplate)),[outputFilePrefix{stack} 'templates.tif']);
            clear templateInfo tempX tempY temp xIdx yIdx;

%         end
        close(hTemplate); % closes template figure
        %% user picks ROI
        % pick region of interest by reading first frame and having user select
        % region
        
        % Compute average image
        avgImg = zeros(imgHeight,imgWidth); 
        avgImgFrames = min(200,length(frames));
%         im = tiffread2([dataPath dataFile{stack}], 1,avgImgFrames);
        for a = 1:avgImgFrames
%             avgImg = avgImg + double(im(frames(a)).data) - darkAvg;
            avgImg = avgImg + double(imread([dataPath dataFile{stack}],frames(a),'Info',fileInfo)) - darkAvg;
        end
        avgImg = avgImg/avgImgFrames;
%         clear im
        
        hROI = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(avgImg,[0 min(avgImg(:))+2*std(avgImg(:))]);axis image;colormap hot;
        if channel == 'g'
            ROI = imrect(gca,[1 1 270 270]);
        elseif channel == 'r'
            ROI = imrect(gca,[243 243 270 270]);
        else
            ROI = imrect(gca,[1 1 size(avgImg,1) size(avgImg,2)]);
        end
        
        % ROI = imrect(gca,[1 1 128 128]);
        title({'Shape box and double-click to choose region of interest for PSF extraction' ...
            ['[xmin ymin width height] = ' mat2str(ROI.getPosition)]...
            ['The displayed image is the average of the first ' num2str(avgImgFrames) ' frames']});
        addNewPositionCallback(ROI,@(p) title({'Shape box and double-click to choose region of interest for PSF extraction' ...
            ['[xmin ymin width height] = ' mat2str(p,3)]...
            'The displayed image is the average of the first 200 frames'}));
        % make sure rectangle stays within image bounds
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(ROI,fcn);
        ROI = round(wait(ROI));
        % make sure ROI is an even number of pixels
        if mod(ROI(3),2)==1
            ROI(3) = ROI(3)-1;
        end
        if mod(ROI(4),2)==1
            ROI(4) = ROI(4)-1;
        end
        %ROI = [84 127 128 130];
        cropWidth = ROI(3);
        cropHeight = ROI(4);
        
            
        close(hROI) % closes ROI selection
        
        if usePolyROI
            hFOVmaskFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            imagesc(avgImg(ROI(2):ROI(2)+ROI(4)-1, ...
                ROI(1):ROI(1)+ROI(3)-1),[0 min(avgImg(:))+2*std(avgImg(:))]);
            axis image;colorbar;colormap hot;
            title('Select ROIpoly of area to keep');
            [FOVmask, maskX, maskY] = roipoly;
            xCenter=(max(maskX)+min(maskX))/2;
            yCenter=(max(maskY)+min(maskY))/2;
            x1=(maskX-xCenter)*0.95+xCenter;
            y1=(maskY-yCenter)*0.95+yCenter;
            FOVmask1=roipoly(avgImg(ROI(2):ROI(2)+ROI(4)-1, ...
                ROI(1):ROI(1)+ROI(3)-1),x1,y1);
            close(hFOVmaskFig);
        end
        
        % a possible alternate way to do this automatically:
%         blah=avgImg(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
%         for i = 5:5:95
%             datapoint(i) = prctile(blah(:),i);
%         end
%         thresh=(mean(datapoint([5 30]))-min(blah(:)))/(max(blah(:))-min(blah(:)));
%         bw=im2bw((blah-min(blah(:)))/(max(blah(:))-min(blah(:))),thresh);
%         figure;imagesc(bw.*blah,[0 20])
        %% prepare template for template matching
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
            
%             if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
%                 templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
%                     [(cropHeight-size(template,2))/2 ...
%                     (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
%             else
                templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
%             end
            
            % multiplying by conjugate of template in FT domain is squivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad temp;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
        
    end % end of the prep that is done only for first file.
    
    %% Identify frames to analyze when limiting by sif log

    if ~isequal(logPath,0)
%         if length(logFile) == length(dataFile)
%             sifLogData =  importdata([logPath logFile{stack}]);
%             sifLogData = sifLogData(1:numFrames,:);
%         else
%             sifLogData =  importdata([logPath logFile{1}]);
%             sifLogData = sifLogData(frameNum:frameNum+numFrames-1,:);
%             frameNum = frameNum + numFrames;
%         end
        if length(logFile) == length(dataFile)
            sifLogData =  importdata([logPath logFile{stack}]);
            if size(sifLogData,2) > 2
                logType = 1;
                sifLogData = sifLogData(1:numFrames,:);
            else
                logType = 3; % zscan done with PI stage/MATLAB
                load([logPath logFile{1}],'fVec_cp','nFrames_per_step')
                
                fVec_cp = fVec_cp.*-1; % flip sign
                sifLogData = zeros(length(fVec_cp).*nFrames_per_step,2);
                sifLogData(:,1) = 1:length(fVec_cp).*nFrames_per_step;
                
                zVecRep = repmat(fVec_cp, [nFrames_per_step, 1]);
                zVecRep = zVecRep(:);
                sifLogData(:,2) = zVecRep-mean(fVec_cp); %fVec_cp is actual z pos (micron units)
                sifLogData(:,2) = sifLogData(:,2)*1000; % convert to nm
            end
        else
            sifLogData =  importdata([logPath logFile{1}]);
            sifLogData = sifLogData(frameNum:frameNum+numFrames-1,:);
            frameNum = frameNum + numFrames;
        end
        if logType == 1
            if channel == 'g'   % use an intersect in case frames are limited by user
                frames = intersect(find(sifLogData(:,2) == 1),frames);
            elseif channel == 'r'
                frames = intersect(find(sifLogData(:,3) == 1),frames);
            end
        elseif logType ==3
            frames = intersect(sifLogData(:,1),frames);
        end
    end
    %% do template matching
    
    hMatchFig = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    totalPSFfits = zeros(10000, 6+15+3);
    numPSFfits = 0;
    startTime = tic;
%     frameNum = 1;
    if size(frames,1) > 1 % make sure it will work in the for loop (need 1xn)
        frames = frames';
    end
    
    dataWindow = nan([2*floor(windowSize/2)+1,ROI(4),ROI(3)]);
    if medianBlurSigma~=0
        medBlurFilt = fspecial('gaussian',100,medianBlurSigma);
    else
        medBlurFilt=1;
    end
    lastFrame = nan;
    
    meanBG = nan(length(frames),1);
    meanSignal = nan(length(frames),1);
    

    if ~medianFilter&&~nhaData&&~waveFilter
        xavg = 45:64;
        yavg = 55:64;
        disp('finding background as average of all frames, may take a while')
        avgImg = zeros([length(yavg),length(xavg)]);
        for a = frames
            newData = double(imread([dataPath dataFile{stack}],a,'Info',fileInfo)) - darkAvg;
            avgImg = avgImg + newData(yavg,xavg);
        end
        avgImg = avgImg/length(frames);
    end
    
    for c = frames(end:-1:1)
        
        currIdx = find(frames==c);
        
%         im = tiffread2([dataPath dataFile{stack}], c);
%         data = double(im.data) - darkAvg;
        data = double(imread([dataPath dataFile{stack}],c,'Info',fileInfo))-darkAvg;
        data = data(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
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
%         clear im
        % subtract the background and continue
        
        if medianFilter
            
            % only run on first frame, or frames divisible by interpVal
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
            % use median of all pixels as BG estimate
            bkgndImg = median(data(:)).*ones(size(data));
        elseif waveFilter
            bkgndImg = f_waveletBackground(data);
        else
            bkgndImg = ones(size(data))*mean(avgImg(:));
%             bkgndImg = zeros(size(data));
        end
        
        if usePolyROI && medianFilter
            bkgndImg(~FOVmask) = data(~FOVmask);
        end
        
        data = data - bkgndImg;
        
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
            peakImg = peakImg.*FOVmask;
            % normalize response of peakImg by dividing by number of pixels in
            % data
            %peakImg = peakImg / (cropHeight*cropWidth);
            maxPeakImg = max(maxPeakImg, peakImg);
            
            % only remember matches that are 3 standard deviations above the
            % mean
            peakThreshold = mean(peakImg(:))+3*std(peakImg(:));
            if nhaData
%                 peakThreshold=peakThreshold/2; % account for not using phase
                peakThreshold = median(peakImg(:));%+std(peakImg(:));
            end
            peakImg(peakImg < peakThreshold) = peakThreshold;
            temp = find(imregionalmax(peakImg));
            % make sure threshold didn't eliminate all peaks and create
            % lots of matches
            if length(temp) < cropHeight*cropWidth/2;
                [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                    [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                numPSFLocs = numPSFLocs+length(temp);
            end
        end

        clear H dataFT peakImg tempX tempY temp
        
        %% filter out extraneous matches due to very strong signals
        numAllLocs = numPSFLocs;
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
        totalPSFfits(numPSFfits+1:numPSFfits+numPSFLocs,1:6) = ...
            [numAllLocs*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];

        %% output an example image for each threshold level
        %  so that user can pick appropriate threshold later
        
        for b=1:numPSFLocs
            moleThreshold = round(PSFLocs(b,4)*10000);
            if nhaData
                moleThreshold=round(PSFLocs(b,4)); %makes number more reasonable
            end
            moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') '.png']; %%f6.4 without scaling
            if isempty(dir([outputFilePrefix{stack} moleFileName]))
                % create indices to isolate image of candidate molecule
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
                
                % output a picture of the
                img = data(yIdx(:,1),xIdx(1,:));
                img = 1+round(255*(img-min(img(:)))/(max(img(:))-min(img(:))));
                imwrite(imresize(ind2rgb(img,hot(256)),3,'nearest'),[outputFilePrefix{stack} moleFileName]);
            end
        end
        numPSFfits = numPSFfits+numPSFLocs;
        
        
        %%  plot results of template matching and fitting
        
%         name = {'bkgndImg', 'data', 'data+bkgndImg'};
%         dataToDraw = {bkgndImg, data, data+bkgndImg};
%         for drawNum = 1:3
%             draw = dataToDraw{drawNum};
%             imwrite(ind2rgb(uint8(256*(draw-min(draw(:)))/(max(draw(:))-min(draw(:)))),hot(256)),[outputFilePrefix{1} name{drawNum} num2str(c) '.png'],'png')
%         end
%         if c/5 == round(c/5)
%             figure('Visible','off');
%             imagesc(bkgndImg); axis image; colormap hot; colorbar;
%             title(['background image for frame ' num2str(c)]);
%             saveas(gcf,[outputFilePrefix{1} 'exampleBG' num2str(c) '.png']);
%             close(gcf)
%             
%             figure('Visible','off');
%             imagesc(data+bkgndImg); axis image; colormap hot; colorbar;
%             title(['raw data for frame ' num2str(c)]);
%             saveas(gcf,[outputFilePrefix{1} 'exampleRaw' num2str(c) '.png']);
%             close(gcf)
%             
%             figure('Visible','off');
%             imagesc(data); axis image; colormap hot; colorbar;
%             title(['bgsub data for frame ' num2str(c)]);
%             saveas(gcf,[outputFilePrefix{1} 'exampleBGsub' num2str(c) '.png']);
%             close(gcf)
%             
%             waveBG = f_waveletBackground(data+bkgndImg);
%             
%             figure('Visible','off');
%             imagesc(waveBG); axis image; colormap hot; colorbar;
%             title(['wavelet BG for frame ' num2str(c)]);
%             saveas(gcf,[outputFilePrefix{1} 'exampleWaveBG' num2str(c) '.png']);
%             close(gcf)
%             
%             figure('Visible','off');
%             imagesc(data+bkgndImg-waveBG); axis image; colormap hot; colorbar;
%             title(['data - wavelet BG for frame ' num2str(c)]);
%             saveas(gcf,[outputFilePrefix{1} 'exampleWaveBGImg' num2str(c) '.png']);
%             close(gcf)
%         end
%         clear name dataToDraw drawNum draw
        
        set(0,'CurrentFigure',hMatchFig);
        subplot('Position',[0.025 0.025 .9/2 .95],'parent',hMatchFig);
        imagesc(maxPeakImg,[0 3*peakThreshold]);axis image;
        title({'Peaks correspond to likely template matches' ...
            [num2str(numPSFLocs) ' matches found']});
        
        subplot('Position',[0.525 0.025 .9/2 .95],'parent',hMatchFig);
        imagesc(data); axis image;colormap hot;
        hold on;
        for b=1:numPSFLocs
            plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                'MarkerSize', 15*PSFLocs(b,4)/peakThreshold, ...
                'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
        end
        hold off;
        title({['Frame ' num2str(c) ': raw data - bkgnd & dark offset'] ...
            ['ROI [xmin ymin width height] = ' mat2str(ROI)]});

        drawnow;
        
        meanSignal(currIdx) = mean(data(FOVmask));
        meanBG(currIdx) = mean(bkgndImg(FOVmask));
        
    end
    elapsedTime = toc(startTime);
    totalPSFfits = totalPSFfits(1:numPSFfits,:);
    clear data bkgnd residual fileInfo maxPeakImg reconstructImg xIdx yIdx temp;
    close(hMatchFig) % closes fitting figure to prevent messiness
    
    fps = length(frames)/elapsedTime
    moleculesPerSec = numPSFfits/elapsedTime
    
    %% output data to external files
    
    textHeader = {'frame number' 'molecule number' 'template x location in ROI (px)' 'template y location in ROI (px)' ...
        'matching template number' 'match confidence (au)' ...
        'amp 1 (counts)' 'amp 2 (counts)' 'x location 1 (px)' 'y location 1 (px)' ...
        'x location 2 (px)' 'y location 2 (px)' 'sigma 1 (px)' 'sigma 2 (px)' 'background mean (counts)' ...
        'total fit error (counts)' 'good fit flag' 'x center (nm)' 'y center (nm)' ...
        'angle (deg)' 'number of photons' 'aberration corrected x location (nm)' ...
        'aberration corrected y location (nm)' 'z location (nm)'};
    % save fit info to MATLAB mat file
    save([outputFilePrefix{stack} 'threshold output.mat']);

    
end
outputFilePrefix{1} = outputFilePrefix{selectedFiles(1)};
end