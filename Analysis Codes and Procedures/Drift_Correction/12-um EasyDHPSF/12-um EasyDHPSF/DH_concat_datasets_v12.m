function DH_concat_datasets_v12()

% 2018-10-31
% v12 rewrote for no fid correction

%clear all;
close all;
numSyncFrames = 25;
useDenoising = 1;
moleFiles = {}; 

% load('C:\UserFiles\UsersNOTBackedUp\Temp\2014-04-09_rib-eYFP-P146\S1R1\S1R1-green.mat')
[easyFile, easyPath] = uigetfile({'*.mat';'*.*'}, 'Open easy DHPSF file')

load([easyPath easyFile]); 

fitFilePrefix = s.fitFilePrefix; 
fidFilePrefix = s.fidFilePrefix; 
channel = s.channel;
calFile = [s.calFilePrefix 'calibration.mat']; 
useFidCorrections = 0; 

%% Spatially dependent calibrations
%Construct a questdlg with three options
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

calFile1 = 0; 
if spatialCorr == true
%     [calFile1, calPath] = uigetfile({'*.mat';'*.*'}, 'Open calibration file'); 
    calFile = [s.calFilePrefix 'calibration.mat']; 
end

if ~isequal(calFile1,0)
    calFile = [calPath calFile1]
end

%% query for where to save data
[saveFilePrefix, savePath] = uiputfile({'*.*'},'Enter a prefix for the concatenated datafiles');
if isequal(saveFilePrefix,0)
    error('User cancelled the program');
end
saveFilePrefix = [savePath saveFilePrefix ' '];


if length(fitFilePrefix)>0

    % load datafiles:  
    for cc = 1:length(fitFilePrefix)

        load([fitFilePrefix{cc} 'molecule fits.mat']);
        
        moleFiles = [moleFiles; {[fitFilePrefix{cc}]}];

        if cc == 1
            tempPSFfits = totalPSFfits(:,1:27);
            numFramesInFiles = numFrames;
        else
            numFramesInFiles = [numFramesInFiles numFrames];
            totalPSFfits(:,1) = totalPSFfits(:,1) + sum(numFramesInFiles(1:cc-1));
            tempPSFfits = [tempPSFfits; totalPSFfits(:,1:27)];
        end
    end
    totalPSFfits = tempPSFfits;
    numFrames = sum(numFramesInFiles);

    % apply spatially dependent 
    if spatialCorr
        totalPSFfits = makeLocalCals(totalPSFfits,calFile,'SMACM',0,1);
    end

    
    clear tempPSFfits; 

    % save concatenated data
    save([saveFilePrefix 'molecule fits.mat']);
 
end

end % end main function 

function [PSFfits] = makeLocalCals(PSFfits,calFile,type,calBeadIdx,undoZangleZero)
%     interpType = 'linear';
    interpType = 'NN';
    load(calFile,'absLocs','goodFit_f','meanAngles','meanX','meanY','z','zAngleZero');
    numCals = size(absLocs,1);
    
    if strcmp(type,'SMACM')
        xCol = 18;
        yCol = 19;
        angCol = 20;
        xOut = 25;
        yOut = 26;
        zOut = 27;
        goodFitCol = 13;
    elseif strcmp(type,'FID')
        xCol = 14;
        yCol = 15;
        angCol = 16;
        xOut = 21;
        yOut = 22;
        zOut = 23;
        goodFitCol = 13;
    end
    
    goodCal = true(numCals,1);
    for i = 1:numCals 
        if sum(goodFit_f(1,i,:)) < 10 %&& sum(goodFit_b(1,i,:)) < 10
            calDists(:,i) = nan;
            goodCal(i) = false;
            disp(['fiducial # ' num2str(i) ' did not have enough good fits!']);
        end
    end    
    
    if numCals >= 3 && strcmp(type,'SMACM')
    goodX = absLocs(goodCal,1); goodY = absLocs(goodCal,2);
    calHull = convhull(goodX,goodY);
    calHull = [goodX(calHull) goodY(calHull)];
    withinHull = inpolygon(PSFfits(:,xCol),PSFfits(:,yCol),calHull(:,1),calHull(:,2));
    
    figure; plot(calHull(:,1),calHull(:,2))
    hold on; plot(PSFfits(withinHull,xCol),PSFfits(withinHull,yCol),'.g')
    hold on; plot(PSFfits(~withinHull,xCol),PSFfits(~withinHull,yCol),'.r')
    title('Well-defined points, that fall within calibrations');
    else
        withinHull = true(size(PSFfits,1),1);
    end
    
    
    if strcmp(interpType,'NN')
    % length(PSFfits) x numCals... is this too big?
    calDists = sqrt(...
                   (repmat(PSFfits(:,xCol),1,numCals)-repmat(absLocs(:,1)',length(PSFfits),1)).^2 +...
                   (repmat(PSFfits(:,yCol),1,numCals)-repmat(absLocs(:,2)',length(PSFfits),1)).^2);
    
    [~,calNN] = min(calDists,[],2);
    
    
    % this should be rewritten to be an interpolation - this is not a good
    % way to do it long-term
    nearnessScores = ones(size(calDists))./(calDists.^2);
    nearnessScores = nearnessScores./repmat(sum(nearnessScores,2),1,size(nearnessScores,2));

        
    goodFits = PSFfits(:,goodFitCol)>0;
    
    if exist('calBeadIdx') && calBeadIdx ~= 0
        calNN(:)=calBeadIdx;
    end
    %goodFit_forward = logical(squeeze(goodFit_f(1,calBeadIdx,:)));
    
    % for fids: this easily tests whether fid drifts between cals
    % figure; hist(calNN(goodFits,:),1:numCals)
    
    for i = 1:numCals % loop over each cal and redo points near that cal
        if ~goodCal(i)
            continue
        end
        goodFit_forward = logical(squeeze(goodFit_f(1,i,:)));
        PSFfits(calNN==i,xOut) = PSFfits(calNN==i,xCol) ...
            - interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
            squeeze(meanX(1,i,goodFit_forward)),PSFfits(calNN==i,angCol),'spline',nan);
        PSFfits(calNN==i,yOut) = PSFfits(calNN==i,yCol) ...
            - interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
            squeeze(meanY(1,i,goodFit_forward)),PSFfits(calNN==i,angCol),'spline',nan);
        PSFfits(calNN==i,zOut) = interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
            squeeze(z(1,i,goodFit_forward)),PSFfits(calNN==i,angCol),'spline',nan);
        if undoZangleZero
            PSFfits(calNN==i,zOut) = PSFfits(calNN==i,zOut) + zAngleZero(i);
        end
    end
    % end nearest-neighbor module
    elseif strcmp(interpType,'linear')
        % re-extract true observed z values; do not make z=0 at theta=0
        z = squeeze(z(1,1,:))+zAngleZero(1);
        % define vectors repeating across each calibration
        calXrpt = repmat(absLocs(:,1),1,length(z));
        calXrpt = calXrpt+squeeze(meanX); % correct for 'proper drift'
        calXrpt = reshape(calXrpt',[],1);
        
        calYrpt = repmat(absLocs(:,2),1,length(z));
        calYrpt = calYrpt+squeeze(meanY); % correct for 'proper drift'
        calYrpt = reshape(calYrpt',[],1);
        
        calZrpt = repmat(z,numCals,1);
        calTrpt = reshape(squeeze(meanAngles)',[],1);
        goodFitRpt = logical(reshape(squeeze(goodFit_f)',[],1));
        F = scatteredInterpolant(calXrpt(goodFitRpt),calYrpt(goodFitRpt),...
                                 calTrpt(goodFitRpt),calZrpt(goodFitRpt),...
                                 'linear','none');        
        PSFfits(:,zOut) = F(PSFfits(:,xOut),PSFfits(:,yOut),PSFfits(:,angCol));
        % no xy interpolation/correction yet!!
    end % end choice of interpolation modules 
    
    
    
    if exist('withinHull')
        PSFfits(~withinHull,[xOut,yOut,zOut]) = nan; % cannot assign these, since outside range of cals
    end

        
        %% if using nearness scores - be careful with this!!
        % Should be rewritten to be a more careful interpolation
%     PSFfits(:,zOut) = 0;
%     for i = 1:numCals
%         goodFit_forward = logical(squeeze(goodFit_f(1,i,:)));
%         PSFfits(:,zOut) = PSFfits(:,zOut) + nearnessScores(:,i).*interp1(squeeze(meanAngles(1,i,goodFit_forward)),...
%             squeeze(z(1,i,goodFit_forward)),PSFfits(:,angCol),'spline');
%         if undoZangleZero
%         PSFfits(:,zOut) = PSFfits(:,zOut) + nearnessScores(:,i)*zAngleZero(i);
%         end
%     end
end