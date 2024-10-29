%% Do a custom transformation of the points
function [tform, cp_channel2_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
    cp_channel1, cp_channel2, tform_mode, nEqn,...
    weight_type, kthNeighbor, smoothnessParameter, figures, toCheck)
% This function computes a custom 2D/3D transformation using the given set of control points
% and also computes the associated FRE and TRE values.
% This is from "Single-Molecule High Resolution Colocalization of Single
% Probes" by L. Stirling Churchman and James A. Spudich in Cold Spring
% Harbor Protocols (2012), doi: 10.1101/pdb.prot067926
% Modified Definitions of FRE and TRE - Andreas Gahlmann, 20110511
% Added custom locally weighted 2D and 3D transformation that allow greater
% flexibility in choosing parameters.


% Evaluate the transformation funtion
if isequal(tform_mode, 'globalquadratic')  
    
    keep = randi(size(cp_channel1,1),[nEqn,1]);    
       
    tform = cp2tform_AG(...
        cp_channel1(keep,:),cp_channel2(keep,:),tform_mode, nEqn, ...
        weight_type, kthNeighbor, smoothnessParameter);
    
    trans_cp_channel2 = transformData(cp_channel2(keep,:),tform);
        scatter3(cp_channel1(keep,1), cp_channel1(keep,2), cp_channel1(keep,3))
        title({'Reference Channel';'Channel 1'})
        hold on
        scatter3(trans_cp_channel2(:,1), trans_cp_channel2(:,2), trans_cp_channel2(:,3), 10, 'filled')
        hold off
        CPdeviation = trans_cp_channel2(:,:) - cp_channel1(keep,:);
    
else
    
    tform = cp2tform_AG(...
        cp_channel1,cp_channel2,tform_mode, nEqn, ...
        weight_type, kthNeighbor, smoothnessParameter);
    
end


% Calculate the metrices to estimate the error associated with this
% Calculate the fiducial registration error (FRE)
trans_cp_channel2 = transformData(cp_channel2,tform);

FRE_full = [sqrt(sum((cp_channel1-trans_cp_channel2).^2,2)),...
    (cp_channel1-trans_cp_channel2)];
FRE = mean(FRE_full(:,1));
% hist(FRE_full, 30)

% Calculate the target registration error (TRE)
number_cp = length(cp_channel1); % find the number of control points
% Loop through the control points
if exist('toCheck')&&length(toCheck) == 1 % toCheck can specify a number of random points to choose
    numCutoff = toCheck;
end

% set up which CP to test TRE for. Default is all, then a random set of
% N=toCheck if toCheck is a scalar. Otherwise use idx specified by toCheck
if ~exist('toCheck') || (exist('numCutoff') && number_cp <= numCutoff)
    toCheck = 1:number_cp;
elseif exist('numCutoff') && number_cp > numCutoff
    toCheck = randperm(number_cp,numCutoff);    
end
% make sure all members of toCheck will be represented in set
if any(~ismember(toCheck,intersect(toCheck,1:number_cp)))
    warning('looks like you are trying to access something out of range')
end

tic
check_trans_cp_channel2 = nan(length(toCheck),3);
parfor j=1:length(toCheck)%i=1:number_cp
% for i=1:number_cp
    i = toCheck(j); % parfor must be increasing consecutive, use dummy
    remain_cp = [1:i-1 i+1:number_cp]; % take out that control point
    
    % Calculate the transformation without the ith control point
    tform = cp2tform_AG(cp_channel1(remain_cp,:),cp_channel2(remain_cp,:),...
        tform_mode, nEqn, weight_type, kthNeighbor, smoothnessParameter);
        
    % Transform left out control point with above found transformation
%     trans_cp_channel2(j,:) = transformData(cp_channel2(i,:),tform);
   check_trans_cp_channel2(j,:) = transformData(cp_channel2(i,:),tform);
end
toc

% TRE_full = [sqrt(sum((cp_channel1-trans_cp_channel2).^2,2)),...
%     (cp_channel1-trans_cp_channel2)];
TRE_full = [sqrt(sum((cp_channel1(toCheck,:)-check_trans_cp_channel2).^2,2)),...
    (cp_channel1(toCheck,:)-check_trans_cp_channel2)];
TRE = mean(TRE_full(TRE_full(:,1)<100000,1));

% Restore the full transform function again
tform = cp2tform_AG(cp_channel1,cp_channel2,tform_mode, nEqn, ...
    weight_type, kthNeighbor, smoothnessParameter);
cp_channel2_trans = transformData(cp_channel2,tform);

if figures
    figure
    distlimit = 30;
    
    if size(cp_channel1,2)==2
        subplot(2,2,1)
        scatter(cp_channel1(:,1), cp_channel1(:,2))
        title({'Reference Channel';'Channel 1'})
        hold on
        scatter(cp_channel2_trans(:,1), cp_channel2_trans(:,2), 10, 'filled')
        hold off
        axis square
        subplot(2,2,2)
        scatter(cp_channel2(:,1), cp_channel2(:,2))
        title({'Target Channel';'Channel 2'})
        axis square
    elseif size(cp_channel1,2)==3
        subplot(2,2,1) % reference channel + tformed points
        scatter3(cp_channel1(:,1), cp_channel1(:,2), cp_channel1(:,3));
        xCenter = (max(cp_channel1(:,1))+min(cp_channel1(:,1)))/2;
        yCenter = (max(cp_channel1(:,2))+min(cp_channel1(:,2)))/2;
        xRad = (max(cp_channel1(:,1))-min(cp_channel1(:,1)))/2;
        yRad = (max(cp_channel1(:,2))-min(cp_channel1(:,2)))/2;
        xlim(xCenter+[-xRad xRad]);
        ylim(yCenter+[-yRad yRad]);
        title({'Reference Channel';'Channel 1'})
        hold on
        scatter3(cp_channel2_trans(:,1), cp_channel2_trans(:,2), cp_channel2_trans(:,3), 10, 'filled')
        hold off
        subplot(2,2,2) % target channel only
        scatter3(cp_channel2(:,1), cp_channel2(:,2), cp_channel2(:,3));
        xCenter = (max(cp_channel2(:,1))+min(cp_channel2(:,1)))/2;
        yCenter = (max(cp_channel2(:,2))+min(cp_channel2(:,2)))/2;
        xRad = (max(cp_channel2(:,1))-min(cp_channel2(:,1)))/2;
        yRad = (max(cp_channel2(:,2))-min(cp_channel2(:,2)))/2;
        xlim(xCenter+[-xRad xRad]);
        ylim(yCenter+[-yRad yRad]);
        title({'Target Channel';'Channel 2'})
    end
    
    subplot(2,2,3)
    hist(FRE_full(FRE_full(:,1)<=distlimit,1), 30)
    title({['Fiducial Registration Error']; [tform_mode ' Transformation']});
    xlabel('Distance (nm)');
    ylabel('Frequency');
    xlim([0 distlimit]);
    legend(['Mean = ' num2str(FRE, 3) ' nm']);
    
    subplot(2,2,4)
    hist(TRE_full(TRE_full(:,1)<distlimit,1),30)
    title({['Target Registration Error']; [tform_mode ' Transformation']});
    xlabel('Distance (nm)');
    ylabel('Frequency');
    xlim([0 distlimit]);
    legend(['Mean = ' num2str(TRE, 3) ' nm']);
    disp(num2str(TRE, 3))
%     saveas(gcf,['Transform_' tform_mode '.fig']);
%     saveas(gcf,['Transform_' tform_mode '.png']);
%     close(gcf); 
end

end
