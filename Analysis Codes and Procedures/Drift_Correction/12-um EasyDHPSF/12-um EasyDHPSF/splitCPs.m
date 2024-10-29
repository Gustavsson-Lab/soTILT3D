% this script allows the user to split up identified control points
% (i.e., output from Identify_CP_3D) by time or column. Run immediately
% before Evaluate_Transform.

% to double check that matched locs are indeed matched correctly...
% x = [matched_cp_transmitted(1:5,5) matched_cp_reflected(1:5,5)];
% y = [matched_cp_transmitted(1:5,6) matched_cp_reflected(1:5,6)];
% z = [matched_cp_transmitted(1:5,7) matched_cp_reflected(1:5,7)];
% figure; scatter(x(1:5),y(1:5),'r');
% hold on; scatter(x(6:10),y(6:10),'b')
% shg
% axis equal

clear all;
%close all;

[CPFile, CPPath] = uigetfile({'*.mat';'*.*'},'Open output from Identify_ControlPoints_3D');
if isequal(CPFile,0)
   error('User cancelled the program');
end

%% split fits into first half and second half in time
load([CPPath CPFile]);

% matched_cp arrays have fields (n/a means not useful, i.e. related to loc# within clusters, etc.):
% [N/A,N/A,frameStart,frameEnd,meanX,...
%  meanY,meanZ,stdX,stdY,stdZ,...
%  meanPhotons,numMeasurements,N/A]

refCP = matched_cp_reflected;
transCP = matched_cp_transmitted;

sizeRef = size(refCP,1);
sizeTrans = size(transCP,1);

refIdx = [true(sizeRef,1); false(sizeTrans,1)];
allCP = [refCP;transCP];

frames = allCP(:,3); % starting frame

earlyFits = frames <= median(frames);

goodRef = refIdx & earlyFits;
goodTrans = ~refIdx & earlyFits;
matched_cp_reflected = allCP(goodRef,:);
matched_cp_transmitted = allCP(goodTrans,:);

save([CPPath CPFile(1:length(CPFile)-4) ' first half.mat']);

goodRef = refIdx & ~earlyFits;
goodTrans = ~refIdx & ~earlyFits;

matched_cp_reflected = allCP(goodRef,:);
matched_cp_transmitted = allCP(goodTrans,:);

save([CPPath CPFile(1:length(CPFile)-4) ' second half.mat']);

%% split fits into random z-columns interspersed in time

if sizeRef ~= sizeTrans
    error('The following step will absolutely not work without correctly matched CP.');
end

goodFits1 = false(sizeRef,1);
goodFits2 = false(sizeRef,1);

xRef = refCP(:,5);
yRef = refCP(:,6);
xTrans = transCP(:,5);
yTrans = transCP(:,6);

odd = true;
for a=1:sizeRef
    if ~goodFits1(a) && ~goodFits2(a)
        inCol = sqrt((xRef-xRef(a)).^2 + (yRef-yRef(a)).^2) <= 100 &...
                sqrt((xTrans-xTrans(a)).^2 + (yTrans-yTrans(a)).^2) <= 100;
        if odd
            goodFits1 = goodFits1 | inCol;
        else
            goodFits2 = goodFits2 | inCol;
        end
        odd = ~odd;
    end
end
clear odd;

matched_cp_reflected = refCP(goodFits1,:);
matched_cp_transmitted = transCP(goodFits1,:);

save([CPPath CPFile(1:length(CPFile)-4) ' odd.mat']);

matched_cp_reflected = refCP(goodFits2,:);
matched_cp_transmitted = transCP(goodFits2,:);

save([CPPath CPFile(1:length(CPFile)-4) ' even.mat']);
