function [filledFids] = f_fill_NHA(fidBounds,nmPerPixel,holePitch,angled)

% n x 2 array of vertices that bound the selection
% nm per pixel in image space
% pitch in NHA coordinate space in nm
if ~exist('angled') % to make 'backwards-compatible'
    % if angled = false, the NHA is at 0 or 90 degrees to lab coords, else the NHA is 45 degrees to lab coords (i.e. 8A back)
    angled = true; 
end

if angled
pixSep = holePitch * sqrt(2) / nmPerPixel; % nominal distance between holes in pix
elseif ~angled
pixSep = holePitch / nmPerPixel; % distance between holes in NHA frame
end
xLoc=fidBounds(:,1);
yLoc=fidBounds(:,2);

cornerPoints = 1:4;

% define square grid in lab coordinates (thus separation is pitch*sqrt(2))
% may need to rotate this if NHA is not parallel to lab axes
[~,NEpoint] = max(xLoc(cornerPoints).^2+yLoc(cornerPoints).^2);
[~,NWpoint] = max(-xLoc(cornerPoints).^2+yLoc(cornerPoints).^2);
[~,SEpoint] = max(xLoc(cornerPoints).^2-yLoc(cornerPoints).^2);
[~,SWpoint] = max(-xLoc(cornerPoints).^2-yLoc(cornerPoints).^2);
% % corners = [NEpoint NWpoint SEpoint SWpoint];
corners = fidBounds([NEpoint;NWpoint;SEpoint;SWpoint],:);

NE = fidBounds(NEpoint,:);
NW = fidBounds(NWpoint,:);
SE = fidBounds(SEpoint,:);
SW = fidBounds(SWpoint,:);

% NE = [max(xLoc(cro),max(yLoc)];
% NW = [min(xLoc),max(yLoc)];
% SE = [max(xLoc),min(yLoc)];
% SW = [min(xLoc),min(yLoc)];

center = [max(xLoc)+min(xLoc),max(yLoc)+min(yLoc)]/2;

if length(xLoc)>4
xLoc = xLoc([5:end,5]); % exclude first four points, and close ring
yLoc = yLoc([5:end,5]);
chosenBound = true;
else
chosenBound = false;
    xLoc = center(1)+pixSep*10*[1 1 -1 -1]';
    yLoc = center(2)+pixSep*10*[1 -1 1 -1]';
end
xD=diff(xLoc); % generate diffs for later
yD=diff(yLoc);

boxWidth = round(sqrt(sum((NW-NE).^2))/pixSep); % width in holes (total #holes-1)
boxHeight = round(sqrt(sum((NW-SW).^2))/pixSep); % height in holes (total #holes-1)

rowStep = ( (NE-NW) + (SE-SW) )/2 / boxWidth;
colStep = ( (NW-SW) + (NE-SE) )/2  / boxHeight;

trueSep = mean([sqrt(sum((NW-NE).^2)) / boxWidth,sqrt(sum((NW-SW).^2))/boxHeight]);

startPos = [-ceil((SW(1)-min(xLoc))/trueSep),-ceil((SW(2)-min(yLoc))/trueSep)];
endPos = [ceil((max(xLoc)-SW(1))/trueSep),ceil((max(yLoc)-SW(2))/trueSep)];

[numRow,numCol] = meshgrid(startPos(1):endPos(1),startPos(2):endPos(2));
% fidLocs1: in line with corners
fidLocs1 = [numRow(1:end)*rowStep(1) + numCol(1:end)*colStep(1);...
           numRow(1:end)*rowStep(2) + numCol(1:end)*colStep(2)]'+repmat(SW,length(numRow(:)),1);
% [numRow2,numCol2] = meshgrid(0:boxWidth-1,0:boxHeight-1);
% fidLocs2 = [numRow2(1:end)*rowStep(1) + numCol2(1:end)*colStep(1);...
%            numRow2(1:end)*rowStep(2) + numCol2(1:end)*colStep(2)]'+repmat(SW+rowStep/2+colStep/2,length(numRow2(:)),1);
if angled
    fidLocs2(:,1) = fidLocs1(:,1) + rowStep(1)/2 + colStep(1)/2;
    fidLocs2(:,2) = fidLocs1(:,2) + rowStep(2)/2 + colStep(2)/2;
    fidLocs = round([fidLocs1; fidLocs2]);
elseif ~angled
    fidLocs = round(fidLocs1);
end
if ~chosenBound
    filledFids = fidLocs;
    return
end
% now select points *within* the ring of points chosen by the user
% to do so, find the set of points on the 'inside side' of each line
above = false(length(xLoc),length(fidLocs));
centerUp = false(length(xLoc),1);
within = false(length(fidLocs),1);

lineTol = 1000/nmPerPixel; % allow 1000 nm outside of 'lasso'

for i = 1:length(xLoc)-1 % loop over each line of perimeter
if xD(i) ~= 0 % test whether line is vertical
    % a point is above the line AB if Y-Ya > (Yb-Ya)/(Xb-Xa) (X-Xa)

    % to call more points above, subtract from line
    % to call more points below, add to line

    % to be more lenient wrt center position,
    % subtract if center is above and add if center is below

    % find whether 'center' is above or below the current line
    centerUp(i)=center(2)-yLoc(i)>(yD(i)/xD(i))*(center(1)-xLoc(i));

    angle = atan(yD(i)/xD(i));
    angleCorr = 1/cos(angle);
    above(i,:)=fidLocs(:,2)-yLoc(i)>(yD(i)/xD(i))*(fidLocs(:,1)-xLoc(i))...
                + angleCorr*(-centerUp(i)*lineTol + ~centerUp(i)*lineTol);
else
    centerUp(i) = center(1)>xLoc(i); % define up as right in case of vertical lines
    above(i,:) = fidLocs(:,1) > xLoc(i);
end % end if statement
end % end loop over lines
% find whether points are on same side as center
for i = 1:length(fidLocs)
within(i) = isequal(above(:,i), centerUp);
end


filledFids = fidLocs(within,:);

end

% for i = 1:length(xLoc-1)
% figure; plot(fidLocs(above(i,:)==centerUp(i),1),fidLocs(above(i,:)==centerUp(i),2),'go');
% hold on; plot(fidLocs(above(i,:)~=centerUp(i),1),fidLocs(above(i,:)~=centerUp(i),2),'ro');
% hold on; plot(xLoc(i:i+1),yLoc(i:i+1),'b')
% title(num2str(i))
% end

% figure; plot(fidLocs(:,1),fidLocs(:,2),'ob')
% axis equal
% hold on; plot(filledFids(:,1),filledFids(:,2),'og')
% hold on; plot(xLoc,yLoc,'or')
