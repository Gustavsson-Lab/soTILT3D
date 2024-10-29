function F = myfuncc(x,XYZ1,XYZ2,nmfromtop)

X1 = XYZ1(:,1);
Y1 = XYZ1(:,2);
Z1 = XYZ1(:,3);
X2 = XYZ2(:,1);
Y2 = XYZ2(:,2);
Z2 = XYZ2(:,3);
maxZ1 = max(Z1);
nearestneb = zeros(length(Z1),1);
parfor i = 1:length(XYZ1)
    thisX = X1(i);
    thisY = Y1(i);
    thisZ = Z1(i) + x;
    if abs(maxZ1-Z1(i)) < nmfromtop
        nearestneb(i) = min(pdist2([thisX thisY thisZ],[X2 Y2 Z2]));
    end
end
F = sum(nearestneb(:));