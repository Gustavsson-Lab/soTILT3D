% easy3Dcc
% easy 3D cross-correlation
% find the X Y Z to move one set of XYZ's to match the other set of XYZ's
% started on 21 April 2017

%% Place the XYZ's into separate vectors
X1 = S5(:,2);
Y1 = S5(:,3);
Z1 = S5(:,4);

X2 = S6(:,2);
Y2 = S6(:,3);
Z2 = S6(:,4);

%% Histogram the Z values
figure
plot(min(Z1):100:max(Z1),hist(Z1,min(Z1):100:max(Z1)),'r')
hold on
plot(min(Z2):100:max(Z2),hist(Z2,min(Z2):100:max(Z2)),'b')
hold off
legend('Z1','Z2')

%% Use lsqnonlin
maxZ1 = max(Z1);
options = optimset('TolFun',1e-6, 'MaxFunEvals', 3000, 'MaxIter', 3000,'Display','off','FunValCheck','on');

lb = -500;
ub = +500;
x0 = +50;
thicknessofslice = 200; % this takes the top 200 nm or so from the dataset with the lower Z positions
[Zcorrection,~,residual,exitflag] = lsqnonlin(@(x) ...
    myfuncc(x,[X1 Y1 Z1],[X2 Y2 Z2],thicknessofslice),x0,lb,ub,options);

%% Plot the after Z correction plots within a Z slice
Z1c = Z1+Zcorrection;
chosenZc_max = max(Z1c)-300; % this chooses the max Z value in the slice to plot
chosenZc_min = chosenZc_max - 500; % this chooses the min Z value in the slice to plot

figure

chosenZ1 = and((Z1c>chosenZc_min),(Z1c<chosenZc_max));
plot( X1(chosenZ1) , Y1(chosenZ1),'bo','MarkerSize',10)

chosenZ2 = and((Z2>chosenZc_min),(Z2<chosenZc_max));
hold on
plot( X2(chosenZ2) , Y2(chosenZ2),'rx','MarkerSize',10)
hold off

axis image
title('After Z correction')
legend('Z1','Z2')
% xlim([8600 10400])
% ylim([1.85E4 2.05E4])

%% Plot the before Z correction plots within a Z slice
figure

chosenZ1 = and((Z1>chosenZc_min),(Z1<chosenZc_max));
plot( X1(chosenZ1) , Y1(chosenZ1),'bo','MarkerSize',10)

chosenZ2 = and((Z2>chosenZc_min),(Z2<chosenZc_max));
hold on
plot( X2(chosenZ2) , Y2(chosenZ2),'rx','MarkerSize',10)
hold off

axis image
title('Before Z correction')
legend('Z1','Z2')
% xlim([8600 10400])
% ylim([1.85E4 2.05E4])











