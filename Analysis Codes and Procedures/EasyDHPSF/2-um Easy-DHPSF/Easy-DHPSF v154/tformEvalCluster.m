%% NOTE! Need to fix keep = randperm to be same w/in all cores...

function [] = tformEvalCluster(iStart,iEnd,iStep,jStart,jEnd,jStep,kStart,kEnd,kStep,nCP,CPPathFile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Ask user for relevant datafiles
% 
% [CPFile, CPPath] = uigetfile({'*.mat';'*.*'},'Open data file');
% if isequal(CPFile,0)
%     error('User cancelled the program');
% end
load(CPPathFile);

% iStart = 18;
% iEnd = 80;
% iStep = 2;
% jStart = 5;
% jEnd = 11;
% jStep = 1;
% kStart = 0.5;
% kEnd = 1.3;
% kStep = 0.2;
% temp=parcluster;        % query local machine for number of usable cores
% nCores = temp.NumWorkers-1; % set default to number of cores - 1
% clear temp;
% nCP = 1000;

% dlg_title = 'User input required';
% prompt = {'Would you like to evaluate the  parameter ranges?'};
% def =       { 'No'  };
% questiondialog = questdlg(prompt,dlg_title, def);
% % Handle response
% switch questiondialog
%     case 'Yes'
%         anotherpass = true;
%     case 'No'
%         anotherpass = false;
%     case 'Cancel'
%         error('User cancelled the program');
% end


% Specific Parameters

% while anotherpass == true
    
%     dlg_title = 'Please input range for parameters';
%     prompt = {  'Number of controlpoints start', ...
%         'Number of controlpoints end', ...
%         'Number of controlpoints stepsize', ...
%         'Weight function range (sigma to nearest neighbor) start', ...
%         'Weight function range (sigma to nearest neighbor) end', ...
%         'Weight function range (sigma to nearest neighbor) stepsize', ...
%         'Global weight fucntion scale factor start', ...
%         'Global weight fucntion scale factor end'...
%         'Global weight fucntion scale factor stepsize',...
%         'How many CPU cores? (Leave one core unused, if you want to still use this machine for other stuff',...
%         'How many control points used for FRE and TRE evaluation'...
%         };
%     def = {     num2str(iStart), ...
%         num2str(iEnd), ...
%         num2str(iStep), ...
%         num2str(jStart), ...
%         num2str(jEnd), ...
%         num2str(jStep), ...
%         num2str(kStart), ...
%         num2str(kEnd), ...
%         num2str(kStep), ...
%         num2str(nCores), ...
%         num2str(nCP) ...
%         };
%     num_lines = 1;
%     inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
%     
%     iStart = str2double(inputdialog{1});
%     iEnd = str2double(inputdialog{2});
%     iStep = str2double(inputdialog{3});
%     jStart = str2double(inputdialog{4});
%     jEnd = str2double(inputdialog{5});
%     jStep = str2double(inputdialog{6});
%     kStart = str2double(inputdialog{7});
%     kEnd = str2double(inputdialog{8});
%     kStep = str2double(inputdialog{9});
%     nCores = str2double(inputdialog{10});
%     nCP = str2double(inputdialog{11});
    
    %% Find the optimal transformation parameters
    
%     step = floor(size(matched_cp_transmitted,1)/1000);
%     keep = 1:step:size(matched_cp_transmitted,1);     % to accelerate the computation, reduce number of points
%     keep = randi(size(matched_cp_transmitted,1),[nCP,1]);   % this generated duplicates      
    keep = randperm(size(matched_cp_transmitted,1), nCP)'; 
%     figure
%     scatter3(matched_cp_reflected(keep,5), matched_cp_reflected(keep,6), matched_cp_reflected(keep,7))
%     title('Inspect averaged bead positions in reflected channel. Press any key to continue')
%     pause
    keep = 1:nCP; %% obviously change this!
    clear FRE TRE FRE_full TRE_full
    
    totalSteps = length(iStart:iStep:iEnd)*...
        length(jStart:jStep:jEnd)*...
        length(kStart:kStep:kEnd)
    
    s=0;
    a=1; b=1; c=1; % indexes progress over parameter space
    % matlabpool('open', nCores)
    
    for i = iStart:iStep:iEnd
        for j = jStart:jStep:jEnd
            for k = kStart:kStep:kEnd
                
                startTime = tic;
                [~, ~, FRE(a,b,c), TRE(a,b,c), ~, ~] = custom_transformation(...
                    matched_cp_reflected(keep,5:7),matched_cp_transmitted(keep,5:7),...
                    'lwquadratic',i,'Gaussian',j,k,false);
                elapsedTime = toc(startTime)
                
                estimatedComputationTimeHours = elapsedTime*totalSteps/3600
                
                c = c + 1;
                s = (s + 1);
                t = s/(((iEnd-iStart)/iStep+1)*((jEnd-jStart)/jStep+1)*((kEnd-kStart)/kStep+1))
                
            end
            c = 1;
            b = b + 1;
        end
        b = 1;
        a = a + 1;
    end
%    matlabpool close
    elapsedTime = toc(startTime)
    
    % j = 5:19;
    % k = 0.2:0.2:2.0;
    
    
    [~,ind]=min(TRE(:));
    [ind1,ind2,ind3]=ind2sub(size(TRE),ind)
    
    save('FRE.mat', 'FRE')
    save('TRE.mat', 'TRE')
    save('EvaluateTransform_3D_workspace.mat');
    disp('done with first stage, yo')
    %%
    
%     showAgain = true;
%     while showAgain == true
%         
%         load('EvaluateTransform_3D_workspace.mat');
%         load('FRE', 'FRE')
%         load('TRE', 'TRE')
%         
%         j = jStart:jStep:jEnd;
%         k = kStart:kStep:kEnd;
%         
%         for i = 1:a-1
%             contours = 40;
%             subplot(1,2,1)
%             contourf(k,j,squeeze(FRE(i,:,:)),min(min(squeeze(FRE(i,:,:)))):...
%                 (max(max(squeeze(FRE(i,:,:))))-min(min(squeeze(FRE(i,:,:)))))/contours:...
%                 max(max(squeeze(FRE(i,:,:)))))
%             title(num2str((i-1)*iStep+iStart))
%             colorbar
%             
%             subplot(1,2,2)
%             contourf(k,j,squeeze(TRE(i,:,:)),min(min(squeeze(TRE(i,:,:)))):...
%                 (max(max(squeeze(TRE(i,:,:))))-min(min(squeeze(TRE(i,:,:)))))/contours:...
%                 max(max(squeeze(TRE(i,:,:)))))
%             colorbar
%             
%             drawnow
%             waitforbuttonpress
%             %     pause(2)
%         end
%         
%         dlg_title = 'Show again';
%         prompt = {'Would you like to see the output again?'};
%         def =       { 'Yes'  };
%         questiondialog = questdlg(prompt,dlg_title, def);
%         % Handle response
%         switch questiondialog
%             case 'Yes'
%                 showAgain = true;
%             case 'No'
%                 showAgain = false;
%             case 'Cancel'
%                 error('User cancelled the program');
%         end
%     end
    
    
%     %% Construct a questdlg with three options
%     
% %     % f = figure;
% %     h = uicontrol('Position',[20 20 200 40],'String','Continue',...
% %         'Callback','uiresume(gcbf)');
% %     % disp('This will print immediately');
% %     uiwait(gcf);
% %     % disp('This will print after you click Continue');
% %     %     close(f);
%     
%     dlg_title = 'Evaluate Again';
%     prompt = {'Would you like to reevaluate with a different parameter range?'};
%     def =       { 'No'  };
%     questiondialog = questdlg(prompt,dlg_title, def);
%     % Handle response
%     switch questiondialog
%         case 'Yes'
%             anotherpass = true;
%         case 'No'
%             anotherpass = false;
%         case 'Cancel'
%             error('User cancelled the program');
%     end
    
%end

%% Ask for user input
% 
% dlg_title = 'Please input the optimal parameters';
% prompt = {  'Number of controlpoints n', ...
%     'kth Neighbor for Gaussian weights sigmas', ...
%     'Global scale factors for Gaussian weights', ...
%     'How many CPU cores? (Leave one core unused, if you want to still use this machine for other stuff'...
%       };
% def = {     '40', ...
%     '7', ...
%     '0.9', ...
%     '3', ...
%     };
% num_lines = 1;
% inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
% 
% nControlPoints = str2double(inputdialog{1});
% kthNeighbor = str2double(inputdialog{2});
% globalScale = str2double(inputdialog{3});
% nCores = str2double(inputdialog{4});

%% Calculate the final transform
% % load('EvaluateTransform_3D_workspace.mat');
% 
% matlabpool('open', nCores)
% startTime = tic;
% 
% % Beads in PVA from 20120814
% [tform, matched_cp_transmitted_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
%     matched_cp_reflected(:,5:7),matched_cp_transmitted(:,5:7),...
%     'lwquadratic',nControlPoints,'Gaussian',kthNeighbor,globalScale,true);
% 
% % Transform from 20120703
% % [tform, matched_cp_transmitted_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
% %     matched_cp_reflected(:,5:7),matched_cp_transmitted(:,5:7),...
% %     'lwquadratic',40,'Gaussian',8,0.6,true);
% 
% % [tform, matched_cp_transmitted_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
% %     matched_cp_reflected(:,5:7),matched_cp_transmitted(:,5:7),...
% %     'lwquadratic',60,'Gaussian',7,0.8,true);
% 
% % [tform, matched_cp_transmitted_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
% %     matched_cp_reflected(:,5:7),matched_cp_transmitted(:,5:7),...
% %     'lwlinear',36,'Gaussian',8,0.8,true);
% 
% 
% % Complete Set from Initial Calibration
% % [tform, matched_cp_transmitted_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
% %     matched_cp_reflected(:,5:7),matched_cp_transmitted(:,5:7),...
% %     'lwquadratic',70,'Gaussian',7,1.2,true);
% 
% % [tform, matched_cp_transmitted_trans, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
% %     matched_cp_reflected(:,5:7),matched_cp_transmitted(:,5:7),...
% %     'lwlinear',27,'Gaussian',14,0.6,true);
% 
% 
% elapsedTime = toc(startTime)
% matlabpool close
% 
% deviation = matched_cp_transmitted_trans - matched_cp_reflected(:,5:7);
% mean_deviation = mean(deviation,1)
% std_deviation = std(deviation,1)
% 
% saveas(gcf,['Transform_' tform.method '_RegistrationError.fig']);
% saveas(gcf,['Transform_' tform.method '_RegistrationError.png']);
% 
% save('EvaluateTransform_3D_workspace.mat');
% save(['3D_Transform_' tform.method '.mat'], 'tform', 'matched_cp_reflected', 'matched_cp_transmitted',...
%     'matched_cp_transmitted_trans','FRE', 'TRE', 'FRE_full', 'TRE_full', 'nCores',...
%     'nControlPoints', 'kthNeighbor', 'globalScale',...
%     'matched_cp_reflected', 'matched_cp_transmitted', 'matched_cp_transmitted_trans');

%% Generate Figures
% load('EvaluateTransform_3D_workspace.mat');
% 
% figure
% subplot(1,2,2)
% scatter3(matched_cp_reflected(:,5),matched_cp_reflected(:,6),matched_cp_reflected(:,7),...
%     20,TRE_full(:,1),'filled')
% title('TRE Spatial Distribution');
% colorbar;
% 
% subplot(1,2,1)
% scatter3(matched_cp_reflected(:,5),matched_cp_reflected(:,6),matched_cp_reflected(:,7),...
%     20,FRE_full(:,1),'filled')
% title('FRE Spatial Distribution');
% colorbar;
% 
% saveas(gcf,['Transform_' tform.method '_DeviationSpatialDistribution.fig']);
% saveas(gcf,['Transform_' tform.method '_DeviationSpatialDistribution.png']);
% 
% % figure
% % markerColors = jet(round(max(TRE_full))-floor(min(TRE_full))+1);
% % hold on
% % for a = 1:size(matched_cp_reflected,1)
% %     scatter3(matched_cp_reflected(a,5),matched_cp_reflected(a,6),matched_cp_reflected(a,7),...
% %         40,'filled',...
% %         'MarkerFaceColor', markerColors(round(TRE_full(a))-floor(min(TRE_full))+1,:),...
% %         'MarkerEdgeColor', markerColors(round(TRE_full(a))-floor(min(TRE_full))+1,:));
% % end
% % hold off
% 
% figure
% subplot(1,3,1)
% hist(deviation(:,1),40)
% xlim([-40 40]);
% title('X Deviation');
% xlabel('Distance (nm)');
% ylabel('Frequency');
% legend(['StdDev = ' num2str(std_deviation(1), 3) ' nm']);
% 
% subplot(1,3,2)
% hist(deviation(:,2),40)
% xlim([-40 40]);
% title('Y Deviation');
% xlabel('Distance (nm)');
% ylabel('Frequency');
% legend(['StdDev = ' num2str(std_deviation(2), 3) ' nm']);
% 
% subplot(1,3,3)
% hist(deviation(:,3),40)
% xlim([-40 40]);
% title('Z Deviation');
% xlabel('Distance (nm)');
% ylabel('Frequency');
% legend(['StdDev = ' num2str(std_deviation(3), 3) ' nm']);
% 
% saveas(gcf,['Transform_' tform.method '_XYZDeviations.fig']);
% saveas(gcf,['Transform_' tform.method '_XYZDeviations.png']);

end


