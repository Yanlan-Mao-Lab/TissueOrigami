% %%
% % This script is to create growth/expansion grid files for each time step. This
% % is of course for instances where limited number of input files are used.
% % If a better time resolution is needed, I have to modify Melda's code in
% % order to be able to read spatio-temporal growth/expansion maps.
% clear all
% close all
% %% This part should be defined by user
% noTimePoints=5; % number of timepoints (e.g. if the user wants to input 5 growth maps, noTimePoints=5)
% timePoints=[0.14,2.14,3.14,5.14,10.14]; % the timepoints that the user whats from the input excel sheet
% noPeriods=3;   % number of periods (gyrus or sulcus) that should appear in the tissue
% noColumns=(4*noPeriods+1)*3;    % number of columns in the expansion grid file. For each square on the expansion rate grid, we need three values (i.e. expansion rate in x, y, z).
% noRows=10;  % number of rows in the expansion grid file
% umToPixelRatio = 10 % I have rescaled the lengths such that 10 um is 1 pixel
% timeStepToSecRatio = 0.25    % Each time step is 15 min. The input files should be pixel per hour
% 
% %% Converting the data into pixel/sec
% ImportExpansionRateRawData  % call the file that imports data
% cd
% AvgGvelocityX =  AvgGvelocityX/(umToPixelRatio*timeStepToSecRatio);
% AvgGvelocityZ = AvgGvelocityZ/(umToPixelRatio*timeStepToSecRatio);
% AvgSGvelocityX = AvgSGvelocityX/(umToPixelRatio*timeStepToSecRatio);
% AvgSGvelocityZ = AvgSGvelocityZ/(umToPixelRatio*timeStepToSecRatio);
% AvgSvelocityX = AvgSvelocityX/(umToPixelRatio*timeStepToSecRatio);
% AvgSvelocityZ = AvgSvelocityZ/(umToPixelRatio*timeStepToSecRatio);
% AvgGSvelocityX = AvgGSvelocityX/(umToPixelRatio*timeStepToSecRatio);
% AvgGSvelocityZ = AvgGSvelocityZ/(umToPixelRatio*timeStepToSecRatio);
% 
% %% Creating the grids and writing them into file
% % Defining the grids where data will be stored temporarily
% growth_mag_X=zeros(noRows,noColumns);   % growth magnitude in x
% growth_mag_Z=zeros(noRows,noColumns);   % growth magntude in z
% growth_mag_X_temp=zeros(noRows,noColumns/3);   % growth magnitude in x
% growth_mag_Z_temp=zeros(noRows,noColumns/3);   % growth magntude in z
% growth_ori_X=zeros(noRows,noColumns/3);   % growth orientation in x, noColumns is 1/3 of the growth magnitude because we don't need x,y,z
% growth_ori_Z=zeros(noRows,noColumns/3);   % growth orientation in z, noColumns is 1/3 of the growth magnitude because we don't need x,y,z
% 
% 
% % For expansion rates, we need to make two matrices, one for magnitude and
% % one for orientation. Since for the case of brain folding, we don't have y
% % growth, we either have orientations as 0 or 180 degrees of x or z axis.
% % So I will first create the magnitude matrix. Then find the ones with
% % negative values, set the orientation of those to 180. And then make all
% % values in the magnitude matrix positive.
% for timePointsID=1:noTimePoints
%     rowID=find(Timeframesec==timePoints(timePointsID));
%     for j=1:noRows
%         for i=1:noPeriods
%             growth_mag_X(j,(i-1)*12+1:(i-1)*12+12)=[AvgSGvelocityX(rowID),0,0,AvgGvelocityX(rowID),0,0,AvgGSvelocityX(rowID),0,0,AvgSvelocityX(rowID),0,0];
%             growth_mag_X_temp(j,(i-1)*4+1:(i-1)*4+4)=[AvgSGvelocityX(rowID),AvgGvelocityX(rowID),AvgGSvelocityX(rowID),AvgSvelocityX(rowID)];
%             growth_mag_Z(j,(i-1)*12+1:(i-1)*12+12)=[0,0,AvgSGvelocityZ(rowID),0,0,AvgGvelocityZ(rowID),0,0,AvgGSvelocityZ(rowID),0,0,AvgSvelocityZ(rowID)];
%             growth_mag_Z_temp(j,(i-1)*4+1:(i-1)*4+4)=[AvgSGvelocityZ(rowID),AvgGvelocityZ(rowID),AvgGSvelocityZ(rowID),AvgSvelocityZ(rowID)];
%         end
%         growth_mag_X(j,noColumns-2:noColumns)=[AvgSGvelocityX(rowID),0,0];
%         growth_mag_X_temp(j,noColumns/3)=AvgSGvelocityX(rowID);
%         growth_mag_Z(j,noColumns-2:noColumns)=[0,0,AvgSGvelocityZ(rowID)];
%         growth_mag_Z_temp(j,noColumns/3)=AvgSGvelocityZ(rowID);
%     end
%     
%     negID_X=zeros(noRows,noColumns/3);
%     negID_X(find(growth_mag_X_temp<0))=1;    % find the positions in magnitude matrix with negative values
%     growth_ori_X=180*negID_X;   % set the orientation of those positions to 180
%     growth_mag_X=abs(growth_mag_X); % make all values in magnitude matrix positive
%     
%     negID_Z=zeros(noRows,noColumns/3);
%     negID_Z(find(growth_mag_Z_temp<0))=1;    % find the positions in magnitude matrix with negative values
%     growth_ori_Z=180*negID_Z;   % set the orientation of those positions to 180
%     growth_mag_Z=abs(growth_mag_Z); % make all values in magnitude matrix positive
%     
%     filename=sprintf('ExpansionRateX_3periods_10rows_t%d',timePointsID-1);
%     fileID = fopen(filename,'w');
%     fprintf(fileID,'%d %d \n',[noColumns/3,noRows]);  % number of columns and rows in the grid matrix
%     % writing growth magnitude in x
%     for i=1:noRows
%         fprintf(fileID,'%f ',growth_mag_X(i,:));
%         fprintf(fileID,'\n');
%     end
%     % writing growth orientation in x
%     for i=1:noRows
%         fprintf(fileID,'%d ',growth_ori_X(i,:));
%         fprintf(fileID,'\n');
%     end
%     fclose(fileID)
%     
%     filename=sprintf('ExpansionRateZ_3periods_10rows_t%d',timePointsID-1);
%     fileID = fopen(filename,'w');
%     fprintf(fileID,'%d %d \n',[noColumns/3,noRows]);  % number of columns and rows in the grid matrix
%     % writing growth magnitude in z
%     for i=1:noRows
%         fprintf(fileID,'%f ',growth_mag_Z(i,:));
%         fprintf(fileID,'\n');
%     end
%     % writing growth orientation in z
%     for i=1:noRows
%         fprintf(fileID,'%d ',growth_ori_Z(i,:));
%         fprintf(fileID,'\n');
%     end
%     fclose(fileID)
%     
%     filename=sprintf('ExpansionRateZOpposite_3periods_10rows_t%d',timePointsID-1);
%     fileID = fopen(filename,'w');
%     fprintf(fileID,'%d %d \n',[noColumns/3,noRows]);  % number of columns and rows in the grid matrix
%     % writing growth magnitude in z
%     for i=1:noRows
%         fprintf(fileID,'%f ',growth_mag_Z(i,:));
%         fprintf(fileID,'\n');
%     end
%     % writing growth orientation in z
%     for i=1:noRows
%         fprintf(fileID,'%d ',growth_ori_Z(i,:)+180);
%         fprintf(fileID,'\n');
%     end
%     fclose(fileID)
%     
% end
%% Creating growth/shape change rate input for Lola's project
noFiles=5; % number of timepoints (e.g. if the user wants to input 5 growth maps, noTimePoints=5)
noColumns=33;    % number of columns in the expansion grid file. For each square on the expansion rate grid, we need three values (i.e. expansion rate in x, y, z).
noRows=4;  % number of rows in the expansion grid file

% Defining the grids where data will be stored temporarily
growth_mag_XY=zeros(noRows,noColumns);   % growth magnitude in xy
growth_mag_Z=zeros(noRows,noColumns);   % growth magntude in z
growth_ori_XY=zeros(noRows,noColumns/3);   % growth orientation in xy, noColumns is 1/3 of the growth magnitude because we don't need x,y,z
growth_ori_Z=zeros(noRows,noColumns/3);   % growth orientation in z, noColumns is 1/3 of the growth magnitude because we don't need x,y,z


% For expansion rates, we need to make two matrices, one for magnitude and
% one for orientation. Since for the case of brain folding, we don't have y
% growth, we either have orientations as 0 or 180 degrees of x or z axis.
% So I will first create the magnitude matrix. Then find the ones with
% negative values, set the orientation of those to 180. And then make all
% values in the magnitude matrix positive.
for fileID=1:noFiles
    for i=5:7
        growth_mag_XY(:,(i-1)*3+1)= -(fileID-1)*0.25;
        growth_mag_XY(:,(i-1)*3+2)= -(fileID-1)*0.25;
        growth_mag_XY(:,(i-1)*3+3)= 0;
        
        growth_mag_Z(:,(i-1)*3+1)= 0;
        growth_mag_Z(:,(i-1)*3+2)= 0;
        growth_mag_Z(:,(i-1)*3+3)= -(fileID-1)*0.25;
    end
        
    growth_ori_XY(:,5:7)=45;
    
    filename=sprintf('ShapeChangeRate96hrRectangleWingXY_Reduction_%d',fileID-1);
    filename = fopen(filename,'w');
    fprintf(filename,'%d %d \n',[noColumns/3,noRows]);  % number of columns and rows in the grid matrix
    % writing growth magnitude in x
    for i=1:noRows
        fprintf(filename,'%f ',growth_mag_XY(i,:));
        fprintf(filename,'\n');
    end
    % writing growth orientation in x
    for i=1:noRows
        fprintf(filename,'%d ',growth_ori_XY(i,:));
        fprintf(filename,'\n');
    end
    fclose(filename)
    
    filename=sprintf('ShapeChangeRate96hrRectangleWingZ_Reduction_%d',fileID-1);
    filename = fopen(filename,'w');
    fprintf(filename,'%d %d \n',[noColumns/3,noRows]);  % number of columns and rows in the grid matrix
    % writing growth magnitude in z
    for i=1:noRows
        fprintf(filename,'%f ',growth_mag_Z(i,:));
        fprintf(filename,'\n');
    end
    % writing growth orientation in z
    for i=1:noRows
        fprintf(filename,'%d ',growth_ori_Z(i,:));
        fprintf(filename,'\n');
    end
    fclose(filename)
    
end
