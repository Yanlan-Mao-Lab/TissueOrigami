%%
% This script is to create stiffness grid files for each time step. This
% is of course for instances where limited number of input files are used.
% If a better time resolution is needed, I have to modify Melda's code in
% order to be able to read spatio-temporal growth/expansion maps.
clear all
close all
%% This part should be defined by user
noTimePoints=3; % number of timepoints (e.g. if the user wants to input 5 growth maps, noTimePoints=5)
noPeriods=3;   % number of periods (gyrus or sulcus) that should appear in the tissue
noColumns=4*noPeriods+1;    % number of columns in the expansion grid file. For each square on the expansion rate grid, we need three values (i.e. expansion rate in x, y, z).
noRows=10;  % number of rows in the expansion grid file

% define stiffnesses for control and folding, for 8hr and 16hr
stiffnessAtG_ctrl_0hr = 37.49;
stiffnessAtS_ctrl_0hr = 37.49;
stiffnessAtSG_ctrl_0hr = (stiffnessAtG_ctrl_0hr+stiffnessAtS_ctrl_0hr)/2;
stiffnessAtGS_ctrl_0hr = (stiffnessAtG_ctrl_0hr+stiffnessAtS_ctrl_0hr)/2;

stiffnessAtG_ctrl_8hr = 37.49;
stiffnessAtS_ctrl_8hr = 37.49;
stiffnessAtSG_ctrl_8hr = (stiffnessAtG_ctrl_8hr+stiffnessAtS_ctrl_8hr)/2;
stiffnessAtGS_ctrl_8hr = (stiffnessAtG_ctrl_8hr+stiffnessAtS_ctrl_8hr)/2;

stiffnessAtG_ctrl_16hr = 82.27;
stiffnessAtS_ctrl_16hr = 82.27;
stiffnessAtSG_ctrl_16hr = (stiffnessAtG_ctrl_16hr+stiffnessAtS_ctrl_16hr)/2;
stiffnessAtGS_ctrl_16hr = (stiffnessAtG_ctrl_16hr+stiffnessAtS_ctrl_16hr)/2;

stiffnessAtG_folding_0hr = 37.49;
stiffnessAtS_folding_0hr = 37.49;
stiffnessAtSG_folding_0hr = (stiffnessAtG_folding_0hr+stiffnessAtS_folding_0hr)/2;
stiffnessAtGS_folding_0hr = (stiffnessAtG_folding_0hr+stiffnessAtS_folding_0hr)/2;

stiffnessAtG_folding_8hr = 458.83;
stiffnessAtS_folding_8hr = 71.71;
stiffnessAtSG_folding_8hr = (stiffnessAtG_folding_8hr+stiffnessAtS_folding_8hr)/2;
stiffnessAtGS_folding_8hr = (stiffnessAtG_folding_8hr+stiffnessAtS_folding_8hr)/2;

stiffnessAtG_folding_16hr = 793.77;
stiffnessAtS_folding_16hr = 59.43;
stiffnessAtSG_folding_16hr = (stiffnessAtG_folding_16hr+stiffnessAtS_folding_16hr)/2;
stiffnessAtGS_folding_16hr = (stiffnessAtG_folding_16hr+stiffnessAtS_folding_16hr)/2;
%% Creating the grids and writing them into file

% For expansion rates, we need to make two matrices, one for magnitude and
% one for orientation. Since for the case of brain folding, we don't have y
% growth, we either have orientations as 0 or 180 degrees of x or z axis.
% So I will first create the magnitude matrix. Then find the ones with
% negative values, set the orientation of those to 180. And then make all
% values in the magnitude matrix positive.
for timePointsID=1:noTimePoints
    stiffness_ctrl=zeros(noRows,noColumns);   % stiffness grid
    stiffness_folding=zeros(noRows,noColumns);   % stiffness grid
    for j=1:noRows
        for i=1:noPeriods
            if timePointsID==1
                stiffness_ctrl(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_ctrl_0hr,stiffnessAtG_ctrl_0hr,stiffnessAtGS_ctrl_0hr,stiffnessAtS_ctrl_0hr];
                stiffness_folding(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_folding_0hr,stiffnessAtG_folding_0hr,stiffnessAtGS_folding_0hr,stiffnessAtS_folding_0hr];
            elseif timePointsID==2
                stiffness_ctrl(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_ctrl_8hr,stiffnessAtG_ctrl_8hr,stiffnessAtGS_ctrl_8hr,stiffnessAtS_ctrl_8hr];
                stiffness_folding(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_folding_8hr,stiffnessAtG_folding_8hr,stiffnessAtGS_folding_8hr,stiffnessAtS_folding_8hr];
            elseif timePointsID==3
                stiffness_ctrl(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_ctrl_16hr,stiffnessAtG_ctrl_16hr,stiffnessAtGS_ctrl_16hr,stiffnessAtS_ctrl_16hr];
                stiffness_folding(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_folding_16hr,stiffnessAtG_folding_16hr,stiffnessAtGS_folding_16hr,stiffnessAtS_folding_16hr];
            end
        end
        if timePointsID==1
            stiffness_ctrl(j,noColumns)=stiffnessAtSG_ctrl_0hr;
            stiffness_folding(j,noColumns)=stiffnessAtSG_folding_0hr;
        elseif timePointsID==2
            stiffness_ctrl(j,noColumns)=stiffnessAtSG_ctrl_8hr;
            stiffness_folding(j,noColumns)=stiffnessAtSG_folding_8hr;
        elseif timePointsID==3
            stiffness_ctrl(j,noColumns)=stiffnessAtSG_ctrl_16hr;
            stiffness_folding(j,noColumns)=stiffnessAtSG_folding_16hr;
        end
    end
    
    
    
    filename=sprintf('StiffnessCtrl_3periods_10rows_t%d',timePointsID-1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%d %d \n',[noColumns,noRows]);  % number of columns and rows in the grid matrix
    % writing stiffness file for ctrl case
    for i=1:noRows
        fprintf(fileID,'%f ',stiffness_ctrl(i,:));
        fprintf(fileID,'\n');
    end
    fclose(fileID)
    
    filename=sprintf('StiffnessFolding_3periods_10rows_t%d',timePointsID-1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%d %d \n',[noColumns,noRows]);  % number of columns and rows in the grid matrix
    % writing stiffness file for folding case
    for i=1:noRows
        fprintf(fileID,'%f ',stiffness_folding(i,:));
        fprintf(fileID,'\n');
    end
    fclose(fileID)
end
%% This is just for plotting the stiffnesses
width = 6;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 2;      % LineWidth
msz = 7; % MarkerSize

fig1=figure()
pos=get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) width*50 height*50]);
set(gca,'FontSize',fsz,'LineWidth',lw)

cmap1=gray(10);
cmap2=autumn(10);

%plot([stiffnessAtSG_ctrl_0hr,stiffnessAtG_ctrl_0hr,stiffnessAtGS_ctrl_0hr,stiffnessAtS_ctrl_0hr,stiffnessAtSG_ctrl_0hr],'*-','color',cmap1(3,:),'Linewidth',lw)
%hold on
%plot([stiffnessAtSG_ctrl_8hr,stiffnessAtG_ctrl_8hr,stiffnessAtGS_ctrl_8hr,stiffnessAtS_ctrl_8hr,stiffnessAtSG_ctrl_8hr],'*-','color',cmap1(5,:),'Linewidth',lw)
%hold on
%plot([stiffnessAtSG_ctrl_16hr,stiffnessAtG_ctrl_16hr,stiffnessAtGS_ctrl_16hr,stiffnessAtS_ctrl_16hr,stiffnessAtSG_ctrl_16hr],'*-','color',cmap1(7,:),'Linewidth',lw)
%hold on
plot([stiffnessAtSG_folding_0hr,stiffnessAtG_folding_0hr,stiffnessAtGS_folding_0hr,stiffnessAtS_folding_0hr,stiffnessAtSG_folding_0hr],'o-','color',cmap1(3,:),'Linewidth',lw)
hold on
plot([stiffnessAtSG_folding_8hr,stiffnessAtG_folding_8hr,stiffnessAtGS_folding_8hr,stiffnessAtS_folding_8hr,stiffnessAtSG_folding_8hr],'o-','color',cmap1(5,:),'Linewidth',lw)
hold on
plot([stiffnessAtSG_folding_16hr,stiffnessAtG_folding_16hr,stiffnessAtGS_folding_16hr,stiffnessAtS_folding_16hr,stiffnessAtSG_folding_16hr],'o-','color',cmap1(7,:),'Linewidth',lw)

set(gcf,'Position',[pos(1) pos(2) width*60 height*60]);
set(gca,'FontSize',fsz,'LineWidth',lw)
%xlim([-5 100])
ylim([0 1000])
xticklabels({'S->G','G','G->S','S','S->G'})
legend('ctrl 0 hr','ctrl 8 hr','ctrl 16 hr','folding 0 hr','folding 8 hr','folding 16 hr','location','NorthEast')
xlabel('Position along X axis');ylabel('Stiffness (Pa)');
%% Extrapolating for 0 hr
width = 6;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 14;      % Fontsize
lw = 2;      % LineWidth
msz = 7; % MarkerSize

fig1=figure()
pos=get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) width*50 height*50]);
set(gca,'FontSize',fsz,'LineWidth',lw)

x=[8,16];
y=[37.49,82.27];

plot(x,y,'bo','color',cmap1(2,:),'Linewidth',lw)
set(gcf,'Position',[pos(1) pos(2) width*60 height*60]);
set(gca,'FontSize',fsz,'LineWidth',lw)
xlim([0 16])
ylim([0 100])
xlabel('time (hr)');ylabel('Stiffness (Pa)');
%% Creating the grids and writing them into file for ventricular surface

% For expansion rates, we need to make two matrices, one for magnitude and
% one for orientation. Since for the case of brain folding, we don't have y
% growth, we either have orientations as 0 or 180 degrees of x or z axis.
% So I will first create the magnitude matrix. Then find the ones with
% negative values, set the orientation of those to 180. And then make all
% values in the magnitude matrix positive.
for timePointsID=1:noTimePoints
    stiffness_folding=zeros(noRows,noColumns);   % stiffness grid
    for j=1:noRows
        for i=1:noPeriods
            if timePointsID==1
                stiffness_folding(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_folding_0hr,stiffnessAtSG_folding_0hr,stiffnessAtGS_folding_0hr,stiffnessAtGS_folding_0hr];
            elseif timePointsID==2
                stiffness_folding(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_folding_8hr,stiffnessAtSG_folding_8hr,stiffnessAtGS_folding_8hr,stiffnessAtGS_folding_8hr];
            elseif timePointsID==3
                stiffness_folding(j,(i-1)*4+1:(i-1)*4+4)=[stiffnessAtSG_folding_16hr,stiffnessAtSG_folding_16hr,stiffnessAtGS_folding_16hr,stiffnessAtGS_folding_16hr];
            end
        end
        if timePointsID==1
            stiffness_folding(j,noColumns)=stiffnessAtSG_folding_0hr;
        elseif timePointsID==2
            stiffness_folding(j,noColumns)=stiffnessAtSG_folding_8hr;
        elseif timePointsID==3
            stiffness_folding(j,noColumns)=stiffnessAtSG_folding_16hr;
        end
    end
    
    
    filename=sprintf('StiffnessFolding_3periods_10rows_VZ_t%d',timePointsID-1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%d %d \n',[noColumns,noRows]);  % number of columns and rows in the grid matrix
    % writing stiffness file for folding case
    for i=1:noRows
        fprintf(fileID,'%f ',stiffness_folding(i,:));
        fprintf(fileID,'\n');
    end
    fclose(fileID)
end
%% Creating stiffness grids for Lola's project


noFiles=5; %
noColumns=11;    %
noRows=4;  %

initialStiffness=1600;


for fileID=1:noFiles
    stiffness=zeros(noRows,noColumns)+1600;   % stiffness grid
    
    stiffness(:,5:7)=(1-(fileID-1)*0.25)*stiffness(:,5:7);
    if fileID==noFiles
        stiffness(:,5:7)=1;
    end
 
    
    filename=sprintf('Stiffness96hrRectangleWing_Reduction_%d',fileID-1);
    filename = fopen(filename,'w');
    fprintf(filename,'%d %d \n',[noColumns,noRows]);  % number of columns and rows in the grid matrix
    % writing stiffness file for folding case
    for i=1:noRows
        fprintf(filename,'%f ',stiffness(i,:));
        fprintf(filename,'\n');
    end
    fclose(filename)
end