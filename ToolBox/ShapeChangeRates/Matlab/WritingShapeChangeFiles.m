%% Creating shape change rate input for Lola's project
noFiles=5; % number of timepoints (e.g. if the user wants to input 5 shape change maps, noTimePoints=5)
noColumns=55;    % number of columns in the expansion grid file.
midColumn=ceil(noColumns/2);
noRows=4;  % number of rows in the expansion grid file

% Defining the grids where data will be stored temporarily
shapeChange_mag_XY=zeros(noRows,noColumns);   % shape change magnitude in xy
shapeChange_mag_Z=zeros(noRows,noColumns);   % shape change magntude in z


% Filling the shape change grids
for fileID=1:noFiles
    
    shapeChange_mag_XY(:,midColumn-1:midColumn+1)=(fileID-1)*0.25;
    shapeChange_mag_Z(:,midColumn-1:midColumn+1)=(fileID-1)*0.25;

    filename=sprintf('ShapeChangeRate96hrRectangleWingXY_FineGrid_Reduction_%d',fileID-1);
    filename = fopen(filename,'w');
    fprintf(filename,'%d %d \n',[noColumns,noRows]);  % number of columns and rows in the grid matrix
    % writing shapeChange magnitude in x
    for i=1:noRows
        fprintf(filename,'%f ',shapeChange_mag_XY(i,:));
        fprintf(filename,'\n');
    end
    fclose(filename)
    
    filename=sprintf('ShapeChangeRate96hrRectangleWingZ_FineGrid_Reduction_%d',fileID-1);
    filename = fopen(filename,'w');
    fprintf(filename,'%d %d \n',[noColumns,noRows]);  % number of columns and rows in the grid matrix
    % writing shapeChange magnitude in z
    for i=1:noRows
        fprintf(filename,'%f ',shapeChange_mag_Z(i,:));
        fprintf(filename,'\n');
    end
    fclose(filename)  
end