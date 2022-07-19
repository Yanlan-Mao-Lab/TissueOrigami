clear;
clc;
close all;
s = 5; %separation between points, absolute distance (s/R) is the actual angle of separation
R = 24.5; %the outer radius of the intestinol ring shape organoid


Node_z = fix(2* pi * R / s); %integer number of nodes in Z-direction

theta2 = linspace(0,pi/2,(Node_z/4)); 
LOOP = numel(theta2);


MATRIX_x = 0;
MATRIX_y = 0;
MATRIX_z = R;
CHECK = zeros(LOOP-1,1);


  for  i = 1:1:(LOOP-1)
        h = sin(theta2(i)) * R;
        R_local = sqrt(R^2 - h^2);
        C_local = pi*2*R_local;
        node_local = fix(C_local/s);
        disp('there are %d nodes in loop %d,',node_local,i);
        CHECK(i) = node_local;
        theta = linspace(0,2*pi,node_local);

        theta(end) = []; %remove overlapping 360 degree and 0 degree
        Value = numel(theta);
        I = (R_local * 2 * pi / node_local) *  (0.03 * rand(1,Value) -0.015);
        %addition of 3% of noise, -1.5% to 1.5%
        %I = s *  (0.1 * rand(1,Value) -0.005);

        theta_n = theta + I;

        x = zeros(1,Value);
        y = zeros(1,Value);

        for j = 1:1:Value
            x(j) = R_local*cos(theta_n(j));
            y(j) = R_local*sin(theta_n(j));
        end
        z = sqrt(abs(R^2-x.^2-y.^2));
        MATRIX_x = [MATRIX_x,x];
        MATRIX_y = [MATRIX_y,y];
        MATRIX_z = [MATRIX_z,z];

  end
%   MATRIX_x = [MATRIX_x, MATRIX_x];
%   MATRIX_y = [MATRIX_y, MATRIX_y];
%   MATRIX_z = [MATRIX_z, -1*MATRIX_z];  

  MATRIX_x = [MATRIX_x, 0, MATRIX_x((CHECK(1)+1):end)];  %remove repeated nodes (first layer where R_local = R
  MATRIX_y = [MATRIX_y, 0, MATRIX_y((CHECK(1)+1):end)];
  MATRIX_z = [MATRIX_z, -1 * R ,-1*MATRIX_z((CHECK(1)+1):end)];
  Original = [MATRIX_x; MATRIX_y; MATRIX_z];
  Tri = transpose(Original);
  L = 1:1:length(MATRIX_z);
% 
 fig1 = figure;
 scatter3(MATRIX_x,MATRIX_y,MATRIX_z)
%  

DT = delaunayTriangulation(transpose(MATRIX_x), transpose(MATRIX_y), transpose(MATRIX_z));
[triHull, vbOutside, vbInside] = AlphaHull(Tri,s);
 fig2 = figure;
plot3(Tri(:,1),Tri(:,2),Tri(:,3),'ro');
hold on;
fprintf(' * Found %d points for the outer surface\n',sum(vbOutside));
if (sum(vbOutside) > 0)
    trisurf(triHull(vbOutside,:),Tri(:,1),Tri(:,2),Tri(:,3),...
        'FaceColor','cyan','FaceAlpha',0.5)
end


% 
 No_of_nodes = length(MATRIX_x);
 Coordinates = zeros(5,length(MATRIX_x));  %%% Might  be able to use, as the connectivity is set in the .m function under special order
 %%% use command unique()  (?)
 No_of_Connectivity = length(triHull); 
 Connectivity = zeros(3,length(triHull));

 %prepare the required input format for the SphericalTriangulation.text file 


 for k = 1:1:(length(MATRIX_x))
     Coordinates(1,k) = k-1;  % C++ starts with 0 
     Coordinates(2,k) = Tri(k,1); %transpose reqiured
     Coordinates(3,k) = Tri(k,2);
     Coordinates(4,k) = Tri(k,3);
     Coordinates(5,k) = 0;
 end

% all index starts with zero in C++
 for m = 1:1:(length(triHull))
     Connectivity(1,m) = triHull(m,1)-1; %transpose reqiured
     Connectivity(2,m) = triHull(m,2)-1;
     Connectivity(3,m) = triHull(m,3)-1; %m columns 3 rows
 end
% 
% 
% 



% 
filespec = fullfile('/Users/jefferywei/Desktop/IPLs_Mesh/', 'MeshPointResult.txt' ); %full file path
fileID = fopen(filespec,'w');
fprintf(fileID,'%d\n',No_of_nodes);
fprintf(fileID, '%d %.4f %.4f %.4f %d\n',Coordinates);
fprintf(fileID, '%d\n',No_of_Connectivity);
fprintf(fileID, '%d %d %d\n',Connectivity);
