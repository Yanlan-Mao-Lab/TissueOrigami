clear;
clc;
% 
no_Columns = 50;
no_Rows = 50;
R = 15;
Stiffness = zeros(no_Columns , no_Rows);

centre_x = no_Columns/2;
centre_y = no_Rows/2;
x_start = centre_x - R;
x_end = centre_x + R;
y_start = centre_y - R;
y_end = centre_y + R;

% circular region of stiffness softening
 for i = 1:1:no_Rows
     if (x_start <= i) && (i <= x_end)
         for j = 1:1:no_Columns
             if sqrt((centre_x - i)^2 + (centre_y - j)^2) > R
                Stiffness(i,j) = 1600;
             else 
                Stiffness(i,j) = 800;
             end
         end
     else
        for j = 1:1:no_Columns
            Stiffness(i,j) = 1600;
        end
     end
 end


                




% S = (no_Rows - D) / 2;
% Start = fix(S) + 1;
% disp(Start);
% End = Start + D -1;
% disp(End);










%  sqaure region stiffness softening
% for i = 1:1:no_Rows
%     if (Start <= i) && (i <= End)
%         for j = 1:1:no_Columns
%             if (Start <= j) && (j <= End)
%                 Stiffness(i,j) = 800;
%             else 
%                 Stiffness(i,j) = 1600;
%             end
%         end
%     else
%         for j = 1:1:no_Columns
%             Stiffness(i,j) = 1600;
%         end
%     end
% end

          

filespec = fullfile('/Users/jefferywei/Desktop/IPLs_Mesh/', 'StiffnessMatrix.txt' ); %full file path
fileID = fopen(filespec,'w');
fprintf(fileID,'%d ',no_Columns);
fprintf(fileID,'%d \n',no_Rows);
for i = 1:1:no_Rows
    for j = 1:1:no_Columns
        if j == no_Columns
            fprintf(fileID,'%d\n',Stiffness(i,j));
        else
            fprintf(fileID, '%d ',Stiffness(i,j));
        end
    end
end