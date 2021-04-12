%
% %%
% % 20200501 - Testing mesh that works
% clar all
% close all
% k=1;
% xLength=100;
% yLength=30;
% desiredSideLength=5;
%
% for i=0:desiredSideLength:xLength
%     for j=0:desiredSideLength:yLength
%         x=i+1*(rand);
%         y=j+1*(rand);
%         if x>0 && x<1
%             x=0.0;
%         elseif x>xLength && x<xLength+1
%             x=xLength;
%         end
%         if y>0 && y<1
%             y=0.0;
%         elseif y>yLength && y<yLength+1
%             y=yLength;
%         end
%         coord(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
% %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         NodeInfo(i,3)=1;
%     end
% end
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('Points.1.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('Points.1.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);

% %%
% %20200504 - Rough experimental measures - small (1/5th) mesh
% clear all
% close all
% k=1;
% xLength=200;
% yLength=40;
% desiredSideLength=5;
%
% for i=0:desiredSideLength:xLength
%     for j=0:desiredSideLength:yLength
%         if i==0 || i==xLength || j==0 || j==yLength
%             x=i;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% scatter(borderNodes(:,1),borderNodes(:,2),'r')
%
% k=1;
% for i=desiredSideLength:desiredSideLength:xLength-desiredSideLength
%     for j=desiredSideLength:desiredSideLength:yLength-desiredSideLength
%         x=i;
%         y=j;
%         surfaceNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
% hold on
% scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')
%
%
% random=randperm(length(surfaceNodes));
% for i=1:length(random)
%     x_noise=desiredSideLength*(0.5*rand-0.25);
%     y_noise=desiredSideLength*(0.5*rand-0.25);
%     if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
%     else
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
%     end
%     if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
%     else
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
%     end
% end
%
% coord=[borderNodes;surfaceNodes];
%
% figure()
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         NodeInfo(i,3)=1;
%     end
% end
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('Points.1.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('Points.1.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);
% %%
% % 20200504 - Rough experimental measures - large mesh
% clear all
% close all
% k=1;
% xLength=200;
% yLength=40;
% desiredSideLength=50;
%
% for i=0:desiredSideLength:xLength
%     for j=0:desiredSideLength:yLength
%         if i==0 || i==xLength || j==0 || j==yLength
%             x=i;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% scatter(borderNodes(:,1),borderNodes(:,2),'r')
%
% k=1;
% for i=desiredSideLength:desiredSideLength:xLength-desiredSideLength
%     for j=desiredSideLength:desiredSideLength:yLength-desiredSideLength
%         x=i;
%         y=j;
%         surfaceNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
% hold on
% scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')
%
%
% random=randperm(length(surfaceNodes));
% for i=1:length(random)
%     x_noise=desiredSideLength*(0.5*rand-0.25);
%     y_noise=desiredSideLength*(0.5*rand-0.25);
%     if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
%     else
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
%     end
%     if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
%     else
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
%     end
% end
%
% coord=[borderNodes;surfaceNodes];
%
% figure()
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         NodeInfo(i,3)=1;
%     end
% end
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('Points.1.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('Points.1.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);
% %%
% % 20200504 - Rough experimental measures - thin tissue - small mesh
% clear all
% close all
% k=1;
% xLength=240;
% yLength=40;
% desiredSideLength=5;
%
% for i=0:desiredSideLength:xLength
%     for j=0:desiredSideLength:yLength
%         if i==0 || i==xLength || j==0 || j==yLength
%             x=i;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% scatter(borderNodes(:,1),borderNodes(:,2),'r')
%
% k=1;
% for i=desiredSideLength:desiredSideLength:xLength-desiredSideLength
%     for j=desiredSideLength:desiredSideLength:yLength-desiredSideLength
%         x=i;
%         y=j;
%         surfaceNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
% if mod(ceil(xLength),desiredSideLength)~=0
%     for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
%             x=floor(xLength/desiredSideLength)*desiredSideLength;
%             y=j;
%             if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%                 surfaceNodes(k,:)=[x,y];
%                 k=k+1;
%             end
%     end
% end
%
%
% if mod(ceil(yLength),desiredSideLength)~=0
%     for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
%             x=i;
%             y=floor(yLength/desiredSideLength)*desiredSideLength;
%             if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%                 surfaceNodes(k,:)=[x,y];
%                 k=k+1;
%             end
%     end
% end
%
% hold on
% scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')
%
%
% random=randperm(length(surfaceNodes));
% for i=1:length(random)
%     x_noise=desiredSideLength*(0.5*rand-0.25);
%     y_noise=desiredSideLength*(0.5*rand-0.25);
%     if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
%     else
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
%     end
%     if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
%     else
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
%     end
% end
%
% coord=[borderNodes;surfaceNodes];
%
% figure()
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
%         %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         NodeInfo(i,3)=1;
%         hold on
%         plot(NodeInfo(i,1),NodeInfo(i,2),'go');
%     end
% end
%
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('thinTissue_smallMesh.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('thinTissue_smallMesh.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);
% %%
% 20200504 - Rough experimental measures - thin tissue - fine mesh
clear all
close all
k=1;
xLength=240;
yLength=40;
desiredSideLength=2;

for i=0:desiredSideLength:xLength
    for j=0:desiredSideLength:yLength
        if i==0 || i==xLength || j==0 || j==yLength
            x=i;
            y=j;
            borderNodes(k,:)=[x,y];
            k=k+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')

k=1;
for i=desiredSideLength:desiredSideLength:xLength-desiredSideLength
    for j=desiredSideLength:desiredSideLength:yLength-desiredSideLength
        x=i;
        y=j;
        surfaceNodes(k,:)=[x,y];
        k=k+1;
    end
end

if mod(ceil(xLength),desiredSideLength)~=0
    for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
            x=floor(xLength/desiredSideLength)*desiredSideLength;
            y=j;
            if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
                surfaceNodes(k,:)=[x,y];
                k=k+1;
            end
    end
end


if mod(ceil(yLength),desiredSideLength)~=0
    for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
            x=i;
            y=floor(yLength/desiredSideLength)*desiredSideLength;
            if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
                surfaceNodes(k,:)=[x,y];
                k=k+1;
            end
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')


random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);

for i=1:length(NodeInfo)
    %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
        %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end

ConnectionInfo = DT.ConnectivityList;

fileID = fopen('thinTissue_fineMesh.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thinTissue_fineMesh.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
%%
% 20200504 - Rough experimental measures - thick tissue - small mesh
clear all
close all
k=1;
xLength=240;
yLength=82.57;
desiredSideLength=5;

for i=0:desiredSideLength:ceil(xLength)
    for j=0:desiredSideLength:ceil(yLength)
        if i==0 || i==ceil(xLength) || j==0 || j==ceil(yLength)
            x=i;
            y=j;
            borderNodes(k,:)=[x,y];
            k=k+1;
        end
    end
end

if mod(ceil(xLength),desiredSideLength)~=0
    for j=0:desiredSideLength:ceil(yLength)
            x=xLength;
            y=j;
            borderNodes(k,:)=[x,y];
            k=k+1;
    end
end

if mod(ceil(yLength),desiredSideLength)~=0
    for i=0:desiredSideLength:ceil(xLength)
            x=i;
            y=yLength;
            borderNodes(k,:)=[x,y];
            k=k+1;
    end
end
scatter(borderNodes(:,1),borderNodes(:,2),'r')

k=1;
for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
    for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
        x=i;
        y=j;
        surfaceNodes(k,:)=[x,y];
        k=k+1;
    end
end


if mod(ceil(xLength),desiredSideLength)~=0
    for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
            x=floor(xLength/desiredSideLength)*desiredSideLength;
            y=j;
            if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
                surfaceNodes(k,:)=[x,y];
                k=k+1;
            end
    end
end


if mod(ceil(yLength),desiredSideLength)~=0
    for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
            x=i;
            y=floor(yLength/desiredSideLength)*desiredSideLength;
            if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
                surfaceNodes(k,:)=[x,y];
                k=k+1;
            end
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')


random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);

for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end

ConnectionInfo = DT.ConnectivityList;

fileID = fopen('thickTissue_smallMesh.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thickTissue_smallMesh.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
% %%
% % 20200504 - Rough experimental measures - thick tissue - fine mesh
% clear all
% close all
% k=1;
% xLength=240;
% yLength=82.57;
% desiredSideLength=2;
%
% for i=0:desiredSideLength:ceil(xLength)
%     for j=0:desiredSideLength:ceil(yLength)
%         if i==0 || i==ceil(xLength) || j==0 || j==ceil(yLength)
%             x=i;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% if mod(ceil(xLength),desiredSideLength)~=0
%     for j=0:desiredSideLength:ceil(yLength)
%             x=xLength;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%     end
% end
%
% if mod(ceil(yLength),desiredSideLength)~=0
%     for i=0:desiredSideLength:ceil(xLength)
%             x=i;
%             y=yLength;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%     end
% end
% scatter(borderNodes(:,1),borderNodes(:,2),'r')
%
% k=1;
% for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
%     for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
%         x=i;
%         y=j;
%         surfaceNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
%
% if mod(ceil(xLength),desiredSideLength)~=0
%     for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
%             x=floor(xLength/desiredSideLength)*desiredSideLength;
%             y=j;
%             if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%                 surfaceNodes(k,:)=[x,y];
%                 k=k+1;
%             end
%     end
% end
%
%
% if mod(ceil(yLength),desiredSideLength)~=0
%     for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
%             x=i;
%             y=floor(yLength/desiredSideLength)*desiredSideLength;
%             if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%                 surfaceNodes(k,:)=[x,y];
%                 k=k+1;
%             end
%     end
% end
%
% hold on
% scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')
%
%
% random=randperm(length(surfaceNodes));
% for i=1:length(random)
%     x_noise=desiredSideLength*(0.5*rand-0.25);
%     y_noise=desiredSideLength*(0.5*rand-0.25);
%     if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
%     else
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
%     end
%     if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
%     else
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
%     end
% end
%
% coord=[borderNodes;surfaceNodes];
%
% figure()
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
%         NodeInfo(i,3)=1;
%         hold on
%         plot(NodeInfo(i,1),NodeInfo(i,2),'go');
%     end
% end
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('thickTissue_fineMesh.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('thickTissue_fineMesh.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);
% %%
% % 20210215 - Rough experimental measures - thin tissue - fine mesh - 5
% % periods
% clear all
% close all
% k=1;
% xLength=120;
% yLength=40;
% desiredSideLength=2;
%
% for i=0:desiredSideLength:xLength
%     for j=0:desiredSideLength:yLength
%         if i==0 || i==xLength || j==0 || j==yLength
%             x=i;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% scatter(borderNodes(:,1),borderNodes(:,2),'r')
%
% k=1;
% for i=desiredSideLength:desiredSideLength:xLength-desiredSideLength
%     for j=desiredSideLength:desiredSideLength:yLength-desiredSideLength
%         x=i;
%         y=j;
%         surfaceNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
% if mod(ceil(xLength),desiredSideLength)~=0
%     for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
%         x=floor(xLength/desiredSideLength)*desiredSideLength;
%         y=j;
%         if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%             surfaceNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
%
% if mod(ceil(yLength),desiredSideLength)~=0
%     for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
%         x=i;
%         y=floor(yLength/desiredSideLength)*desiredSideLength;
%         if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%             surfaceNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% hold on
% scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')
%
%
% random=randperm(length(surfaceNodes));
% for i=1:length(random)
%     x_noise=desiredSideLength*(0.5*rand-0.25);
%     y_noise=desiredSideLength*(0.5*rand-0.25);
%     if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
%     else
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
%     end
%     if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
%     else
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
%     end
% end
%
% coord=[borderNodes;surfaceNodes];
%
% figure()
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
%         %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
%         NodeInfo(i,3)=1;
%         hold on
%         plot(NodeInfo(i,1),NodeInfo(i,2),'go');
%     end
% end
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('thinTissue_fineMesh_5periods.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('thinTissue_fineMesh_5periods.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);
% %%
% % 20210215 - Rough experimental measures - thick tissue - fine mesh - 5
% % periods
% clear all
% close all
% k=1;
% xLength=120;
% yLength=82.57;
% desiredSideLength=2;
%
% for i=0:desiredSideLength:ceil(xLength)
%     for j=0:desiredSideLength:ceil(yLength)
%         if i==0 || i==ceil(xLength) || j==0 || j==ceil(yLength)
%             x=i;
%             y=j;
%             borderNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% if mod(ceil(xLength),desiredSideLength)~=0
%     for j=0:desiredSideLength:ceil(yLength)
%         x=xLength;
%         y=j;
%         borderNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
% if mod(ceil(yLength),desiredSideLength)~=0
%     for i=0:desiredSideLength:ceil(xLength)
%         x=i;
%         y=yLength;
%         borderNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
% scatter(borderNodes(:,1),borderNodes(:,2),'r')
%
% k=1;
% for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
%     for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
%         x=i;
%         y=j;
%         surfaceNodes(k,:)=[x,y];
%         k=k+1;
%     end
% end
%
%
% if mod(ceil(xLength),desiredSideLength)~=0
%     for j=desiredSideLength:desiredSideLength:ceil(yLength)-desiredSideLength
%         x=floor(xLength/desiredSideLength)*desiredSideLength;
%         y=j;
%         if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%             surfaceNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
%
% if mod(ceil(yLength),desiredSideLength)~=0
%     for i=desiredSideLength:desiredSideLength:ceil(xLength)-desiredSideLength
%         x=i;
%         y=floor(yLength/desiredSideLength)*desiredSideLength;
%         if xLength-x>=desiredSideLength && yLength-y>desiredSideLength
%             surfaceNodes(k,:)=[x,y];
%             k=k+1;
%         end
%     end
% end
%
% hold on
% scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')
%
%
% random=randperm(length(surfaceNodes));
% for i=1:length(random)
%     x_noise=desiredSideLength*(0.5*rand-0.25);
%     y_noise=desiredSideLength*(0.5*rand-0.25);
%     if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<xLength
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
%     else
%         surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
%     end
%     if surfaceNodes(random(i),2)+y_noise>0 && surfaceNodes(random(i),2)+y_noise<yLength
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
%     else
%         surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
%     end
% end
%
% coord=[borderNodes;surfaceNodes];
%
% figure()
% scatter(coord(:,1),coord(:,2))
%
% figure()
% DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
% triplot(DT);
%
% NodeInfo=zeros(length(DT.Points),3);
% NodeInfo(:,1)=DT.Points(:,1);
% NodeInfo(:,2)=DT.Points(:,2);
%
% for i=1:length(NodeInfo)
%     if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
%         NodeInfo(i,3)=1;
%         hold on
%         plot(NodeInfo(i,1),NodeInfo(i,2),'go');
%     end
% end
%
% ConnectionInfo = DT.ConnectivityList;
%
% fileID = fopen('thickTissue_fineMesh_5periods.node','w');
% fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
% for i=1:length(NodeInfo)
%     fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
% end
% fclose(fileID);
%
% fileID = fopen('thickTissue_fineMesh_5periods.ele','w');
% fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
% for i=1:length(ConnectionInfo)
%     fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
% end
% fclose(fileID);
% %%
%% 20210217 Rough experimental measures - thin tissue fine mesh - 1 period
clear all
close all

xLength=24;
yLength=40;
desiredSideLength=2;
symmetricX=0;
symmetricY=0;  %I added this option but don't need to hardcode this, as the FEM model can handle this separately. I will leave it here for now but remove it from the rest of the code that generates symmetric files.

borderNodeCounter=1;

for i=-xLength:desiredSideLength:0
    for j=-yLength:desiredSideLength:0
        if i==0 || i==-xLength || j==0 || j==-yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<0 && j+desiredSideLength>0
        if borderNodes(borderNodeCounter-1,2)<-desiredSideLength/2
            x=i;
            y=0;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=0;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<0 && i+desiredSideLength>0
    if i<-desiredSideLength/2
        for j=-yLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==-mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=-yLength+desiredSideLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=0;
else
    maxSurfaceNodeX=-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=0;
else
    maxSurfaceNodeY=-desiredSideLength;
end

surfaceNodeCounter=1;
for i=-xLength+desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=-yLength+desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise>-xLength && surfaceNodes(random(i),1)+x_noise<0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise>-yLength && surfaceNodes(random(i),2)+y_noise<0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);

if symmetricY==1
    for i=1:length(NodeInfo)
        if NodeInfo(i,1)==0 || NodeInfo(i,1)==-xLength || NodeInfo(i,2)==-yLength
            NodeInfo(i,3)=1;
            hold on
            plot(NodeInfo(i,1),NodeInfo(i,2),'go');
        end
    end
else
    for i=1:length(NodeInfo)
        if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==-xLength || NodeInfo(i,2)==-yLength
            NodeInfo(i,3)=1;
            hold on
            plot(NodeInfo(i,1),NodeInfo(i,2),'go');
        end
    end
end

ConnectionInfo = DT.ConnectivityList;

if symmetricY==1
    fileID = fopen('thinTissue_fineMesh_1period_symmetricY.node','w');
    fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
    for i=1:length(NodeInfo)
        fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
    end
    fclose(fileID);
    
    fileID = fopen('thinTissue_fineMesh_1period_symmetricY.ele','w');
    fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
    for i=1:length(ConnectionInfo)
        fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
    end
    fclose(fileID);
else
    fileID = fopen('thinTissue_fineMesh_1period.node','w');
    fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
    for i=1:length(NodeInfo)
        fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
    end
    fclose(fileID);
    
    fileID = fopen('thinTissue_fineMesh_1period.ele','w');
    fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
    for i=1:length(ConnectionInfo)
        fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
    end
    fclose(fileID);
end
%% 20210217 Rough experimental measures - thin tissue fine mesh - 1 period symmetricY
clear all
close all

xLength=24;
yLength=40/2;
desiredSideLength=2;

borderNodeCounter=1;

for i=-xLength:desiredSideLength:0
    for j=-yLength:desiredSideLength:0
        if i==0 || i==-xLength || j==0 || j==-yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<0 && j+desiredSideLength>0
        if borderNodes(borderNodeCounter-1,2)<-desiredSideLength/2
            x=i;
            y=0;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=0;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<0 && i+desiredSideLength>0
    if i<-desiredSideLength/2
        for j=-yLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==-mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=-yLength+desiredSideLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=0;
else
    maxSurfaceNodeX=-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=0;
else
    maxSurfaceNodeY=-desiredSideLength;
end

surfaceNodeCounter=1;
for i=-xLength+desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=-yLength+desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise>-xLength && surfaceNodes(random(i),1)+x_noise<0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise>-yLength && surfaceNodes(random(i),2)+y_noise<0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);


for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==-xLength || NodeInfo(i,2)==-yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;

fileID = fopen('thinTissue_fineMesh_1period_symmetricY.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thinTissue_fineMesh_1period_symmetricY.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
%% 20210217 Rough experimental measures - thick tissue fine mesh - 1 period
clear all
close all

xLength=24;
yLength=82.57;
desiredSideLength=2;
symmetricX=0;
symmetricY=0;

borderNodeCounter=1;

for i=-xLength:desiredSideLength:0
    for j=-yLength:desiredSideLength:0
        if i==0 || i==-xLength || j==0 || j==-yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<0 && j+desiredSideLength>0
        if borderNodes(borderNodeCounter-1,2)<-desiredSideLength/2
            x=i;
            y=0;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=0;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<0 && i+desiredSideLength>0
    if i<-desiredSideLength/2
        for j=-yLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==-mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=-yLength+desiredSideLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=0;
else
    maxSurfaceNodeX=-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=0;
else
    maxSurfaceNodeY=-desiredSideLength;
end

surfaceNodeCounter=1;
for i=-xLength+desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=-yLength+desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise>-xLength && surfaceNodes(random(i),1)+x_noise<0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise>-yLength && surfaceNodes(random(i),2)+y_noise<0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);

if symmetricY==1
    for i=1:length(NodeInfo)
        if NodeInfo(i,1)==0 || NodeInfo(i,1)==-xLength || NodeInfo(i,2)==-yLength
            NodeInfo(i,3)=1;
            hold on
            plot(NodeInfo(i,1),NodeInfo(i,2),'go');
        end
    end
else
    for i=1:length(NodeInfo)
        if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==-xLength || NodeInfo(i,2)==-yLength
            NodeInfo(i,3)=1;
            hold on
            plot(NodeInfo(i,1),NodeInfo(i,2),'go');
        end
    end
end

ConnectionInfo = DT.ConnectivityList;

if symmetricY==1
    fileID = fopen('thickTissue_fineMesh_1period_symmetricY.node','w');
    fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
    for i=1:length(NodeInfo)
        fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
    end
    fclose(fileID);
    
    fileID = fopen('thickTissue_fineMesh_1period_symmetricY.ele','w');
    fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
    for i=1:length(ConnectionInfo)
        fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
    end
    fclose(fileID);
else
    fileID = fopen('thickTissue_fineMesh_1period.node','w');
    fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
    for i=1:length(NodeInfo)
        fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
    end
    fclose(fileID);
    
    fileID = fopen('thickTissue_fineMesh_1period.ele','w');
    fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
    for i=1:length(ConnectionInfo)
        fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
    end
    fclose(fileID);
end
%% 20210217 Rough experimental measures - thick tissue fine mesh - 1 period symmetricY
clear all
close all

xLength=24;
yLength=82.57/2;
desiredSideLength=2;

borderNodeCounter=1;

for i=-xLength:desiredSideLength:0
    for j=-yLength:desiredSideLength:0
        if i==0 || i==-xLength || j==0 || j==-yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<0 && j+desiredSideLength>0
        if borderNodes(borderNodeCounter-1,2)<-desiredSideLength/2
            x=i;
            y=0;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=0;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<0 && i+desiredSideLength>0
    if i<-desiredSideLength/2
        for j=-yLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==-mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=-yLength+desiredSideLength:desiredSideLength:0
            x=0;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=0;
else
    maxSurfaceNodeX=-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=0;
else
    maxSurfaceNodeY=-desiredSideLength;
end

surfaceNodeCounter=1;
for i=-xLength+desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=-yLength+desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise>-xLength && surfaceNodes(random(i),1)+x_noise<0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise>-yLength && surfaceNodes(random(i),2)+y_noise<0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);


for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==-xLength || NodeInfo(i,2)==-yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;


fileID = fopen('thickTissue_fineMesh_1period_symmetricY.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thickTissue_fineMesh_1period_symmetricY.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
%% 20210217 Rough experimental measures - thin tissue fine mesh - 1 period symmetricY - positive coordinates
clear all
close all

xLength=24;
yLength=40/2;
desiredSideLength=2;

borderNodeCounter=1;

for i=0:desiredSideLength:xLength
    for j=0:desiredSideLength:yLength
        if i==0 || i==xLength || j==0 || j==yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<yLength && j+desiredSideLength>yLength
        if borderNodes(borderNodeCounter-1,2)<yLength-desiredSideLength/2
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<xLength && i+desiredSideLength>xLength
    if i<xLength-desiredSideLength/2
        for j=0:desiredSideLength:yLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=0:desiredSideLength:yLength-desiredSideLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=xLength;
else
    maxSurfaceNodeX=xLength-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=yLength;
else
    maxSurfaceNodeY=yLength-desiredSideLength;
end

surfaceNodeCounter=1;
for i=desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise<xLength && surfaceNodes(random(i),1)+x_noise>0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise<yLength && surfaceNodes(random(i),2)+y_noise>0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);


for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;


fileID = fopen('thinTissue_fineMesh_1period_symmetricY_posCoord.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thinTissue_fineMesh_1period_symmetricY_posCoord.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
%% 20210217 Rough experimental measures - thick tissue fine mesh - 1 period symmetricY - positive coordinates
clear all
close all

xLength=24;
yLength=82.57/2;
desiredSideLength=2;

borderNodeCounter=1;

for i=0:desiredSideLength:xLength
    for j=0:desiredSideLength:yLength
        if i==0 || i==xLength || j==0 || j==yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<yLength && j+desiredSideLength>yLength
        if borderNodes(borderNodeCounter-1,2)<yLength-desiredSideLength/2
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<xLength && i+desiredSideLength>xLength
    if i<xLength-desiredSideLength/2
        for j=0:desiredSideLength:yLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=0:desiredSideLength:yLength-desiredSideLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=xLength;
else
    maxSurfaceNodeX=xLength-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=yLength;
else
    maxSurfaceNodeY=yLength-desiredSideLength;
end

surfaceNodeCounter=1;
for i=desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise<xLength && surfaceNodes(random(i),1)+x_noise>0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise<yLength && surfaceNodes(random(i),2)+y_noise>0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);


for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;


fileID = fopen('thickTissue_fineMesh_1period_symmetricY_posCoord.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thickTissue_fineMesh_1period_symmetricY_posCoord.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
%% 20210217 Rough experimental measures - thin tissue fine mesh - 1 period symmetricY - positive coordinates
clear all
close all

xLength=72;
yLength=40/2;
desiredSideLength=2;

borderNodeCounter=1;

for i=0:desiredSideLength:xLength
    for j=0:desiredSideLength:yLength
        if i==0 || i==xLength || j==0 || j==yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<yLength && j+desiredSideLength>yLength
        if borderNodes(borderNodeCounter-1,2)<yLength-desiredSideLength/2
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<xLength && i+desiredSideLength>xLength
    if i<xLength-desiredSideLength/2
        for j=0:desiredSideLength:yLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=0:desiredSideLength:yLength-desiredSideLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=xLength;
else
    maxSurfaceNodeX=xLength-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=yLength;
else
    maxSurfaceNodeY=yLength-desiredSideLength;
end

surfaceNodeCounter=1;
for i=desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise<xLength && surfaceNodes(random(i),1)+x_noise>0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise<yLength && surfaceNodes(random(i),2)+y_noise>0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);


for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;


fileID = fopen('thinTissue_fineMesh_3periods_symmetricY_posCoord.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thinTissue_fineMesh_3periods_symmetricY_posCoord.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
%% 20210217 Rough experimental measures - thick tissue fine mesh - 1 period symmetricY - positive coordinates
clear all
close all

xLength=72;
yLength=82.57/2;
desiredSideLength=2;

borderNodeCounter=1;

for i=0:desiredSideLength:xLength
    for j=0:desiredSideLength:yLength
        if i==0 || i==xLength || j==0 || j==yLength
            x=i;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
    if j<yLength && j+desiredSideLength>yLength
        if borderNodes(borderNodeCounter-1,2)<yLength-desiredSideLength/2
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        else
            x=i;
            y=yLength;
            borderNodes(borderNodeCounter-1,:)=[x,y];
        end
    end
end

if i<xLength && i+desiredSideLength>xLength
    if i<xLength-desiredSideLength/2
        for j=0:desiredSideLength:yLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    else
        id=find(borderNodes(:,1)==mod(xLength,desiredSideLength));
        borderNodes(id,1)=0;
        for j=0:desiredSideLength:yLength-desiredSideLength
            x=xLength;
            y=j;
            borderNodes(borderNodeCounter,:)=[x,y];
            borderNodeCounter=borderNodeCounter+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')


if mod(xLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeX=xLength;
else
    maxSurfaceNodeX=xLength-desiredSideLength;
end

if mod(yLength,desiredSideLength)>desiredSideLength/2
    maxSurfaceNodeY=yLength;
else
    maxSurfaceNodeY=yLength-desiredSideLength;
end

surfaceNodeCounter=1;
for i=desiredSideLength:desiredSideLength:maxSurfaceNodeX
    for j=desiredSideLength:desiredSideLength:maxSurfaceNodeY
        x=i;
        y=j;
        surfaceNodes(surfaceNodeCounter,:)=[x,y];
        surfaceNodeCounter=surfaceNodeCounter+1;
    end
end

random=randperm(length(surfaceNodes));
for i=1:length(random)
    x_noise=desiredSideLength*(0.5*rand-0.25);
    y_noise=desiredSideLength*(0.5*rand-0.25);
    if surfaceNodes(random(i),1)+x_noise<xLength && surfaceNodes(random(i),1)+x_noise>0
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+y_noise<yLength && surfaceNodes(random(i),2)+y_noise>0
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+y_noise;
    else
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2);
    end
end

hold on
scatter(surfaceNodes(:,1),surfaceNodes(:,2),'b')

coord=[borderNodes;surfaceNodes];

figure()
scatter(coord(:,1),coord(:,2))

figure()
DT = delaunayTriangulation(coord(:,1),coord(:,2)); %2D delaunay triangulation
triplot(DT);

NodeInfo=zeros(length(DT.Points),3);
NodeInfo(:,1)=DT.Points(:,1);
NodeInfo(:,2)=DT.Points(:,2);


for i=1:length(NodeInfo)
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==xLength || NodeInfo(i,2)==yLength
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;


fileID = fopen('thickTissue_fineMesh_3periods_symmetricY_posCoord.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('thickTissue_fineMesh_3periods_symmetricY_posCoord.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);
