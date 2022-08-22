close all
clear all
% input parameters
r = 5;                    %sphere radius
desiredSideLength = 3;      %desired side length 

nNodesPerSide = ceil((2*pi*r)/desiredSideLength)+1;
desiredSideAngle = (2*pi)/(nNodesPerSide-1);
j=1;
for i=1:nNodesPerSide-2
    x_middleNodes(j,1)= r * cos(desiredSideAngle*i);
    y_middleNodes(j,1)= r * sin(desiredSideAngle*i);
    z_middleNodes(j,1)= 0;
    j=j+1;
end

x_middleNodes=[x_middleNodes;r];
y_middleNodes=[y_middleNodes;0];



j=1;
for i=1:nNodesPerSide-2
    arcRadius= x_middleNodes(i,1);
    arcLength= (2* pi)*arcRadius;
    nNodesPerArc=ceil(arcLength/desiredSideLength)+1;
    
    foundDesiredLength = false;
    while (~foundDesiredLength) && (nNodesPerArc>0)
        effectiveArcLength=arcLength/nNodesPerArc;
        if(effectiveArcLength < 0.75*desiredSideLength)
            nNodesPerArc = nNodesPerArc -1;
        else
            foundDesiredLength=true;
        end
    end
    if nNodesPerArc>0
        desiredArcAngle = (2* pi)/(nNodesPerArc);
        for i=1:(nNodesPerArc)
            x_surfaceNodes(j,1)= arcRadius * cos(desiredArcAngle*i);
            y_surfaceNodes(j,1)= arcRadius * sin(desiredArcAngle*i);
            j=j+1;
        end
    end
end

surfaceNodes_all=[round(x_surfaceNodes,2) round(y_surfaceNodes,2)];
surfaceNodes_all_unique=unique(surfaceNodes_all,'rows');

random=randperm(length(surfaceNodes_all_unique));
for i=1:length(random)
    x_noise=desiredSideLength*(0.02*rand-0.01);
    y_noise=desiredSideLength*(0.02*rand-0.01);
%     if x_surfaceNodes(random(i))+x_noise>0
%         x_rand(i,1)=x_surfaceNodes(random(i))+x_noise;
%     else
%         x_rand(i,1)=x_surfaceNodes(random(i));
%     end
%     if y_surfaceNodes(random(i))+y_noise>0
%         y_rand(i,1)=y_surfaceNodes(random(i))+y_noise;
%     else
%         y_rand(i,1)=y_surfaceNodes(random(i));
%     end
    x_rand(i,1)=surfaceNodes_all_unique(i,1)+x_noise;
    y_rand(i,1)=surfaceNodes_all_unique(i,2)+y_noise;
end


x = [0;x_middleNodes;x_rand];
y = [0;y_middleNodes;y_rand];
z = sqrt(abs((r.^2)-(x.^2)-(y.^2)));

NodeInfo=[x,y,z];
radius=desiredSideLength;
[triHull, vbOutside, vbInside] = AlphaHull(NodeInfo,radius)
plot3(NodeInfo(:,1),NodeInfo(:,2),NodeInfo(:,3),'ro');
hold on;
fprintf(' * Found %d points for the outer surface\n',sum(vbOutside));
if (sum(vbOutside) > 0)
    trisurf(triHull(vbOutside,:),NodeInfo(:,1),NodeInfo(:,2),NodeInfo(:,3),...
        'FaceColor','cyan','FaceAlpha',0.5)
end

%% Generate a file with connectivity data
% [number of nodes]
% [node ID (start indexing from 0] [x] [y] [z] [is it at the border? 0 for No, 1 for yes. In case of sphere, it is all 0 ]
% ...
% [number of connections]
% [corner0] [corner1] [corner2] .  %start indexing from 0
% ...

for i=1:length(NodeInfo)
    if abs(NodeInfo(i,1))<0.01*desiredSideLength
        NodeInfo(i,1)=0;
    elseif abs(NodeInfo(i,2))<0.01*desiredSideLength
        NodeInfo(i,2)=0;
    elseif abs(NodeInfo(i,3))<0.01*desiredSideLength
        NodeInfo(i,3)=0;
    end
end

for i=1:length(NodeInfo)
    %if ((x(i)==0 && y(i)==0)||(x(i)==0 && z(i)==0)||(y(i)==0 && z(i)==0))
    if  NodeInfo(i,3)==0
        NodeInfo(i,4)=1;
    end
end
ConnectionInfo = triHull(vbOutside,:);

fileID = fopen('SphericalTriangulation','w');
fprintf(fileID,'%d\n',length(NodeInfo)); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fprintf(fileID,'%d\n',length(ConnectionInfo)); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d \n',[ConnectionInfo(i,2)-1,ConnectionInfo(i,1)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);



%%
for i=1:length(ConnectionInfo)
    node1=ConnectionInfo(i,1);
    node2=ConnectionInfo(i,2);
    node3=ConnectionInfo(i,3);
    sideLength1=sqrt((NodeInfo(node1,1)-NodeInfo(node2,1))^2+(NodeInfo(node1,2)-NodeInfo(node2,2))^2+(NodeInfo(node1,3)-NodeInfo(node2,3))^2);
    sideLength2=sqrt((NodeInfo(node1,1)-NodeInfo(node3,1))^2+(NodeInfo(node1,2)-NodeInfo(node3,2))^2+(NodeInfo(node1,3)-NodeInfo(node3,3))^2);
    sideLength3=sqrt((NodeInfo(node2,1)-NodeInfo(node3,1))^2+(NodeInfo(node2,2)-NodeInfo(node3,2))^2+(NodeInfo(node2,3)-NodeInfo(node3,3))^2);
    sideLength(i,:)=[sideLength1,sideLength2,sideLength3];
end
figure()
histogram(sideLength);
%%
% NodesAtBorder(:,1)=[0;r/sqrt(2);r;0;x_middleNodes];
% NodesAtBorder(:,2) = [0;r/sqrt(2);0;r;y_middleNodes];
% NodesAtBorder(:,3)= sqrt(abs((r.^2)-((NodesAtBorder(:,1)).^2)-((NodesAtBorder(:,2)).^2)));
% length(NodesAtBorder)
% length(find(NodeInfo(:,4)==1))
% %scatter3(NodeInfo(:,1),NodeInfo(:,2),NodeInfo(:,3))
% scatter3(NodesAtBorder(:,1),NodesAtBorder(:,2),NodesAtBorder(:,3),'r')
% hold on
% borderID=find(NodeInfo(:,4)==1);
% scatter3(NodeInfo(borderID,1),NodeInfo(borderID,2),NodeInfo(borderID,3),'b')