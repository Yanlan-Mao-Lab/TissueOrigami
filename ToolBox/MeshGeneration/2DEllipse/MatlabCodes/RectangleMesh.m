clear all
close all
k=1;
maxLength=8;
minLength=4;
desiredSideLength=2;

for i=0:desiredSideLength:maxLength
    for j=0:desiredSideLength:minLength
        if i==0 || i==maxLength || j==0 || j==minLength
            x=i;
            y=j;
            borderNodes(k,:)=[x,y];
            k=k+1;
        end
    end
end

scatter(borderNodes(:,1),borderNodes(:,2),'r')

k=1;
for i=desiredSideLength:desiredSideLength:maxLength-desiredSideLength
    for j=desiredSideLength:desiredSideLength:minLength-desiredSideLength
        x=i;
        y=j;
        surfaceNodes(k,:)=[x,y];
        k=k+1;
    end
end

if mod(ceil(maxLength),desiredSideLength)~=0
    for j=desiredSideLength:desiredSideLength:ceil(minLength)-desiredSideLength
            x=floor(maxLength/desiredSideLength)*desiredSideLength;
            y=j;
            if maxLength-x>=desiredSideLength && minLength-y>desiredSideLength
                surfaceNodes(k,:)=[x,y];
                k=k+1;
            end
    end
end


if mod(ceil(minLength),desiredSideLength)~=0
    for i=desiredSideLength:desiredSideLength:ceil(maxLength)-desiredSideLength
            x=i;
            y=floor(minLength/desiredSideLength)*desiredSideLength;
            if maxLength-x>=desiredSideLength && minLength-y>desiredSideLength
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
    if surfaceNodes(random(i),1)+x_noise>0 && surfaceNodes(random(i),1)+x_noise<maxLength
        surfaceNodes(random(i),1)=surfaceNodes(random(i),1)+x_noise;
    else
        surfaceNodes(random(i),1)=surfacceNodes(random(i),1);
    end
    if surfaceNodes(random(i),2)+x_noise>0 && surfaceNodes(random(i),2)+x_noise<minLength
        surfaceNodes(random(i),2)=surfaceNodes(random(i),2)+x_noise;
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
    if NodeInfo(i,1)==0 || NodeInfo(i,2)==0 || NodeInfo(i,1)==maxLength || NodeInfo(i,2)==minLength
        %     if NodeInfo(i,1)==0 || NodeInfo(i,1)==10 || NodeInfo(i,2)==10
        NodeInfo(i,3)=1;
        hold on
        plot(NodeInfo(i,1),NodeInfo(i,2),'go');
    end
end


ConnectionInfo = DT.ConnectivityList;

fileID = fopen('smallRectangle.node','w');
fprintf(fileID,'%d %d %d %d\n',[length(NodeInfo),2,0,1]); %number of nodes
for i=1:length(NodeInfo)
    fprintf(fileID,'%d %.4f %.4f %d\n',[i-1,NodeInfo(i,:)]);
end
fclose(fileID);

fileID = fopen('smallRectangle.ele','w');
fprintf(fileID,'%d %d %d\n',[length(ConnectionInfo),3,0]); %number of connections
for i=1:length(ConnectionInfo)
    fprintf(fileID,'%d %d %d %d \n',[i-1, ConnectionInfo(i,1)-1,ConnectionInfo(i,2)-1,ConnectionInfo(i,3)-1]);
end
fclose(fileID);