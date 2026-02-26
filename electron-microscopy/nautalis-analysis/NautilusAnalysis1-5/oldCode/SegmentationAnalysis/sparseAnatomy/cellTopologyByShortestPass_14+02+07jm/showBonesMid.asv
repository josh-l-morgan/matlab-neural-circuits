function[col] = showBonesMid(skel,projectDim)


%%
if ~exist('projectDim','var')
    projectDim = 3;
end

dims = [1:3];
dims = dims(dims ~= projectDim);

subs = round(skel.surf.mid);

if ~exist('dims','var')
    dims = [1 2];
end

for i = 1:3
    subs(:,i) = subs(:,i) - skel.offset(i) + 1;
end

nodeIm  = zeros(skel.fsize(dims(1)),skel.fsize(dims(2)));

nodeInd = sub2ind(size(nodeIm),subs(:,dims(1)),subs(:,dims(2)));
nodeIm(nodeInd) = 1;
image(nodeIm*100)


%% draw edges
edgeIm = nodeIm * 0;    

for b = 1:length(skel.bones)
    
    edges = skel.bones(b).edges;
drawNum = 10;
steps = [0:drawNum]/drawNum;
stepE = zeros(length(steps),3);

for i = 1:size(edges,1);
    startE = subs(edges(i,1),:);
    stopE = subs(edges(i,2),:);
    difE = [(stopE(1)-startE(1))  (stopE(2)-startE(2)) ...
        (stopE(3)-startE(3))];
    %lengthE = sqrt(sum(difE.^2));
    
    stepE(:,1) = startE(1) + difE(1) * steps ;
    stepE(:,2) = startE(2) + difE(2) * steps ;
    stepE(:,3) = startE(3) + difE(3) * steps ;
    stepE = round(stepE);
    
    indE = sub2ind(size(edgeIm),stepE(:,dims(1)),stepE(:,dims(2)));
    for s = 1:length(indE)
        edgeIm(indE(s)) = edgeIm(indE(s)) + 1;
    end
end

end
image(edgeIm*100),pause(.01)

%% Draw bridges

bridgeIm = nodeIm * 0;
edges = skel.bridges;
drawNum = 100;
steps = [0:drawNum]/drawNum;
stepE = zeros(length(steps),3);
for i = 1:size(edges,1);
    startE = subs(edges(i,1),:);
    stopE = subs(edges(i,2),:);
    difE = [(stopE(1)-startE(1))  (stopE(2)-startE(2)) ...
        (stopE(3)-startE(3))];
    %lengthE = sqrt(sum(difE.^2));
    
    stepE(:,1) = startE(1) + difE(1) * steps ;
    stepE(:,2) = startE(2) + difE(2) * steps ;
    stepE(:,3) = startE(3) + difE(3) * steps ;
    stepE = round(stepE);
    
    indE = sub2ind(size(edgeIm),stepE(:,dims(1)),stepE(:,dims(2)));
    for s = 1:length(indE)
        bridgeIm(indE(s)) = bridgeIm(indE(s)) + 1;
    end
    
end

image(bridgeIm*100)

%%

SE = strel('disk',3);
bridgeIm = imdilate(bridgeIm,SE);
bridgeIm = bridgeIm * 200;

SE = strel('disk',2);
nodeIm = imdilate(nodeIm,SE);
nodeIm = nodeIm * 300;

SE = strel('disk',1);
edgeIm = imdilate(edgeIm,SE);
edgeIm = edgeIm * 150/median(edgeIm(edgeIm>0));

col(:,:,2) = edgeIm + bridgeIm - nodeIm/2;
col(:,:,1) = nodeIm;
col(:,:,3) = nodeIm*0;

image(uint8(col)),pause(.010)



