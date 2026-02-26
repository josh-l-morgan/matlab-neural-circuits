function[col] = showBones3D(skel,dSamp)

%%


if ~exist('dSamp','var')
    dSamp = 10;
end

%%Get Subs
subs = round(skel.surf.mid);

%%reshape subs
for i = 1:3
    subs(:,i) = ceil(subs(:,i)/dSamp - skel.offset(i)/dSamp);
end

%%Make Volume
skelVol = zeros(ceil(skel.fsize/dSamp),'uint8');





%% draw edges

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
        
        indE = sub2ind(size(skelVol),stepE(:,1),stepE(:,2),stepE(:,3));
        skelVol(indE) = 1;
    end
    
end


%isosurface(skelVol,0), axis equal,view(3), camlight, lighting none,

%% Add nodes
nodeInd = sub2ind(size(skelVol),subs(:,1),subs(:,2),subs(:,3));
skelVol(nodeInd) = 3;


%% Draw bridges

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
    
     indE = sub2ind(size(skelVol),stepE(:,1),stepE(:,2),stepE(:,3));
        skelVol(indE) = 2;
    
end

%clf
%isosurface(skelVol,'FaceColor','red'), view(3), camlight, lighting flat, 





%%
clf
[faces,verts,colors] = isosurface(skelVol,.1,double(skelVol)); 
% patch('Vertices', verts, 'Faces', faces, ... 
%     'FaceVertexCData', [colors*0  colors * 1 colors * 0], 'BackFaceLighting', 'lit')
% patch('Vertices', verts, 'Faces', faces,  'BackFaceLighting', 'lit',...
%     'SpecularStrength',10,'FaceLighting','flat',...
%     'FaceVertexCData', colors)

patch('Vertices', verts, 'Faces', faces, ... 
    'FaceVertexCData', colors, ... 
    'FaceColor','flat', ... 
    'edgecolor', 'flat');

camlight headlight

axis square


cmap = [ 0 0 0; 1 0 0; 0 1 0; 0 0 1]
colormap(cmap)

%%
% 
% SE = strel('ball',3,3);
% skelVol = imdilate(skelVol,SE);