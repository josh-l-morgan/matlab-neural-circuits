%function[SPN] = readVast(SPN)

%% find images
SPN = GetMyDir
dSPN = dir(SPN);
iNam = {};
for i =1:length(dSPN)
    if sum(regexp(dSPN(i).name,'.tif')) |...
        sum(regexp(dSPN(i).name,'.png'));
       iNam{length(iNam)+1,1} = dSPN(i).name;
    end
end

%% Read Data
imageInfo = imfinfo([SPN iNam{1}]);
ys = imageInfo.Height;
xs = imageInfo.Width;
zs = length(iNam);

S = zeros(ys,xs,zs);
parfor i = 1:length(iNam)
    S(:,:,i) = imread([SPN iNam{i}]);
end


%% Display data
segNum = max(S(:));
colormap colorcube(256)

for i = 1:zs
   image(S(:,:,i)),pause(.1) 
end


%% Down Sample
dsTest = imresize(S(:,:,1),1/8,'nearest');

[dys dxs] = size(dsTest);
dS = zeros(dys,dxs,zs);
parfor i = 1:zs
   dS(:,:,i) = imresize(S(:,:,i),1/8,'nearest'); 
end

%% render

uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}
for s = 1:segNum
    s
    data = dS==s;
    data = smooth3(data,'box',5);
    p1{s} = patch(isosurface(data,.5), ...
       'FaceColor',uCol{mod(s,6)+1},'EdgeColor','none',...
       'FaceAlpha',.3);
    isonormals(data,p1{s})
    view(3); axis vis3d tight
    camlight; lighting gouraud
    pause(1)
end

%