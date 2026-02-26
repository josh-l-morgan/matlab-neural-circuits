

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\exportSyn\'

dOPN = dir([OPN '*.txt']);
        fileName = [OPN dOPN(1).name];
ids = readVastColors(fileName);

allO = readVast(OPN);
obNum = max(allO(:));
if obNum<256
    allO = uint8(allO);
elseif obNum<(2^16)
    allO = uint16(allO);
else
    allO = single(allO);
end


allS = uint8(readVast(SPN));

%%
cells = double(1:max(allO(:)))
numCell = length(cells);

[typeLab typeKey] = getTypes(ids);
typeLab = typeLab(cells);
[sortType typePos] = sort(typeLab);



%% find pre and post cell for each synapse
[syn sProps] = getPrePost(allS,allO);

for i = 1:size(syn,1)
    conDat.syn(i,1) = find(typePos == syn(i,1));
    conDat.syn(i,2) = find(typePos == syn(i,2));
end
%conDat.syn = typePos(syn);
conDat.sProps = sProps;
conDat.type = sortType;
conDat.cellID = typePos;



%% graphCon
colCon = graphCon(syn, typeLab)
background = colCon;

%pause

%% Assign weights

%Make list of targets and synatic signs
weights = zeros(numCell,numCell);
for i = cells;
   ID = i; 
   if sortCells(ID)==1
       sign = 1;
   elseif sortCells(ID) == 2;
       sign = -1;
   else 
       sign = 0;
   end
   
   targs = sortSyn(sortSyn(:,1) ==ID,2);
   weights(i,:) = (hist(targs,cells))*sign;
   
end
preList = unique(sortSyn(:,1));


%% Activate


exciteCells = find(sortCells==1);
background = double(background);


frameNum = 0;    
  clear recMov
for p = 1:length(exciteCells)
    active = weights * 0; 
    stim = active;
    stim(exciteCells(p),:) = 1;
    active = active + stim;  
    
            colCon = uint8(background);

        frameNum = frameNum+1;
    recMov(:,:,:,frameNum) = colCon;
            frameNum = frameNum+1;
    recMov(:,:,:,frameNum) = colCon;
    
     
    
    colCon = colCon + uint8(repmat(stim*1000,[1 1 3]));
%     image(colCon)
%     pause(.2)
    
        frameNum = frameNum+1;
    recMov(:,:,:,frameNum) = colCon;
        colCon = uint8(background);

    for r = 1:3
    colCon(:,:,1) = background(:,:,1) - active * 100;
    colCon(:,:,2) = background(:,:,2)  + active * 100;
%     image(colCon),pause(.2)
      frameNum = frameNum+1;
    recMov(:,:,:,frameNum) = colCon;
    activate = active.*weights;
    colCon(:,:,1) = background(:,:,1) - activate * 100;
    colCon(:,:,2) = background(:,:,1) + activate * 100;
%     image(colCon)
%     pause(.2)
      frameNum = frameNum+1;
    recMov(:,:,:,frameNum) = colCon;
    summed = sum(activate,1); %sum all inputs
    active = repmat(summed',[1 length(summed)]);
    end
    %pause(.5)
    
end
% 
% for i = 1:size(recMov,4)
%     image(recMov(:,:,:,i)),pause(.2)
% end
% 
% waveDir = [TPN 'wave\'];
% mkdir(waveDir);
% for i = 1:size(recMov,4)
%     tifName = sprintf('%swave_%05.0f.tif',waveDir,i);
%     imwrite(recMov(:,:,:,i),tifName);
%     
% end

%% resize
scaleF = 4;
recMov2 = zeros(size(recMov,1)*scaleF,size(recMov,2)*scaleF,3,size(recMov,3),'uint8');
for i = 1:size(recMov,4)
    for c = 1:3
       recMov2(:,:,c,i) = imresize(recMov(:,:,c,i),scaleF);
        
    end   
    
end


%% Make video
    aviName = sprintf('%swdsfaveM2.avi',TPN);

myVideo = VideoWriter(aviName,'Uncompressed AVI');
%myVideo = VideoWriter(aviName, 'Uncompressed AVI');
myVideo.FrameRate = 5;  % Default 30

%myVideo.Quality = 100;    % Default 75
open(myVideo);
writeVideo(myVideo, recMov2);
close(myVideo);
%%
%{
%% Make type image
clear lookUpType
lookUpType(idNums+1)= typeLab;
typeS = lookUpType(allO+1);

writeDir = [TPN 'typeField\'];
mkdir(writeDir)
for i = 1:size(typeS,3)
    image(uint8(typeS(:,:,i))*10),pause(.01)
    imwrite(uint8(typeS(:,:,i)),sprintf('%splane_%04.f.tif',writeDir,i));
end

%%
for t= 1:max(typeLab)
    subS = allO*0;
    subS(typeS==t) = allO(typeS==t);
    uLab = unique(subS(:));
    
    
end


%}


%%




%%

showCells = 0;
if showCells
clf



showS = {[1]}
for s = 1:length(showS)
    
    O = allO * 0;
    useSeg = showS{s};
    for u = 1:length(useSeg)
        O(allS==useSeg(u)) = 1;
    end
    
    [ys xs zs] = size(O);
    dsTest = imresize(O(:,:,1),1/8,'nearest');
    
    [dys dxs] = size(dsTest);
    dS = zeros(dys,dxs,zs);
    parfor i = 1:zs
        dS(:,:,i) = imresize(O(:,:,i),1/8,'bicubic');
    end
    
    %%% render
    
    uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}
    
    data = dS>.1;
    %data = smooth3(data,'box',1);
    p1(s) = patch(isosurface(data,.1), ...
        'FaceColor',uCol{mod(s,6)+1},'EdgeColor','none');
    isonormals(data,p1(s))
    view(2); 
    axis vis3d tight
    camlight; lighting gouraud
    pause(1)
end
end

for i = 1:360
    
    view
end

%}
%%

%{

%% render all

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
    
%}

%% Render back from seed

seed = 41;

showSeed = 1;
if showSeed
%% set up objects
    clf
Ob = logical(allS*0);
[ys xs zs] = size(Ob);
dsOb = imresize(Ob(:,:,1),1/8,'nearest');
[dys dxs] = size(dsOb);
dS = zeros(dys,dxs,zs);
uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}


%% Make show list
[sNum inout] = find(syn == 41);
out = sNum(inout==1);
in = sNum(inout == 2);

outID = synID(out);
inID = synID(in);

inputs = syn(in,1);
outputs = syn(out,2);



%% Show seed 
Ob = Ob * 0;
    for i = 1:length(in)
    Ob(allO==seed) = 1;
end
parfor i = 1:zs
        dS(:,:,i) = imresize(Ob(:,:,i),1/8,'bicubic');
end
 p1(1) = patch(isosurface(dS,.1), ...
        'FaceColor','magenta','EdgeColor','none');    

    view(2); 
    axis vis3d image
    camlight; 
    lighting gouraud
    pause(.1)


%% Show synapses


Ob = Ob * 0;
%---Label---%

if ~isempty(in)
for i = 1:length(in)
    Ob(sProps(synId(in(i))).PixelIdxList) = 1;
end

parfor i = 1:zs
        dS(:,:,i) = imresize(Ob(:,:,i),1/8,'bicubic');
end
 p1(2) = patch(isosurface(dS,.1), ...
        'FaceColor','yellow','EdgeColor','none');
end
   
if ~isempty(out)
Ob = Ob * 0;
    for i = 1:length(out)
    Ob(sProps(synId(out(i))).PixelIdxList) = 1;
end
parfor i = 1:zs
        dS(:,:,i) = imresize(Ob(:,:,i),1/8,'bicubic');
end
 p1(3) = patch(isosurface(dS,.1), ...
        'FaceColor','yellow','EdgeColor','none');    
end
    view(2); 
    axis vis3d image
    camlight right; 
    lighting gouraud
    
    pause(.1)
    
    
    %% show partner cells
    
    typeCol = {'green','red','blue','cyan','cyan','cyan','cyan'}
    
    contacted = [inputs outputs];
    
    contactType = typeLab(contacted);
    uniqueTypes = unique(contactType);
    for t = 1:length(uniqueTypes)
        useContacted = contacted(contactType == uniqueTypes(t));
        
        Ob = Ob * 0;
        for i = 1:length(useContacted);
            Ob(allO==useContacted(i)) = 1;
        end
        
        parfor i = 1:zs
            dS(:,:,i) = imresize(Ob(:,:,i),1/8,'bicubic');
        end
        p1(1) = patch(isosurface(dS,.1), ...
            'FaceColor',typeCol{uniqueTypes(t)},'EdgeColor','none');
        
        view(2);
        axis vis3d image
        camlight right;
        lighting gouraud
        pause(.1)
    end %run all inputs
    
    
    set(p1(1),'backfacelighting','unlit');

    

   %% Set viewing
    %isonormals(dS,p1(s))
    %
    
    
    
    
    
    
    

for s = 1:length(showS)
    
    
    useSeg = showS{s};
    for u = 1:length(useSeg)
        O(allS==useSeg(u)) = 1;
    end
    
    [ys xs zs] = size(O);
    dsOb = imresize(Ob(:,:,1),1/8,'nearest');
    
    [dys dxs] = size(dsTest);
    dS = zeros(dys,dxs,zs);
    parfor i = 1:zs
        dS(:,:,i) = imresize(O(:,:,i),1/8,'bicubic');
    end
    
    %%% render
    
    uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}
    
    data = dS>.1;
    %data = smooth3(data,'box',1);
    p1(s) = patch(isosurface(data,.1), ...
        'FaceColor',uCol{mod(s,6)+1},'EdgeColor','none');
    isonormals(data,p1(s))
    view(2); 
    axis vis3d tight
    camlight; lighting gouraud
    pause(1)
end
end

for i = 1:360
    
    view
end






