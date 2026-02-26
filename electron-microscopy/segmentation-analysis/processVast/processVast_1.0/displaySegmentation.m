

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'

fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
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

I = uint8(readVast(IPN));









%%


typeLab = [];
for i = 1:length(ids)
    nam = ids{i,2}
    ID = ids{i,1}
    if ID>0;
    if sum(regexp(lower(nam),'background'))
        typeLab(ID) = 0;
    elseif sum(regexp(lower(nam),'rgc'))
        typeLab(ID) = 1;
    elseif sum(regexp(lower(nam),'lin'))
        typeLab(ID) = 2;
    elseif sum(regexp(lower(nam),'tcp'))
        typeLab(ID) = 3;
    elseif sum(regexp(lower(nam),'gli'))
        typeLab(ID) = 4;
    else
        typeLab(ID) = 5;
    end
    end
    
end



typeKey = [' RGC,       Inhibitory,       thalamocortical,       glia,       other'];
[sortType typePos] = sort(typeLab);

cells = double(1:max(allO(:)));


%% color types
typeCol = [0 0 0; 0 1 0 ; 1 0 0; 0 0 1; .75 .75 0; 0 .75 .75]
typeSeg = allO*0;
typeSeg(allO>0) = typeLab(allO(allO>0));
[ys xs zs ] = size(I);
typeI = zeros(ys, xs, 3 ,zs,'uint8');
planeI = typeI(:,:,1,1);
for i = 1:size(typeI,4)
    for c = 1:3
        planeI(:) = typeCol(typeSeg(:,:,i)+1,c);
        %image(planeI*100),pause
        typeI(:,:,c,i) = planeI;
    end
end

for i = 1:size(typeI,4)
    image(typeI(:,:,:,i)*50+repmat(I(:,:,i),[1 1 3])),pause(.1)
end

%%

labSyn = bwlabeln(allS,6);
sProps = regionprops(labSyn,allO,'PixelIdxList','PixelValues');
clear labSyn

%%
% 
% showS = allO;
% 
% for i = 107
%     showS(sProps(i).PixelIdxList) = 1000;%sProps(i).PixelValues;
% end
% 
% 
% showS(sProps(86).PixelIdxList) = 1000;%sProps(i).PixelValues;
% 
% for i = 1:size(showS,3)
%     image(showS(:,:,i)),pause
% end
% 

%% count connections
binOb = [0:obNum];
allSyn = zeros(length(sProps),2);
numSyn = 0;
clear syn
for i = 1:length(sProps)
    vals = sProps(i).PixelValues;
    histVals = hist(vals,binOb);
    histVals = histVals(2:end);
    hits = find(histVals>0);
    
    if length(hits)>=2
        
        hitVal = histVals(hits);
        [sortVal hIdx] = sort(hitVal,'descend');
        
        allSyn(i,:) = [hits(hIdx(1)) hits(hIdx(2))];
        numSyn = numSyn+1;
        syn(numSyn,1:2) = [hits(hIdx(1)) hits(hIdx(2))];
        synId(numSyn) = i;
        if typeLab(hits(hIdx(1))) == 3
            pause
        end
    elseif length(hits) == 1
        allSyn(i,:) = [i i];
        
    else
        allSyn(i,:) = [0 0];
    end
    
end


%% reorder for syn
synTypes = typeLab(syn);
useCells = unique(syn(:));
cellTypes = typeLab(useCells);

[sortCells cellIdx] = sort(cellTypes); %use fopr types
sortIDs = useCells(cellIdx);

lookUp = zeros(1,max(syn(:)));
lookUp(sortIDs) = 1:length(sortIDs);

sortSyn = lookUp(syn); %new sorted synapse matrix


numSyn = length(sortSyn);
numCell = length(useCells);
cells = 1:numCell;
typeKey = [' RGC,       Inhibitory,       thalamocortical,       glia,       other'];


%%  Show connections

conG = zeros(numCell,numCell);
conG = conG + repmat(mod(sortCells+1,2),[length(sortCells) 1]);
conG = conG + repmat(mod(sortCells+1,2),[length(sortCells) 1]);
conG = xor(conG,repmat(mod(sortCells'+1,2),[1 length(sortCells)]));

conG = conG ;

conS = conG*0;

colCon = zeros(numCell,numCell,3,'uint8');

for i = 1:size(sortSyn,1)

    if ~sum(sortSyn(i,:)==0) & abs(diff(sortSyn(i,:)))
        prePos = find(typePos == sortSyn(i,1));
        postPos = find(typePos == sortSyn(i,2));
        conS(sortSyn(i,1),sortSyn(i,2)) = conG(sortSyn(i,1),sortSyn(i,2)) +1;
        image(conS)
    end
end

conG = repmat(conG,[1 1 3]);

%conS = repmat((conS)*100,[1 1 3]);

colCon = uint8(conG*20);
colCon(:,:,3) = colCon(:,:,3) + uint8(conS*100);
background = colCon;
 
image(colCon)
title(typeKey)

%% Assign weights

%Make list of targets and synatic signs
preList = unique(sortSyn(:,1));
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
[dys dxs] = size(dsTest);
dS = zeros(dys,dxs,zs);
uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}


%% Make show list
[sNum outin] = find(syn == 41);
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






