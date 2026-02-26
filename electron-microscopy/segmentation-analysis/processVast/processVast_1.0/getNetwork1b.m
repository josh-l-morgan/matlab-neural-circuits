

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'

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
    elseif sum(regexp(lower(nam),'xun'))
        typeLab(ID) = 4;
    else
        typeLab(ID) = 5;
    end
    end
    
end



typeKey = [' RGC,       Inhibitory,       thalamocortical,       glia,       other'];
[sortType typePos] = sort(typeLab);

cells = double(1:max(allO(:)))

%%

labSyn = bwlabeln(allS,6);
sProps = regionprops(labSyn,allO,'PixelIdxList','PixelValues');

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


synTypes = typeLab(syn);

datCon1.synTypes = synTypes;
datCon1.syn = syn;
datCon1.typeLab = typeLab;

%%  Show connections

getPos = 1:max(typePos);
getPos(typePos)
conG = zeros(max(typePos),max(typePos));
conG = conG + repmat(mod(sortType,2),[length(sortType) 1]);
conG = conG + repmat(mod(sortType,2),[length(sortType) 1]);
conG = xor(conG,repmat(mod(sortType',2),[1 length(sortType)]));

conG = conG * 40;
conG = conG * 0;
for i = 1:size(syn,1)

    if ~sum(syn(i,:)==0) & abs(diff(syn(i,:)))
        syn(i,:)

        prePos = find(typePos == syn(i,1));
        postPos = find(typePos == syn(i,2));
       % if max([typeLab(prePos) typeLab(postPos)])<=3
        conG(prePos,postPos) = conG(prePos,postPos)+100;
       % end
    end
end


image(conG)
image(conG(1:334,1:100))
title(typeKey)

%% Activate connections

%Make list of targets and synatic signs
preList = unique(syn(:,1));
for i = cells;
   ID = i; 
   if typeLab(ID)==1
       sign = 1;
   elseif typeLab(ID) == 2;
       sign = -1;
   else 
       sign = 0;
   end
   
   targs = syn(syn(:,1) ==ID,2);
   weights(i,:) = (hist(targs,cells))*sign;
   
end

%% Activate


for p = 1:length(preList)
    active = weights * 0; 
    stim = active;
    activate = preList(p);
    stim(activate,:) = 1;
    active = active + stim;   

    for r = 1:3
       %show active
    image(active*30+128),
    pause(.1)
    %spread active
    activate = active.*weights;
    image(activate*30+128)
    pause(.1)
    summed = sum(activate,1); %sum all inputs
    active = repmat(summed',[1 length(summed)]);
    end
    pause(1)
    
end




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



showS = {[56]}
for s = 1:length(showS)
    
    O = allO * 0;
    useSeg = showS{s};
    for u = 1:length(useSeg)
        O(allO==useSeg(u)) = 1;
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
    view(3); axis vis3d tight
    camlight; lighting gouraud
    pause(1)
end
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

