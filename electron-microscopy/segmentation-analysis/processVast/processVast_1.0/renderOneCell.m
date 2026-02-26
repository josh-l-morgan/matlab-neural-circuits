

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export2\'
%TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
OPN = GetMyDir
TPN = GetMyDir
%SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
%IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'


dOPN = dir(OPN);
for i = 1:length(dOPN)
    if sum(regexp(dOPN(i).name,'.txt'))
        fileName = [OPN dOPN(i).name];
    end
end
        
%fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
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

%allS = uint8(readVast(SPN));
%I = uint8(readVast(IPN));


%%

displayName = 'TCP 19'
foundCell = 0;
displayID = [];
for i = 1:length(ids)
    nam = ids{i,2};
    ID = ids{i,1};
    if ID>0;
    if sum(regexp(nam,displayName))
        foundCell = foundCell+1;
        displayID(foundCell) = ID;
    end
    end
    
end



%%
maxCells = 1;
if maxCells
   
    sumCell = zeros(size(allO,1),size(allO,2),3,'uint8');
    for s = 1:3;%length(displayID)
       tempSum = sum(allO == displayID(s),3); 
       tempSum = tempSum * 256/double(max(tempSum(:)));
        sumCell(:,:,s) = tempSum;
    end
    
    


image(sumCell)
imwrite(sumCell,[TPN 'sumCell.tif'])
end

%%

renderCells = 0;
if showCells
clf


for s = 1:length(displayID)
    
    O = allO * 0;
    useSeg = displayID(s);
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
        dS = smooth3(dS,'box',3);

    %%% render
    
    uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}
    
    data = dS>.1;
    %data = smooth3(data,'box',1);
    p1(s) = patch(isosurface(data,.1), ...
        'FaceColor',uCol{mod(s,6)+1},'EdgeColor','none');
    isonormals(data,p1(s))
    view(2); 
    axis([1 dxs 1 dys 1 zs])
    axis tight
    axis square
    
    camlight; lighting gouraud
    pause(1)
end
end

axis normal

% 
% for i = 1:360
%     
%     view
% end

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






