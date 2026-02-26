

clear all
OPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
%TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'
%OPN = GetMyDir
%TPN = GetMyDir
SPN = 'D:\LGNs1\segmentation\VAST\Joshm\rayAlign1Synapse\eportSyn\'
%IPN = 'D:\LGNs1\HP_processing\rayAlign01\images\'


dOPN = dir([OPN '*.txt']);
        fileName = [OPN dOPN(1).name];
        
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

allS = uint8(readVast(SPN));
%I = uint8(readVast(IPN));


%%

displayName = '78'
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
maxCells = 0;
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

renderCells = 1;
if renderCells
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
%% Grab synapses

renderCells = 1;
if renderCells

[synLab numSyn] = bwlabeln(allS,26);
useSyn = [];
for s = 1:length(displayID)
    
    O = allO * 0;
    useSeg = displayID(s);
    for u = 1:length(useSeg)
        
      useSyn = cat(2,useSyn, unique(synLab(allO==useSeg(u))));
    end
end
    
useSyn = useSyn(useSyn>0);
useSyn = unique(useSyn);

useS = allS * 0;
for i = 1:length(useSyn)
   useS(synLab == useSyn(i)) = 1; 
end

clear synLab

%%
s = length(p1)+1;
  
    [ys xs zs] = size(O);
    dsTest = imresize(useS(:,:,1),1/8,'nearest');

    [dys dxs] = size(dsTest);
    dS = zeros(dys,dxs,zs);
    parfor i = 1:zs
        dS(:,:,i) = imresize(useS(:,:,i),1/8,'bicubic');
    end
        dS = smooth3(dS,'box',3);

    %%% render
    
    uCol = {'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'}
    
    data = dS>.1;
    %data = smooth3(data,'box',1);
    p1(s) = patch(isosurface(data,.1), ...
        'FaceColor','r','EdgeColor','none','FaceAlpha',.6);
    isonormals(data,p1(s))
    view(2); 
    axis([1 dxs 1 dys 1 zs])
   
    
    camlight; lighting gouraud
    pause(1)

axis normal

%%
end



