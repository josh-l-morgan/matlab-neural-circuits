


SPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\export\'
TPN = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\matOut\'

fileName = 'D:\LGNs1\segmentation\VAST\LindaXu\rayAlign01\ColorNames.txt'
ids = readVastColors(fileName);

allS = readVast(SPN);



idNums = [ids{:,1}];

typeLab = [];
for i = 1:length(ids)
    nam = ids{i,2};
       if sum(regexp(lower(nam),'background'))
        typeLab(i) = 0;
       elseif sum(regexp(lower(nam),'rgc'))
        typeLab(i) = 1;
    elseif sum(regexp(lower(nam),'lin'))
        typeLab(i) = 2;
    elseif sum(regexp(lower(nam),'tcp'))
        typeLab(i) = 3;
    elseif sum(regexp(lower(nam),'gli'))
        typeLab(i) = 4;
    else
        typeLab(i) = 5;
    end
    
end



%%
clear lookUpType
lookUpType(idNums+1)= typeLab;
typeS = lookUpType(allS+1);

writeDir = [TPN 'typeField\'];
mkdir(writeDir)
for i = 1:size(typeS,3)
    image(uint8(typeS(:,:,i))*10),pause(.01)
    imwrite(uint8(typeS(:,:,i)),sprintf('%splane_%04.f.tif',writeDir,i)); 
end

%%
for t= 1:max(typeLab)
    subS = allS*0;
    subS(typeS==t) = allS(typeS==t);
    uLab = unique(subS(:));
    
    
end





%%


%{

%%
clf



showS = {[56] [13]}
for s = 1:length(showS)
    
    S = allS * 0;
    useSeg = idNums(showS{s});
    for u = 1:length(useSeg)
        S(allS==useSeg(u)) = 1;
    end
    
    
    
    
    %%%
    
    
    [ys xs zs] = size(S);
    dsTest = imresize(S(:,:,1),1/8,'nearest');
    
    [dys dxs] = size(dsTest);
    dS = zeros(dys,dxs,zs);
    parfor i = 1:zs
        dS(:,:,i) = imresize(S(:,:,i),1/8,'bicubic');
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

