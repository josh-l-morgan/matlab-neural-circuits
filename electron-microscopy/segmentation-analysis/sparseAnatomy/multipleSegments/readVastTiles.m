
SPN = 'C:\Users\joshm\Documents\myWork\LGNs1\jlmHomeSeg\export1\';
TPN = 'C:\Users\joshm\Documents\myWork\LGNs1\jlmHomeSeg\export1Mat\';
mkdir(TPN)

dSPN = dir([SPN '*.png']);

inams = {dSPN.name};


%% Parse images
clear r c s
for i = 1:length(inams)
    nam = inams{i};
    
    Yp = regexp(nam,'_Y');
    Xp = regexp(nam,'_X');
    sp = regexp(nam,'_s');
    dotp = regexp(nam,'.png');
    
    r(i) = str2num(nam(Yp(end)+2:Xp(end)-1));
    c(i) = str2num(nam(Xp(end)+2:dotp(end)-1));
    s(i) = str2num(nam(sp(end)+2:Yp(end)-1))+1;
end


%% find locations

imageInfo = imfinfo([SPN inams{1}]);
width = imageInfo.Width;
height = imageInfo.Height;
offset = [(r-1)*height; (c-1)*width]';

maxDim = [max(offset,[],1) + [height width] max(s)];

%% Read in
%%Create o , in which o{i} lists all the pixels y x z for object number i

clear o
colormap colorcube(256)
bufOb = 100; %factor to add to object holder
objSize = zeros(10000,1); %start size for all objects;
for i = 1:length(inams)
    disp(sprintf('Grabing objects from tile %d of %d',i,length(inams)))
    I = imread([SPN inams{i}]);
    inds = find(I>0);
    [y x ] = ind2sub([height width],inds);
    v = I(inds);
    uniqueObjects = unique(v);
    for u = 1:length(uniqueObjects)
        oId = uniqueObjects(u);
        objPos = find(v==oId);
        objTileSubs = [y(objPos)+offset(i,1) x(objPos) + offset(i,2)...
            ones(length(objPos),1)* s(i)];
        objSize(oId) = objSize(oId) + length(objPos);
        addSize = length(objPos);
        
        try
           oSize =  size(o{oId});
        catch err
            o{oId} = zeros(bufOb,3);
            oSize = size(o{oId});
        end
        
        if objSize(oId)>oSize(1)
            o{oId}(oSize(1)+addSize*bufOb,1:3) = 0;
        end
            
        o{oId}(objSize(oId)- addSize +1:objSize(oId),1:3) = objTileSubs;
    end
end

for i = 1:length(o) %% clip off remaining zeros
    if objSize(i)
        o{i} = o{i}(1:objSize(i),1:3);
    end
end


save([TPN 'o.mat'],'o')

%% Create sum image
downSamp = 5;
colormap gray(256)
xyDim = ceil([maxDim(1) maxDim(2)]/downSamp);
allSum = zeros([xyDim length(o)]);

for i = 1:length(o)
    disp(sprintf('Getting sums of object %d of %d',i,length(o)))
     subs = o{i};
     tempSum = zeros(xyDim);
     if ~isempty(subs)
         xyInd = sub2ind(size(tempSum),ceil(subs(:,1)/downSamp),ceil(subs(:,2)/downSamp));
         histVal = hist(xyInd,[.5:1:length(tempSum(:))]);
         tempSum(:) =histVal;
         allSum(:,:,i) = tempSum;
%          image(tempSum*256/max(tempSum(:)))
%          pause(.001)
     end
end

%%combine colors

getText = dir([SPN '*.txt']);
ids = readVastColors([SPN getText(1).name]);

typeLab = [];
for i = 1:length(ids)
    nam = ids{i,2};
    ID = ids{i,1};
    if ID>0;
    if sum(regexp(lower(nam),'background'))
        colTab(ID,:) = [0 0 0];
    elseif sum(regexp(lower(nam),'rgc'))
        colTab(ID,:) = [0 1 0] + rand(1,3)*.2;
    elseif sum(regexp(lower(nam),'fil'))
        colTab(ID,:) = [1 1 0]  + rand(1,3)*.2;
    elseif sum(regexp(lower(nam),'unk'))
        colTab(ID,:) = [1 0 0]  + rand(1,3)*.2;
    elseif sum(regexp(lower(nam),'segment 52'))
        colTab(ID,:) = [0 0 1]+ rand(1,3)*.2;
    else
        colTab(ID,:) = rand(1,3)*.2+.01;
    end
    end
    
end

colTab(colTab>1) = 1;




%% Paint Colors
col = zeros(size(allSum,1), size(allSum,2), 3);

for i = 1:size(allSum,3)
    disp(sprintf('drawing %d of %d',i,size(allSum,3)))
    grabTemp = allSum(:,:,i);
    if sum(grabTemp(:))
        scaleUp = min(40,200/max(grabTemp(:)));
        scaleUp = max(scaleUp,.1);
        grabTemp = (grabTemp>0)*50 + (grabTemp* scaleUp);
      % grabTemp = grabTemp.^1.3;
        for c = 1:3
            col(:,:,c) = col(:,:,c) + grabTemp * colTab(i,c);
        end
        pause(.01)
    end
end

image(uint8(col.^1))
imwrite(uint8(col),[TPN 'sumObjsY.tif'],'Compression','none')









