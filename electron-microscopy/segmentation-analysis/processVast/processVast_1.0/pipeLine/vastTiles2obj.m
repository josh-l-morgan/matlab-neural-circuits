function[subs] = vastTiles2obj(SPN);


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
%%Create subs , in which subs{i} lists all the pixels y x z for object number i

clear subs
colormap colorcube(256)
bufOb = 100; %factor to add to object holder
objSize = zeros(10000,1); %start size for all objects;
for i = 1:length(inams)
    if ~mod(i,1000)
    disp(sprintf('Grabing objects from tile %d of %d',i,length(inams)))
    end
    
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
           oSize =  size(subs{oId});
        catch err
            subs{oId} = zeros(bufOb,3);
            oSize = size(subs{oId});
        end
        
        if objSize(oId)>oSize(1)
            subs{oId}(oSize(1)+addSize*bufOb,1:3) = 0;
        end
            
        subs{oId}(objSize(oId)- addSize +1:objSize(oId),1:3) = objTileSubs;
    end
end

for i = 1:length(subs) %% clip off remaining zeros
    if objSize(i)
        subs{i} = subs{i}(1:objSize(i),1:3);
    end
end


