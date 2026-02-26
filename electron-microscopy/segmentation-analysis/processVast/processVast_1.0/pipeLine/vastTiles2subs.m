function[oSubs] = vastTiles2subs(SPN);

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
    
    try
    r(i) = str2num(nam(Yp(end)+2:Xp(end)-1));
    catch err
        r(i) = 1;
    end
    try
    c(i) = str2num(nam(Xp(end)+2:dotp(end)-1));
    catch err
        c(i) = 1;
    end
    try
    s(i) = str2num(nam(sp(end)+2:Yp(end)-1))+1;
    catch err
        s(i) = 1;
    end
end


%% find locations

imageInfo = imfinfo([SPN inams{1}]);
width = imageInfo.Width;
height = imageInfo.Height;
offset = [(r-1)*height; (c-1)*width]';

bitDepth = imageInfo.BitDepth;
maxDim = [max(offset,[],1) + [height width] max(s)];

%% Read in
%%Create oSubs , in which oSubs{i} lists all the pixels y x z for object number i

clear oSubs
colormap colorcube(256)
bufOb = 100; %factor to add to object holder
objSize = zeros(10000,1); %start size for all objects;

for i = 1:length(inams)
    
    if ~mod(i,1000)
        disp(sprintf('Grabing objects from tile %d of %d',i,length(inams)))
    end
    
    if bitDepth == 24
        Iraw = double(imread([SPN inams{i}]));
        I = Iraw(:,:,1) * 256^2 + Iraw(:,:,2) * 256  + Iraw(:,:,3);
    else
        I = double(imread([SPN inams{i}]));
    end
    
    inds = find(I>0);
    [y x ] = ind2sub(size(I),inds);
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
           oSize =  size(oSubs{oId});
        catch err
            oSubs{oId} = zeros(bufOb,3);
            oSize = size(oSubs{oId});
        end
        
        if objSize(oId)>oSize(1)
            oSubs{oId}(oSize(1)+addSize*bufOb,1:3) = 0;
        end
            
        oSubs{oId}(objSize(oId)- addSize +1:objSize(oId),1:3) = objTileSubs;
    end
end

for i = 1:length(oSubs) %% clip off remaining zeros
    if objSize(i)
        oSubs{i} = oSubs{i}(1:objSize(i),1:3);
    end
end


