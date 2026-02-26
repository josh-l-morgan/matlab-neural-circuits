function[sub] = stereoPoints_PS(viewProps);

%%

dim = viewProps.dim;
fsize = ceil(viewProps.fsize);
viewWindow = viewProps.viewWindow;


if isfield(viewProps,'resamp')
    resamp = viewProps.resamp;
    viewWindow = viewWindow .* [resamp;resamp];
end

rotDim = 2;

if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end
changeDim = dims(rotDim);



viewSize = viewWindow(2,:) - viewWindow(1,:);
center = viewSize/2;
maxDist = sqrt(center(changeDim).^2 + center(dim).^2);
fsize = ceil(viewSize);
shiftRotDim = maxDist - fsize(changeDim)/2;
% fsize(changeDim) = ceil(maxDist*2);
center(changeDim) = center(changeDim)+shiftRotDim/2;

% anchorScale = [.0184 0.016 0.030];
% voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];
% subScale = voxelScale/max(voxelScale(:));



if ~isfield(viewProps,'degRot')
    degRot = 0;
else
    degRot = viewProps.degRot;
end

if ~isfield(viewProps,'minInt')
    minInt = 20;
else
    minInt = viewProps.minInt;
end

if ~isfield(viewProps,'perspective')
    perspective = 1;
else
    perspective = viewProps.perspective;
end

if isfield(viewProps,'mask')
    mask = viewProps.mask;
    
    mask(:,1) = mask(:,1) - viewWindow(1,1)+1;
    mask(:,2) = mask(:,2) - viewWindow(1,2)+1;
    mask(:,3) = mask(:,3) - viewWindow(1,3)+1;
    
    
    
    if isfield(viewProps,'useMask')
        useMask = viewProps.useMask;
    else
        useMask = ones(length(viewProps.cellId),1);
    end
else
    useMask = zeros(length(viewProps.cellId),1);
end

if ~isfield(viewProps,'contrastFactor')
    contrastFactor = 5;
else
    contrastFactor = viewProps.contrastFactor;
end

if isfield(viewProps,'sumScaleFactor')
    sumScaleFactor = viewProps.sumScaleFactor;
else
    sumScaleFactor = 1; %3
end

if isfield(viewProps,'maxScaleFactor')
    maxScaleFactor = viewProps.maxScaleFactor;
else
    maxScaleFactor = .3;
end

if viewProps.viewWindow == 0;
    viewWindow = [1 1 1; fsize];
end



%%
%fsize = [1700 1700 1300];
I = zeros(fsize(dims));
newHeight = I;
maxHeight = I;
Ic = cat(3,I,I,I);
IcSum = Ic;



sub = viewProps.points;

if isfield(viewProps,'resamp')
    sub = scaleSubsUnique(sub,resamp);
end

if ~isempty(sub)
    
    
    %%Apply window
    useSub = (sub(:,1)>=viewWindow(1,1)) & (sub(:,1)<=viewWindow(2,1)) & (sub(:,2)>=viewWindow(1,2)) & (sub(:,2)<=viewWindow(2,2)) & (sub(:,3)>=viewWindow(1,3)) & (sub(:,3)<=viewWindow(2,3));
    sub = sub(useSub,:);
    sub(:,1) = sub(:,1) - viewWindow(1,1)+1;
    sub(:,2) = sub(:,2) - viewWindow(1,2)+1;
    sub(:,3) = sub(:,3) - viewWindow(1,3)+1;
    
    if useMask(1)>0
        %%Apply mask
        maskInd = sub2ind(fsize,mask(:,1),mask(:,2),mask(:,3));
        subInd = sub2ind(fsize,sub(:,1),sub(:,2),sub(:,3));
        [sameInd idxA idxB] = intersect(subInd,maskInd);
        sub = sub(idxA,:);
    end
    
    
    sub(:,changeDim) = sub(:,changeDim) + shiftRotDim;
    
    if 1
        %%rotate
        otherDim = setdiff([1 2 3],dims);
        oldSub = sub;
        %                 center = [viewWindow(2,dims(2)) - viewWindow(1,dims(2)) ...
        %                     viewWindow(2,otherDim) - viewWindow(1,otherDim)]/2;
        
        dif1  = oldSub(:,dims(rotDim)) - center(dims(rotDim));
        dif2 = oldSub(:,otherDim) - center(otherDim);
        dists = sqrt(dif1.^2 + dif2.^2);
        
        rads = atan2(dif1,dif2);
        %rads(isnan(rads)) = 0;
        rads = rads - (pi/180)*degRot;
        
        new2 = cos(rads) .* dists;
        new1 = sin(rads) .* dists;
        
        %                     new2 = sin(rads) .* dists;
        %                     new1 = cos(rads) .* dists;
        
        sub(:,dims(rotDim)) = new1+center(dims(rotDim));
        sub(:,otherDim) = new2 + center(otherDim);
    end
    
    
    if perspective
        
        %%add perspective
        %distTo = 1;
        
        windowSize = (viewWindow(2,:) - viewWindow(1,:));
        reCenter = windowSize/2;
        distTo = max(windowSize);
        distTo = 1;

        sub1 = sub(:,dims(1))-reCenter(dims(1));
        sub2 = sub(:,dims(2))-reCenter(dims(2));
        
        
        zDepth = windowSize(dim) -(sub(:,dim)); %distance to from front
        s2 = atan2(reCenter(dim(1)),windowSize(dim)*distTo); %angle to top at distance
        bigF = tan(s2).*(windowSize(dim)*distTo + zDepth); % size of field at angle
        
        smallF = tan(s2).*(windowSize(dim)*distTo);
        zScale = smallF./bigF ;
        sub(:,dims(1)) = sub1.* zScale + reCenter(dims(1));
        sub(:,dims(2)) = sub2.* zScale + reCenter(dims(2));
        
        sub = round(sub);
    end
    
    
    for c = 1:3
        badSub = sub(:,c)<1;
        sub(badSub,c) = 1;
        badSub = sub(:,c)>fsize(c);
        sub(badSub,c) = fsize(c);
    end
    sub = round(sub);
    
    inds = sub2ind(fsize(dims),sub(:,dims(1)),sub(:,dims(2)));
    uinds = unique((inds));
    
    if length(uinds)>1
        hinds = hist(inds,uinds);
    else
        hinds = length(inds);
    end
    
    try I(uinds) = I(uinds) + hinds';
    catch err
        'something went horribly wrong. -stereo cells and more 234'
    end
    newHeight(inds) = max(newHeight(inds),sub(:,dim));
    
    
    I = I * contrastFactor;
    I((I>0)) = I(I>0) + minInt;
    %I(I>255) = 255;
    I((I>0) & (I<minInt)) = minInt;
    
    
    useI = newHeight>maxHeight;
    maxHeight(useI) = newHeight(useI);
    
    for c = 1:3
        Itemp = Ic(:,:,c);
        Itemp(useI) = I(useI)*1;
        Ic(:,:,c) = Itemp;
    end
    
    IcSum(:,:,1) = IcSum(:,:,1) + I;
    IcSum(:,:,2) = IcSum(:,:,2) + I;
    IcSum(:,:,3) = IcSum(:,:,3) + I;
    
    
    %image(uint8(IcSum)*300),    pause(.001)
    I = I * 0;
    
end


I_topSum = Ic * maxScaleFactor + IcSum * sumScaleFactor;
%image(uint8(I_topSum*1000))
%Ic = uint8(Ic);


end