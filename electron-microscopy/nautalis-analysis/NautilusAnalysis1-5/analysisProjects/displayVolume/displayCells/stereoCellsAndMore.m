function[I_topSum,newHeight] = showCellsAndMore(viewProps);

%%
cellId = viewProps.cellId;
obI = viewProps.obI;
dsObj = viewProps.dsObj;
col = viewProps.col;
dim = viewProps.dim;
fsize = viewProps.fsize;
viewWindow = viewProps.viewWindow;

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

if isfield(viewProps,'mask')
    useMask = 1;
    mask = viewProps.mask;
else
    useMask = 0;
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

fsize = viewWindow(2,:) - viewWindow(1,:) + 1;

%%

if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end


%%
%fsize = [1700 1700 1300];
I = zeros(fsize(dims));
newHeight = I;
maxHeight = I;
Ic = cat(3,I,I,I);
IcSum = Ic;

if isempty(cellId)
    uCellId = 1:length(dsObj);
    
else
    uCellId = cellId;
end

if length(uCellId) ~= size(col,1)
    disp('Bad color size')
%     colMap = hsv(256);
%     col = colMap(ceil((1:length(uCellId))*256/length(uCellId)),:);
%     col = col(randperm(size(col,1)),:);
end

for i = 1:length(uCellId)
    
    if strcmp(class(uCellId),'cell')
        term = uCellId{i};
        if strcmp(class(term),'char')
            for n = 1:length(obI.nameProps.names)
                isTerm(n) = sum(regexp(lower(obI.nameProps.names{n}),lower(term)));
            end
            obTarg = find(isTerm);
        else
            targ = find(obI.cell.name==term);
            if isempty(targ)
                obTarg = [];
                disp(sprintf('Could not find %d',term))
            else
            obTarg = obI.cell.obIDs{targ};
            end
        end
    else
        obTarg = uCellId(i);
    end
    
    for o = 1:length(obTarg)
        if obTarg(o)<=length(dsObj)
            sub = double(dsObj(obTarg(o)).subs);
            %sub = scaleSubs(sub,subScale);
            
            
            if ~isempty(sub)
                
                if useMask
                %%Apply mask
                maskInd = sub2ind(fsize,mask(:,1),mask(:,2),mask(:,3));
                subInd = sub2ind(fsize,sub(:,1),sub(:,2),sub(:,3));
                [sameInd idxA idxB] = intersect(subInd,maskInd);
                sub = sub(idxA,:);
                end
                
                %%Apply window
                useSub = (sub(:,1)>=viewWindow(1,1)) & (sub(:,1)<=viewWindow(2,1)) & (sub(:,2)>=viewWindow(1,2)) & (sub(:,2)<=viewWindow(2,2)) & (sub(:,3)>=viewWindow(1,3)) & (sub(:,3)<=viewWindow(2,3));
                sub = sub(useSub,:);
                sub(:,1) = sub(:,1) - viewWindow(1,1)+1;
                sub(:,2) = sub(:,2) - viewWindow(1,2)+1;
                sub(:,3) = sub(:,3) - viewWindow(1,3)+1;
                
                if 1
                    %%rotate
                    otherDim = setdiff([1 2 3],dims);
                    oldSub = sub;
                    %                 center = [viewWindow(2,dims(2)) - viewWindow(1,dims(2)) ...
                    %                     viewWindow(2,otherDim) - viewWindow(1,otherDim)]/2;
                    center = (viewWindow(2,:) - viewWindow(1,:))/2;
                    
                    dif1  = oldSub(:,dims(2)) - center(dims(2));
                    dif2 = oldSub(:,otherDim) - center(otherDim);
                    dists = sqrt(dif1.^2 + dif2.^2);
                    
                    rads = atan2(dif1,dif2);
                    %rads(isnan(rads)) = 0;
                    rads = rads - (pi/180)*degRot;
                    
                    new2 = cos(rads) .* dists;
                    new1 = sin(rads) .* dists;
                    
%                     new2 = sin(rads) .* dists;
%                     new1 = cos(rads) .* dists;
                    
                    sub(:,dims(2)) = new1+center(dims(2));
                    sub(:,otherDim) = new2 + center(otherDim);
                    
                    sub = round(sub);
                    for c = 1:3
                        badSub = sub(:,c)<1;
                        sub(badSub,c) = 1;
                        badSub = sub(:,c)>fsize(c);
                        sub(badSub,c) = fsize(c);
                    end
                    
                    
                end
                
                inds = sub2ind(fsize(dims),sub(:,dims(1)),sub(:,dims(2)));
                uinds = unique(inds);
                
                if length(uinds)>1
                    hinds = hist(inds,uinds);
                else
                    hinds = length(inds);
                end
                
                I(uinds) = I(uinds) + hinds';
                newHeight(inds) = max(newHeight(inds),sub(:,dim));
                
            end
            
        end
        
    end
    
    I = I * contrastFactor;
    I((I>0)) = I(I>0) + minInt;
    %I(I>255) = 255;
    I((I>0) & (I<minInt)) = minInt;
    
    
    useI = newHeight>maxHeight;
    maxHeight(useI) = newHeight(useI);
    
    for c = 1:3
        Itemp = Ic(:,:,c);
        Itemp(useI) = I(useI)*col(i,c);
        Ic(:,:,c) = Itemp;
    end
    
    IcSum(:,:,1) = IcSum(:,:,1) + I*col(i,1);
    IcSum(:,:,2) = IcSum(:,:,2) + I*col(i,2);
    IcSum(:,:,3) = IcSum(:,:,3) + I*col(i,3);
    
    
    image(uint8(IcSum)*3),    pause(.001)
    I = I * 0;
    
end


I_topSum = Ic * maxScaleFactor + IcSum * sumScaleFactor;
%image(uint8(I_topSum))
%Ic = uint8(Ic);


