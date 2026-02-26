
MPN = GetMyDir;
load([MPN 'dsObj.mat'])
load([MPN 'vastSubs.mat'])
load([MPN 'obI.mat'])

anchorScale = [.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];
vastScale = [1 4/4.6 1];
pixum = .0184 * 8;

for i = 1:length(dsObj)
    dsObj(i).subs = scaleSubs(double(dsObj(i).subs),voxelScale);
end

allAnchors = scaleSubs(obI.cell.anchors,anchorScale);
cellRefList = obI.cell.name;


names = obI.nameProps.names;
IDs = []; obRefs = [];
for i = 1:length(names)
    nam = names{i};
    if sum(regexp(nam,'circ '));
        fragPos = regexp(nam,'frag');
        numStr = nam(fragPos(1)+5:end);
        ID = str2num(numStr);
        IDs = [IDs ID];
        obRefs = [obRefs i];
        
    end
end

[IDs idx] = sort(IDs,'ascend');
obRefs = obRefs(idx);


%%
clear props
for i = 1:length(IDs)
   
    cSubs = vastSubs{obRefs(i)};
    cSubs = scaleSubs(cSubs,vastScale);
    for d = 1:3
        cSubs(:,d) = ceil(cSubs(:,d)-floor(min(cSubs(:,d),[],1)))+1;
    end
    imax = max(cSubs(:,[1 2]),[],1);
    iraw = zeros(imax(1:2));
    iind = sub2ind(imax,cSubs(:,1),cSubs(:,2));
    iraw(iind) = 1;
    
    %Isc = imscale(
    
    
    Ifill = regionprops(iraw>0,'FilledImage','ConvexImage','Image')
    
    I = Ifill.FilledImage;
    Icon = Ifill.ConvexImage;
    Icol = cat(3,I,Icon,Ifill.Image);
    image(uint8(Icol*300))
    if sum(Icon(:))>(2*sum(I(:)))
        I = Icon;
        image(I*100);
    end
    pause(.01)

    
    
    props(i) = regionprops(I,'Area','Orientation','Extent', ...
        'Extrema','ConvexArea','ConvexHull','Solidity',...
        'Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength');
    
end

cb2d.IDs = IDs;
cb2d.obRefs = obRefs;
cb2d.props = props;
cb2d.areaUM = (sqrt([props(:).Area]) * pixum).^2;
cb2d.longRat = [props(:).MajorAxisLength] ./ [props(:).MinorAxisLength];
cb2d.pixum = pixum; 
cb2d.orientation = [props(:).Orientation];
cb2d

save([MPN 'cb2d.mat'],'cb2d')







