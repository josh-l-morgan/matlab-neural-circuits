function[cellSubs] = obj2conBranch(obI,dsObj,cid)

cellTarg = obI.cell.name == cid;
obIds = obI.cell.obIDs{cellTarg};


indiTag = {'brh'}; %tags used to identify objects as independent branches
isIndi = [];
for t= 1:length(indiTag)
    isIndi = cat(1,isIndi,eval(sprintf('%s.%s','obI.nameProps.tag',indiTag{t})));
end
indis = find(sum(isIndi,1));

clear indiAnc indID subs indID
indiOb = intersect(obIds,indis);
mainOb = setdiff(obIds,indis);
subs = uniqueSubs(cat(1,dsObj(mainOb).subs));

while isempty(subs) & ~isempty(indiOb)
    mainOb = indiOb(1);
    indiOb = setdiff(indiOb,mainOb);
    subs = uniqueSubs(cat(1,dsObj(mainOb).subs));
    
end

if ~isempty(subs)
    
    indID = zeros(size(subs,1),1);
    [conMat conDat] = obj2conS(subs);
    subs = double(subs);
    reSamp = obI.em.res./(1000.*obI.em.dsRes);
    
    indiNum = length(indiOb);
    
    for i = 1:indiNum
        subNum = size(subs,1);
        iSub = double(uniqueSubs(dsObj(indiOb(i)).subs));
        
        if ~isempty(iSub) %% Avoid crash
            subs = cat(1,subs,iSub);
            iConMat = obj2conS(iSub);
            iConMat(iConMat>0) = iConMat(iConMat>0) + subNum;
            conMat = cat(1,conMat,iConMat);
            indID = cat(1,indID,ones(size(iSub,1),1)*i);
        end
        indiAnc(i,:) = double(obI.colStruc.anchors(indiOb(i),[2 1 3]));
        
    end
    
    if indiNum>0
        indiAnc = indiAnc .* repmat(reSamp,[size(indiAnc,1) 1]);
    end
    %
    % [subs ia ic ] = uniqueSubs(double(subs));
    % indID = indID(ia);
    %
    % conMat = obj2conS(subs);
    
    %% Get branch seeds
    for i = 1:indiNum
        isSurf = sum(conMat,2)>0;
        surfInd = find(isSurf & (indID == i));
        surfSubs = (subs(surfInd,:));
        dist = sqrt((surfSubs(:,1)-indiAnc(i,1)).^2 + (surfSubs(:,2)-indiAnc(i,2)).^2 + ...
            (surfSubs(:,3)-indiAnc(i,3)).^2);
        surfTarg = find(dist==min(dist));
        branchSeed = surfInd(surfTarg);
        
        dist = sqrt((subs(:,1)-indiAnc(i,1)).^2 + (subs(:,2)-indiAnc(i,2)).^2 + ...
            (subs(:,3)-indiAnc(i,3)).^2);
        
        minDist = min(dist(~(indID == i)));
        otherTarg = find((dist==minDist) & ~(indID == i),1);
        
        %     isZero = find(conMat(branchSeed,:)==0,1);
        %     conMat(branchSeed,isZero) = otherTarg;
        %     isZero = find(conMat(otherTarg,:)==0,1);
        %     conMat(otherTarg,isZero) = branchSeed;
    end
    
    
    %
    % %clean redundants
    % [subs ia ic] = uniqueSubs(subs);
    % indID = indID(ia);
    %
    % conMat
    
    
    
    
    
    
    cellSubs.subs = subs;
    cellSubs.brhID = indID;
    cellSubs.conMat = conMat;
    cellSubs.conDat = conDat;
else
    cellSubs = [];
    
end


























