function[allAnchors useSynID] = getSynAnchors(obI,prePop,postPop,synType)


if ~exist('synType','var')
    synType = 1:10000;
end
typeVec = zeros(10000,1);
typeVec(synType+1) = 1;

anchors = double(obI.colStruc.anchors);

dSamp =  (obI.em.res )./1000./obI.em.dsRes;
%dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
%dSamp = dSamp ./ [4 4 2];

anchors(:,1) = anchors(:,1)*dSamp(1);
anchors(:,2) = anchors(:,2)*dSamp(2);
anchors(:,3) = anchors(:,3)*dSamp(3);
%anchors = round(anchors);
anchors(anchors<1) = 1;
if ~isempty(obI.nameProps.edges)
    synapses = obI.nameProps.edges(:,1:3);
    
    %
    %         dim = viewProps.dim;
    %         if dim == 1
    %             dims = [3 2];
    %         elseif dim == 2
    %             dims = [3 1];
    %         elseif dim == 3
    %             dims = [1 2];
    %         end
    
    
    
    usePre = [];
    for i = 1:length(prePop)
        usePre = cat(1,usePre,find(synapses(:,2)== prePop(i)));
    end
    
    usePost = [];
    for i = 1:length(postPop)
        usePost = cat(1,usePost,find(synapses(:,1)== postPop(i)));
    end
    
    foundSyn = intersect(usePre,usePost);
    
    if size(obI.nameProps.edges,2) >3
        foundType = obI.nameProps.edges(foundSyn,4);
        goodFind = typeVec(foundType+1);
        foundSyn = foundSyn(goodFind>0);
    end
    
    useSynID = synapses(foundSyn,3);
    allAnchors = anchors(useSynID,:);
else
    
    useSynID = [];
    allAnchors = [];
    
end







