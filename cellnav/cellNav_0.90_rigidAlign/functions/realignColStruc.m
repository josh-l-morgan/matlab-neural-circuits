function[obI] = realignColStruc(obI);

global glob


load("MPN.mat");

if exist([MPN 'shiftZ.mat'],'file')
    
    disp('shifting planes')
    load([MPN 'shiftZ.mat']);
    anchors = double(obI.colStruc.anchors);
    useAnchors = anchors(:,3)>0;
    pts = anchors(useAnchors,[2 1 3]);
   
    if ~isempty(pts)
         a2v = obI.em.vRes(1,:) ./ obI.em.res;
         ptsVS = pts ./ repmat(a2v,[size(pts,1) 1]);
   
        newPtsVS= rigidPtsStack(ptsVS,shiftZ.As);
        newPts = newPtsVS .* repmat(a2v,[size(pts,1) 1]);
        anchors(useAnchors,[2 1 3]) = newPts;
        obI.colStruc.anchors = anchors;
    end


end
