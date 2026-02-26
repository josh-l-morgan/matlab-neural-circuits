function[pos] = dsAnchors(pos,obI,order)

if ~exist('order','var'); order = [2 1 3]; end
if ~exist('obI','var'); 
    dsRes = [0.2 0.2 0.2]; 
    res = [6 4 30];
else
    res = obI.em.res; 
    dsRes = obI.em.dsRes
end

anchors = pos;
dSamp =  (res .* [4 4 1])./1000./dsRes;
anchors(:,1) = anchors(:,1)/dSamp(1);
anchors(:,2) = anchors(:,2)/dSamp(2);
anchors(:,3) = anchors(:,3)/dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
pos(:,order) = anchors;