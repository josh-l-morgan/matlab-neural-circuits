function newFV=trimFV(fv,holeCenters,holeRadii)
badVertMat=zeros(length(fv.vertices(:,1)),1);
allVertRA=[1:length(fv.vertices(:,1))]';
faceLocs=zeros(length(fv.faces(:,1)),3);
for i=1:size(holeCenters,1)
    curCent=holeCenters(i,:);
    curRad=holeRadii(i,:);
    for curFace=1:length(fv.faces(:,1))
        faceLocs(curFace,:)=mean(fv.vertices(fv.faces(curFace,:),:),1);
    end
    % facePosShp=reshape(facePos,length(fv.faces(:,1)),9);
    %dists=pdist2(curCent,faceLocs);
    %badFaceList{i}=find(dists<curRad);
    % ellipse formula ((x-xo)/a)^2 + ((y-yo)/b)^2 +((z-zo)/c)^2 < 1
    badFaceList{i}=find(((faceLocs(:,1)-curCent(1))/curRad(1)).^2 + ...
        ((faceLocs(:,2)-curCent(2))/curRad(2)).^2 + ...
        ((faceLocs(:,3)-curCent(3))/curRad(3)).^2 <1);
        
end
allBadFace=vertcat(badFaceList{:});
allGoodFace=true(size(fv.faces,1),1);
allGoodFace(allBadFace)=0;
%allBadVert=vertcat(badVertAll{:});
%allBadFace=vertcat(badFaceAll{:});
%allGoodVert=allVertRA(~badVertMat);
%allGoodFace=fv.faces(setdiff(1:length(fv.faces(:,1)),allBadFace),:);
%for j=1:length(allGoodVert)
%    curVert=allGoodVert(j);
%    if curVert~=j
%    allGoodFace(find(allGoodFace==curVert))=j;
%    end
%end
%newFV.vertices=fv.vertices(allGoodVert,:);
newFV.vertices=fv.vertices;
newFV.faces=fv.faces(allGoodFace,:);

end