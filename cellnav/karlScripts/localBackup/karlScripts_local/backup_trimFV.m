function newFV=trimFV(fv,holeCenters,holeRadii)
badVertMat=zeros(length(fv.vertices(:,1)),1);
allVertRA=[1:length(fv.vertices(:,1))]';
for i=1:size(holeCenters,1)
    curCent=holeCenters(i,:);
    curRad=holeRadii(i,:);
    badVert=abs(fv.vertices(:,1)-curCent(1))<curRad(1) & ...
        abs(fv.vertices(:,2)-curCent(2))<curRad(2) & ...
        abs(fv.vertices(:,3)-curCent(3))<curRad(3);
    vertRA=[1:length(fv.vertices(:,1))]';
    vertRAtrm=vertRA(~badVert);
    badFace=any(ismember(fv.faces,find(badVert)),2);
    badVertMat=max(horzcat(badVertMat,badVert),[],2);
    badVertAll{i}=find(badVert);
    badFaceAll{i}=find(badFace);
end
allBadVert=vertcat(badVertAll{:});
allBadFace=vertcat(badFaceAll{:});
allGoodVert=allVertRA(~badVertMat);
allGoodFace=fv.faces(setdiff(1:length(fv.faces(:,1)),allBadFace),:);
for j=1:length(allGoodVert)
    curVert=allGoodVert(j);
    if curVert~=j
    allGoodFace(find(allGoodFace==curVert))=j;
    end
end
newFV.vertices=fv.vertices(allGoodVert,:);
newFV.faces=allGoodFace;



% newFV.vertices=fv.vertices(~badFace,:);
   % newFV.faces=fv.faces(~badFace,:);
    %newFace=changem(fv.faces,vertRAtrm,1:length(vertRAtrm));
   % newFV.faces=newFace;
%end
end