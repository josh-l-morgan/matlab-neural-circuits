
figure();
curSeed=2008; 
curPart=1116;
%figure();
hold on

seedVxls=bpcVxls{bpcCids==curSeed};
seedPoly=bpcBoundPolys{bpcCids==curSeed};

compIdxs=find(offDat(:,5)==curSeed|offDat(:,6)==curSeed);
curDat=offDat(compIdxs,:);
[x,srtIdx]=sort(curDat(:,2),'ascend');
curDatSrtd=curDat(srtIdx,:);
testGrp=curDatSrtd(curDatSrtd(:,2)<300,5:6);
testGrp=testGrp(:);
testGrp=testGrp(testGrp~=curSeed);
testGrp=testGrp(~ismember(testGrp,gg));
testGrp=testGrp(~ismember(testGrp,bg));
y=1;
r='';

clf
scatter3(seedVxls(:,1),seedVxls(:,2),seedVxls(:,3),1,'k.');
hold on
partVxls=bpcVxls{bpcCids==curPart};
scatter3(partVxls(:,1),partVxls(:,2),partVxls(:,3),1,'c.');
partPoly=bpcBoundPolys{bpcCids==curPart};
curVerts=partPoly.Vertices;
p2=patch(curVerts(:,1),curVerts(:,2),rand([1 3]));
p2.FaceColor=[1 1 .8];

curVerts=seedPoly.Vertices;
p1=patch(curVerts(:,1),curVerts(:,2),rand([1 3]));
%p1.FaceColor=(rand([1 3])+.5)./1.5;
p1.FaceColor=[.6 1 1];

interPoly=intersect(seedPoly,partPoly);
%interBound=boundary(interPoly);
curVerts=interPoly.Vertices;
goodVerts=sum(isnan(curVerts),2)==0;
if ~isempty(goodVerts)
    upatch=patch(curVerts(goodVerts,1),curVerts(goodVerts,2),rand([1 3]));
    upatch.FaceColor=[1 0 1];
end
view(0,90);
drawnow