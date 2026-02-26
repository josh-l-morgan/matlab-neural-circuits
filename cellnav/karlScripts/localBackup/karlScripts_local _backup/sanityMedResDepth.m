
f1=figure();
hold on;
f2=figure();
hold on;
cidList=[1110,1110];
fvDirList={highFV,medFV};
colMap=colormap(turbo);
for i=1:length(cidList)
    %curFVcell=getFV(cidList(i),fvDirList{i});
    curFV=loadFV([fvDirList{i} num2str(cidList(i)) '.mat']);
    %curFV=curFVcell{1};
    %curHistDat=showHist(f2,{cidList(i)},[],fvDirList{i},'ind');
    [g,h,curDepths]=getIPLdepth(curFV.vertices(:,1),curFV.vertices(:,2), ...
        curFV.vertices(:,3),[],[]);
    %colDat=curHistDat{1,2};
    %colDat=colDat{1};
    colDat=curDepths;
    colScld=(colDat-0.3).*3;
    colScld(colScld<0.01)=0.01;
    colScld(colScld>1)=1;
    colInds=ceil(colScld.*254);
    colTrips=colMap(colInds,:);
    figure(f1);
    patches{i}=patch(curFV);
    %axis equal;
    patches{i}.EdgeColor='interp';
    patches{i}.FaceVertexCData=colTrips;
    patches{i}.FaceAlpha=0.5;
    
end



d = datacursormode;
xyz = getfield(getCursorInfo(d),'Position');
highResCoords = xyz([3 2 1]).*[250 250 25];
uint16(highResCoords)

medResCoords=highResCoords;
fvDir=fvDirList{2};
dim1 = [2 3];
dim2 = [3 1];
if exist([fvDir 'volTransform.mat'])
    load([fvDir 'volTransform.mat']);
    transType = volTransform.type;
    switch transType
        case 'YX XZ'
            medResCoords(:,dim1) = transformPointsForward(volTransform.tform1,highResCoords(:,dim1));
            medResCoords(:,dim2) = transformPointsForward(volTransform.tform2,highResCoords(:,dim2));
    end
    %scatter3(fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),'.','r')
    
end
uint16(medResCoords)
