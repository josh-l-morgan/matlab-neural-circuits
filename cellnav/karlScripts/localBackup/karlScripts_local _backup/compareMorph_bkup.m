function graphObj = compareMorph(figNam,cidList,fvDir)
if isempty(fvDir)
    fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
end
if isempty(figNam)
    figNam=figure();
end
colPal=[[0.537254902	0.192156863	0.937254902]; ...
    [0.949019608	0.792156863	0.098039216]; ...
    [1	0	0.741176471]; ...
    [0	0.341176471	0.91372549]; ...
    [0.529411765	0.91372549	0.066666667]; ...
    [0.882352941	0.094117647	0.270588235]];
colPal2=colPal(:,[3 1 2]);
colPal3=colPal(:,[2 3 1]);
colPal4=colPal(:,[3 2 1]);
colPal5=colPal(:,[2 1 3]);
colPal6=colPal(:,[3 1 2]);
colPal=vertcat(colPal,colPal2,colPal3,colPal4,colPal5,colPal6);
figure(figNam);
clf
hold on
bbox=zeros(2,3,length(cidList));
for i=1:length(cidList)
    curFV=load([fvDir num2str(cidList(i)) '.mat']);
    curFV=curFV.fv;
    if 1
        maxes=max(curFV.vertices);
        mins=min(curFV.vertices);%(curFV.vertices(:,2)>121,:));
        bbox(1,:,i)=maxes;
        bbox(2,:,i)=mins;
    end
    if 0
        badVerts=find(curFV.vertices(:,2)<121|curFV.vertices(:,3)<121);
        badFaces=find(ismember(badVerts,curFV.faces));
        vertTrack=1:size(curFV.vertices,1);
        vertTrack(badVerts)=[];
        curFV.faces(badFaces,:)=[];
        curFV.vertices(badVerts,:)=[];
        for j=1:length(vertTrack)
            curFV.faces(find(curFV.faces==vertTrack(j)))=j;
        end
    end
    vertHistDat=histcounts(curFV.vertices(:,1),[-5:10:55]);
    
    labLocPts=curFV.vertices(find(curFV.vertices(:,1)==max(curFV.vertices(:,1))),:);
    labLoc=mean(labLocPts,1);
    text(labLoc(1),labLoc(2),labLoc(3),num2str(cidList(i)));
    curP=patch(curFV);
    curP.EdgeColor='none';
    %curP.EdgeAlpha=0;
    %curP.EdgeColor=colPal(i,:);
    curP.FaceColor=colPal(i,:);
    curP.FaceAlpha=.25;
end
view(90,0)
axis equal
%axis vis3d
xlim([min(bbox(2,1,:)) max(bbox(1,1,:))]);
ylim([min(bbox(2,2,:)) max(bbox(1,2,:))]);
zlim([min(bbox(2,3,:)) max(bbox(1,3,:))]);
axis vis3d


end