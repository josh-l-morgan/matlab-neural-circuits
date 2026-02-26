function graphObj = compareMorph(ax,cidList,fvDir,labLoc)

if ~exists('labLoc','var')
    labLoc='INL';
end

if isempty(fvDir)
    fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
end
figNam=gcf;
if isempty(ax)
    ax=gca;
end
%if isempty

colPal=[[0.537254902	0.192156863	0.937254902]; ...
    [0.949019608	0.792156863	0.098039216]; ...
    [1	0	0.741176471]; ...
    [0	0.341176471	0.91372549]; ...
    [0.529411765	0.91372549	0.066666667]; ...
    [0.882352941	0.094117647	0.270588235]];
colPal=colPal(randperm(size(colPal, 1)), :);
colPal2=colPal(:,[3 1 2]);
colPal3=colPal(:,[2 3 1]);
colPal4=colPal(:,[3 2 1]);
colPal5=colPal(:,[2 1 3]);
colPal6=colPal(:,[3 1 2]);
colPal=vertcat(colPal,colPal2,colPal3,colPal4,colPal5,colPal6);

figure(figNam);
%clf
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
    if mean(labLoc=='GCL')==1
    labLocPts=curFV.vertices(find(curFV.vertices(:,1)==min(curFV.vertices(:,1))),:);
    labLocPts=labLocPts-[2 0 0];
    elseif mean(labLoc=='INL')==1
    labLocPts=curFV.vertices(find(curFV.vertices(:,1)==max(curFV.vertices(:,1))),:);
    labLocPts=labLocPts+[2 0 0];
    else
    labLocPts=mean(curFV.vertices);
    end
    labLocC=mean(labLocPts,1);
    curT=text(labLocC(1),labLocC(2),labLocC(3),num2str(cidList(i)));
    curP=patch(curFV);
    curP.EdgeColor='none';
    %curP.EdgeAlpha=0;
    %curP.EdgeColor=colPal(i,:);
    curP.FaceColor=colPal(i,:);
    curP.FaceAlpha=.25;
    patches(i)=curP;
    texts(i)=curT;
end
%axes('Units', 'normalized', 'Position', [0 0 1 1])
    %ax=gca;
view(90,0)
%axis equal
%axis vis3d
%axis vis3d equal
xlim([min([squeeze(bbox(2,1,:));ax.XLim(1)]) max([squeeze(bbox(1,1,:));ax.XLim(2)])]);
xlim([min([squeeze(bbox(2,2,:));ax.YLim(1)]) max([squeeze(bbox(1,2,:));ax.YLim(2)])]);
xlim([min([squeeze(bbox(2,3,:));ax.ZLim(1)]) max([squeeze(bbox(1,3,:));ax.ZLim(2)])]);
%ylim([min([bbox(2,2,:),ax.YLim(1)]) max([bbox(1,2,:),ax.YLim(2)])]);
%zlim([min([bbox(2,3,:),ax.ZLim(1)]) max([bbox(1,3,:),ax.ZLim(2)])]);
%axis vis3d
axis vis3d equal
ax.Clipping='off';
ax.Visible='off';
set(gcf,'color',[.9 .9 .9]);

graphObj.figNam=figNam;
graphObj.patches=patches;
graphObj.axes=ax;
graphObj.texts=texts;
end