
%% part A
allPreTypes=cid2type(allSynEdges(:,2),curTis);
allPostTypes=cid2type(allSynEdges(:,1),curTis);

bpcInputTypes=zeros(length(vgcCidList),length(curTis.cells.type.subTypeNames{7}));
for i=1:length(vgcCidList)
    curV=vgcCidList(i);
    bpcInIDs=find(allPreTypes{1}==7 & allSynEdges(:,1)'==curV);
    bpcInHistDat=histcounts(allPreTypes{3}(bpcInIDs),1:length(curTis.cells.type.subTypeNames{7})+1);
    bpcInputTypes(i,:)=bpcInHistDat;
end

fA=figure();
hold on
title('BPC subtype inputs to VG3');
image(bpcInputTypes*5);
yticks([1:length(vgcCidList)]);
yticklabels({'2','3','4','5','13','14'});
xticks([1:length(curTis.cells.type.subTypeNames{7})]);
xticklabels(curTis.cells.type.subTypeNames{7});
ylabel('VGC cell #');

%% part B

allBpcInIDs=find(allPreTypes{1}'==7 & ismember(allSynEdges(:,1),vgcCidList));
bpcOFFtypes=[3 4 5];
bpcONtypes=[6 7 8 9];
bpcBothtype=[bpcOFFtypes bpcONtypes];
bpcLocPlotDat=cell(1,length(bpcBothtype));
for i=1:length(bpcBothtype)
    curType=bpcBothtype(i);
    bpcLocPlotDat{i}=allSynPos((allPreTypes{1}'==7 & ismember(allSynEdges(:,1),vgcCidList) & allPreTypes{3}'==curType),:);
end


bpcMat=zeros(2000,2000,length(bpcLocPlotDat));
emptySlice=zeros(2000,2000);
for i=1:length(bpcLocPlotDat)
    curSlice=emptySlice;
    curLocs=bpcLocPlotDat{i};
    curPix=round(curLocs(:,1:2).*[10 10]);
    curPix=curPix((curPix(:,1)>0 & curPix(:,2)>0),:);
    %offMat(curPix(:,1),curPix(:,2))=1;
    curSlice(sub2ind(size(curSlice),curPix(:,1),curPix(:,2)))=1;
    bpcMat(:,:,i)=curSlice;
end

offMat=zeros(2000,2000,3);
offMat(:,:,[1 3])=offMat(:,:,[1 3])+bpcMat(:,:,1);
offMat(:,:,[2 3])=offMat(:,:,[2 3])+bpcMat(:,:,2);
offMat(:,:,[1 2])=offMat(:,:,[1 2])+bpcMat(:,:,3);

offMatBlr=imgaussfilt(offMat,15);
offMatBlr=offMatBlr/max(offMatBlr(:));
fB=figure();
%hold on;
%title('Where are the BPC input synapses in XY?');
tL=tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
%subplot(1,2,1);
nexttile
imshow(offMatBlr);
xlim([800 2000]);
ylim([800 2000]);
text(1700,1000,'bc3a','Color','magenta','FontSize',16);
text(1700,1050,'bc3b','Color','cyan','FontSize',16);
text(1700,1100,'bc4','Color','yellow','FontSize',16);

% same for the ON bpcs
onMat=zeros(2000,2000,3);
onMat(:,:,1)=bpcMat(:,:,5).*.85+bpcMat(:,:,6).*.929+bpcMat(:,:,7)*.494;
onMat(:,:,2)=bpcMat(:,:,4).*.447+bpcMat(:,:,5).*.325+bpcMat(:,:,6).*.694+bpcMat(:,:,7)*.184;
onMat(:,:,3)=bpcMat(:,:,4).*.741+bpcMat(:,:,5).*.098+bpcMat(:,:,6).*.125+bpcMat(:,:,7)*.556;

onMatBlr=imgaussfilt(onMat,15);
onMatBlr=onMatBlr/(max(onMatBlr(:)*.9));
%subplot(1,2,2);
nexttile
imshow(onMatBlr);
xlim([800 2000]);
ylim([800 2000]);
text(1700,1000,'bc5i','Color','#0072BD','FontSize',16);
text(1700,1050,'bc5o','Color','#D95319','FontSize',16);
text(1700,1100,'bc5t','Color','#EDB120','FontSize',16);
text(1700,1150,'bc6','Color','#7E2F8E','FontSize',16);

%% part C




