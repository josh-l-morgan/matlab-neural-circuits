b2three=find(curTis.syn.edges(:,1)==13&alPre{1}'==7);

figA=figure();
morf=compareMorph(gca,13,fvDir);
sins=plotSyns(b2three,curTis,[]);

figB=figure();
axB=gca;
morf=compareMorph(gca,1022,fvDir);
%sins=plotSyns(b2three,curTis,[]);
morf.patches.FaceAlpha=.4;
morf.patches.FaceColor=[0 1 1];
axB.DataAspectRatio=[.1 1 1]
set(figB,'Color',[1 1 1]);

figB=figure();
axB=gca;
morf=compareMorph(gca,1074,fvDir);
%sins=plotSyns(b2three,curTis,[]);
morf.patches.FaceAlpha=.4;
morf.patches.FaceColor=[1 0 1];
axB.DataAspectRatio=[.1 1 1]
set(figB,'Color',[1 1 1]);



for i = 1:360
   camorbit(1,0,'data',[1 0 0])
   drawnow
end

framCount=360;
%mov(framCount)= struct('cdata',[],'colormap',[]);
u=1;
while u<361
    %view(u,0);
    %drawnow
    camorbit(1,0,'data',[1 0 0])
    drawnow
    %filname=sprintf('%smov\\allVG3_rotate_%03d.png',localDir,u);
    %exportgraphics(gcf,filname,'Resolution',150);
    frameIt=uint16(u);
    mov2(frameIt)=getframe(gcf);
    u=u+1
end
movNam=sprintf('%smov\\bc3aRot.avi',localDir);
v = VideoWriter(movNam,'MPEG-4');
open(v);
writeVideo(v,mov2)
close(v);




%%
figure();
colDat=cool();
scatter(1:256,repmat(1,[1 256]),100,colDat,'o','filled')


postVtype=alPost{1}(allVout);
postVsubtype=alPost{4}(allVout);

postVrgcSubtype=postVsubtype(postVtype==1);
rgcTargTab=tabulate(postVrgcSubtype);
rgcTargTab=sortrows(rgcTargTab,2,'descend');
rgcTargTab(3,:)=[];
rgcTargTab(1,:)=[];
rgcTargTab(8,:)=[];
rgcTargTab(13,:)=[];

pieCol=colDat([128 10 245 134 140 30 20 110 250 135 105 50 80 25 170 160 140 230 200 40 70 144],:); 
%pieCol([1 2 3],:)=[[1 1 1];[1 1 1];[1 1 1]];
names=rgcTargTab(:,1);
%blanks=repmat({' '},[1 11])';
%names=vertcat(names{:},blanks{:});
names(10:21)={' '};
vals=cell2mat( rgcTargTab(:,2));
figure();
newOrder=[2 6 7 12 13 1 4 5 8 11 3 9 10 14 15:21];
piA=pie(vals(newOrder),names(newOrder));
colormap(pieCol(newOrder,:))
ax=gca;
set(piA(2:2:end),'FontSize',16)
set(gcf,'Color',[1 1 1])
'total syns to IDd RGC'
sum(cell2mat(rgcTargTab(:,2)))
'total # of IDd RGC'
v2r=find(alPre{1}'==8&alPre{3}'==1&alPost{1}'==1&ismember(alPost{4}',rgcTargTab(:,1)));
v2allR=find(alPre{1}'==8&alPre{3}'==1&alPost{1}'==1);
postRGCcids=curTis.syn.edges(v2r,1);
uniquePostRGCcids=unique(postRGCcids);
length(uniquePostRGCcids)
'syns to unidentified RGC'
length(v2allR)-length(v2r)
'IDd syn ratio'
length(v2r)/length(v2allR)


v26=getClassConn({'amc'},{'vgc'},{'rgc'},{'6sw'},curTis);
pltV26=plotSyns(v26,curTis,[]);
pltV26.MarkerFaceAlpha=1;
pltV26.MarkerFaceColor=[0 0 0]



%% histogram with depths
bpcCidsOFF=type2cid({'bpc','bpc','bpc'},{'bc5o','bc5i','bc5t'},curTis);
bpcCidsON=type2cid({'bpc','bpc','bpc'},{'bc3a','bc3b','bc4'},curTis);
bpcCidsOFF=bpcCidsOFF{:};
bpcCidsON=bpcCidsON{:};
inLoxON=curTis.syn.pos(find(ismember(curTis.syn.edges(:,2),bpcCidsON) ...
    &ismember(curTis.syn.edges(:,1),allVG3cids)),:);

inLoxOFF=curTis.syn.pos(find(ismember(curTis.syn.edges(:,2),bpcCidsOFF) ...
    &ismember(curTis.syn.edges(:,1),allVG3cids)),:);

[~,~,depsON]=getIPLdepth(inLoxON(:,3),inLoxON(:,1),inLoxON(:,2),[],[]);
[~,~,depsOFF]=getIPLdepth(inLoxOFF(:,3),inLoxOFF(:,1),inLoxOFF(:,2),[],[]);

histDatON=histcounts(depsON,[0:0.02:1]);
histDatOFF=histcounts(depsOFF,[0:0.02:1]);

figure();
arON=area([.01:.02:.99],histDatON);
hold on
arOFF=area([.01:.02:.99],histDatOFF);
arON.FaceAlpha=0.3;
arON.FaceColor=[0 1 1];
arOFF.FaceAlpha=0.3;
arOFF.FaceColor=[1 0 1];



%%

testDat=load('Y:\karlsRetina\EmilyImages\CorrectRegistratedFiles\Ai148_129SVG3_Translation_122618_2003.mat')
avgImg=mean(testDat.I,3);
testLoc=[135,15];
corrDat=reshape(testDat.I,[256*32 1600]);
goodPix=sum(corrDat,2);
pixCutoff=5000;
pixCut=find(goodPix>pixCutoff);
compList=nchoosek(1:8192,2);
corrMat=zeros(size(corrDat,1));
for j=1:size(compList,1)
    if sum(ismember(compList(j,:),pixCut))==2
    curCorr=corrcoef(corrDat(compList(j,1),:),corrDat(compList(j,2),:));
    corrMat(compList(j,1),compList(j,2))=curCorr(1,2);
    end
end

testPix=6876;
corrImg=corrMat(:,testPix)+corrMat(testPix,:)';
testImg=reshape(corrImg,32,256);
figure(); image(testImg.*350);
figure(); histogram(testImg);