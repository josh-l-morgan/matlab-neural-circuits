%% load up the regular suspects
tisRaw=load([localDir 'Volumes\Final\Analysis\fvLibrary\tis.mat']);
objRaw=load([localDir 'Volumes\Final\Merge\dsObj.mat']);
obiRaw=load([localDir 'Volumes\Final\Analysis\fvLibrary\obI.mat']);
dsObj=objRaw.dsObj;
tis=tisRaw.tis;
obI=obiRaw.obI;
clear tisRaw objRaw obiRaw;

fvDir='C:\work\localDat\Volumes\Final\Analysis\fvLibrary\';
medFVdir='C:\work\localDat\Volumes\medRes\Analysis\fvLibrary\';
addpath('C:\work\localDat\karlScripts\')
curTis=tis;

%% get the input ratios for VG3
connDat=struct();
connNames={'a2v','p2v','b2v','v2p','v2a','v2r'};
classConnMat={'amc','unk','amc','vgc'; ...
    'amc','xac','amc','vgc'; ...
    'bpc','all','amc','vgc'; ...
    'amc','vgc','amc','xac'; ...
    'amc','vgc','amc','unk'; ...
    'amc','vgc','rgc','all'; ...
    };
for i=1:length(connNames)
    curCell=getClassConn(classConnMat(i,1),classConnMat(i,2),classConnMat(i,3),classConnMat(i,4),curTis);
	curList=curCell;
    connDat(i).name=connNames(i);
    connDat(i).synList=curList;
end

%% picking up the stragglers
allah2v=getClassConn([8,8,8,8],[0,2,4,6],[8,8,8,8],[1,1,1,1],curTis);
allVG3cids=type2cid({'amc'},{'vgc'},curTis);
allVG3cids=allVG3cids{1};
alPre=cid2type(curTis.syn.edges(:,2),curTis);
alPost=cid2type(curTis.syn.edges(:,1),curTis);
amc02v=find(ismember(curTis.syn.edges(:,1),allVG3cids)&alPre{1}'==8&alPre{3}'==0);
test=alPre{3}';
test2=alPre{4}';
connDat(1).synList=allah2v;

allVin=find(ismember(curTis.syn.edges(:,1),allVG3cids));
allVout=find(ismember(curTis.syn.edges(:,2),allVG3cids));

connDat(7).name='all input to VG3';
connDat(7).synList=allVin;
connDat(8).name='all VG3 outputs';
connDat(8).synList=allVout;


%% make a test plot
%while 0
TF=figure('Position',[100 100 512 512]); 
TA=axes();
hold on
params=struct();
params.color=[1 0 1];
params.size=10;
[~,paramA]=plotSyns(1,curTis,[]);
synPlts=struct();
cols=[[1 0 1];[0 1 1];[1 1 0];[.25 .25 1];[1 .25 .25];[.25 1 .25];[1 0 1];[0 1 1]];
for i=[7:8]%1:length(connNames)
    curp=plotSyns(connDat(i).synList,curTis,params);
    curp.CData=cols(i,:);
    curp.MarkerFaceColor=cols(i,:);
    synPlts(i).sct=curp;
end




morphs=compareMorph([],allVG3cids,fvDir,'NON');

for j=1:length(morphs.patches)
morphs.patches(j).EdgeColor='none';
morphs.patches(j).FaceAlpha=0.4;
end
legNames={'amacrine to VG3','PACv to VG3','BPC to VG3','VG3 to PACv','VG3 to amacrine','VG3 to RGC'};
%legend(legNames);
%end
set(gcf,'color',[1 1 1]);
%curF=gcf;
%curA=gca;
TA.XLim=lims(1,:);
TA.YLim=lims(2,:);
TA.ZLim=lims(3,:);

axs=gca;
[cursX,cursY,cursZ]=deal(20,144,127);
[curLimX,curLimY,curLimZ]=deal(axs.XLim,axs.YLim,axs.ZLim);
Xspan=max(abs(curLimX-cursX));
Yspan=max(abs(curLimY-cursY));
Zspan=max(abs(curLimZ-cursZ));
axs.XLim=[cursX-(Xspan/2) cursX+(Xspan/2)];
axs.YLim=[cursY-(Yspan/2) cursY+(Yspan/2)];
axs.ZLim=[cursZ-(Zspan/2) cursZ+(Zspan/2)];

for k=1:length(synPlts)
    synPlts(k).sct.MarkerFaceAlpha=0;
    synPlts(k).sct.MarkerEdgeAlpha=0;
end

%lims=[curA.XLim;curA.YLim;curA.ZLim];
framCount=360;
%mov(framCount)= struct('cdata',[],'colormap',[]);
u=1;
while u<361
    view(u-270,0);
    %camroll(-90);
    drawnow
    %filname=sprintf('%smov\\allVG3_rotate_%03d.png',localDir,u);
    %exportgraphics(gcf,filname,'Resolution',150);
    frameIt=uint16(u);
    mov(frameIt)=getframe(TF);
    u=u+1
end
movNam=sprintf('%smov\\allVG3_rotate_overlay',localDir);
v = VideoWriter(movNam);
open(v);
writeVideo(v,mov)
close(v);


