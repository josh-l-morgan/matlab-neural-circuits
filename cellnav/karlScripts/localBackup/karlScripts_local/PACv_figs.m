%% PACv analyses
vgCids=type2cid({'amc'},{'vgc'},curTis);
vgCids=vgCids{1};
xaCids=type2cid({'amc'},{'xac'},curTis);
xaCids=xaCids{1};
rgCids=type2cid({'rgc'},{'all'},curTis);
rgCids=rgCids{1};

v2x=getClassConn({'amc'},{'vgc'},{'amc'},{'xac'},curTis)';
b2x=getClassConn({'bpc'},{'all'},{'amc'},{'xac'},curTis)';
xin=find(ismember(curTis.syn.edges(:,1),xaCids));
xacCid=1026;
all2x1026=find(curTis.syn.edges(:,1)==xacCid);
b2x1026=intersect(all2x1026,b2x);
v2x1026=intersect(all2x1026,v2x);
o2x1026=setdiff(all2x1026,[b2x1026;v2x1026]);

allVout=find(ismember(curTis.syn.edges(:,2),vgCids));
v2r=getClassConn({'amc'},{'vgc'},{'rgc'},{'all'},curTis)';
v2x=getClassConn({'amc'},{'vgc'},{'amc'},{'xac'},curTis)';
v2a=getClassConn({'amc'},{'vgc'},{'amc'},{'all'},curTis)';
v2o=find(ismember(curTis.syn.edges(:,2),vgCids)&curTis.syn.edges(:,1)==0);
v2amc=union(setdiff(v2a,v2x),v2o);
%v2r=getClassConn({'amc'},{'vgc'},{'rgc'},{'all'},curTis)';

allVin=find(ismember(curTis.syn.edges(:,1),vgCids));
x2v=getClassConn({'amc'},{'xac'},{'amc'},{'vgc'},curTis)';
a2v=getClassConn({'amc'},{'all'},{'amc'},{'vgc'},curTis)';
o2v=find(ismember(curTis.syn.edges(:,1),vgCids)&curTis.syn.edges(:,2)==0);
amc2v=union(setdiff(a2v,x2v),o2v);
b2v=getClassConn({'bpc'},{'all'},{'amc'},{'vgc'},curTis)';

x2all=getClassConn({'amc'},{'xac'},{'amc'},{'all'},curTis)';

x2o=setdiff(x2all,x2v);
curTis.syn.edges(x2o,1);

%more pies! PIES PIES PIES PIES
%VG3 output pie
voutpie=[length(v2r),length(v2x),length(v2a)];
vinpie=[length(b2v),407,969];
xinpie=[length(b2x1026),length(v2x1026),length(o2x1026)];
xoutpie=[6,length(x2v),6];

all2two=find(curTis.syn.edges(:,1)==2);
a2two=find(alPre{3}(all2two)~=7);
457*.295
%% 
fpi=figure();
tlo=tiledlayout('flow');
nexttile

p=pie(vinpie,[0 2 0],'%.1f%%');
hold on
title(['VGC inputs: BPC=',num2str(length(b2v)),',PACv=',num2str(length(x2v)),',amacrine=',num2str(length(amc2v))]);
nexttile

p=pie(voutpie,[0 2 0],'%.1f%%');%,{'RGC','PACv','amacrine'});
hold on
title(['VGC outputs: RGC=',num2str(length(v2r)),',PACv=',num2str(length(v2x)),',amacrine=',num2str(length(v2a))]);
nexttile

p=pie(xinpie,[0 2 0],'%.1f%%');%,{'BPC','VGC','amacrine'});
hold on
title(['PACv inputs: BPC=',num2str(length(b2x1026)),',VG3=',num2str(length(v2x1026)),',amacrine=',num2str(length(o2x1026))]);
nexttile

p=pie(xoutpie,[2 0 0],'%.1f%%');%,{'amacrine','VGC','???'});
hold on
title(['PACv outputs: amacrine=',num2str(6),',PACv=',num2str(length(x2v)),',???=',num2str(6)]);

pieCol=ones(10,3);
colormap(pieCol);

%% 
xacCid=1026;
xIn=find(curTis.syn.edges(:,1)==xacCid);
b2x=getClassConn({'bpc'},{'all'},{'amc'},{'xac'},curTis)';
b2xCid=intersect(xIn,b2x);
%b2xCid=b2x;
b2xBCid=curTis.syn.edges(b2xCid,2);

b2v=getClassConn({'bpc'},{'all'},{'amc'},{'vgc'},curTis)';
b2vBCid=curTis.syn.edges(b2v,2);

sharedBPC=intersect(b2xBCid,b2vBCid);
b2xShared=b2xCid(find(ismember(b2xBCid,sharedBPC)));
b2vShared=b2v(find(ismember(b2vBCid,sharedBPC)));
b2xNotShared=b2xCid(find(~ismember(b2xBCid,sharedBPC)));
b2vNotShared=b2v(find(~ismember(b2vBCid,sharedBPC)));

XVsharedBfig=figure();
crax=gca;
hold on
xacMorph=compareMorph(1026,fvDir,labLoc='INL');
vgcMorph=compareMorph([2 3 4 5 13 14 10],fvDir,labLoc='INL');
xacMorph.patches(:).FaceColor=[1 .3 1];
for p=1:length(vgcMorph.patches)
    vgcMorph.patches(p).FaceColor=[0,.7,.7];
    vgcMorph.patches(p).FaceAlpha=.15;
end

scatter3(curTis.syn.pos(b2xNotShared,3),curTis.syn.pos(b2xNotShared,1),curTis.syn.pos(b2xNotShared,2),75,[0 1 0],'o','filled')
scatter3(curTis.syn.pos(b2xShared,3),curTis.syn.pos(b2xShared,1),curTis.syn.pos(b2xShared,2),75,[0 0 1],'o','filled')
scatter3(curTis.syn.pos(b2vNotShared,3),curTis.syn.pos(b2vNotShared,1),curTis.syn.pos(b2vNotShared,2),75,[0 1 0],'d','filled')
scatter3(curTis.syn.pos(b2vShared,3),curTis.syn.pos(b2vShared,1),curTis.syn.pos(b2vShared,2),75,[0 0 1],'d','filled')

[~,~,b2xDeps]=getIPLdepth(curTis.syn.pos(b2x,3),curTis.syn.pos(b2x,1),curTis.syn.pos(b2x,2),[],[]);
[~,~,b2vDeps]=getIPLdepth(curTis.syn.pos(b2v,3),curTis.syn.pos(b2v,1),curTis.syn.pos(b2v,2),[],[]);

bnam=curTis.cells.type.subTypeNames{7}'

OFFlist=[1 2 3 4 5 15 18 19];
ONlist=[6 7 8 9 10 11 12 14 17 21 20 24];
b2xOFF=b2x(ismember(alPre{3}(b2x),OFFlist));
b2xON=b2x(ismember(alPre{3}(b2x),ONlist));
b2vOFF=b2v(ismember(alPre{3}(b2v),OFFlist));
b2vON=b2v(ismember(alPre{3}(b2v),ONlist));

synLists={b2x,b2xOFF,b2xON,b2v,b2vOFF,b2vON};
 %%
 
areaCols=[[1 0 1];[.6 .6 .6];[0 1 1]];
areaLinStyl={'--','--','--',':',':',':'};
divs=[99 99 99 499 499 499];
figure(); hold on
selektor=[2 3 5 6];
for t=1:length(selektor)
    nt=selektor(t);
    cSyn=synLists{nt};
    [~,~,cDeps]=getIPLdepth(curTis.syn.pos(cSyn,3),curTis.syn.pos(cSyn,1),curTis.syn.pos(cSyn,2),[],[]);
    cHist=histcounts(cDeps,[0:.025:1]);
    %cHist=cHist./divs(t);
    cArea=area(cHist);
    cArea.FaceAlpha=.1;
    cArea.LineStyle=areaLinStyl{nt};
    cArea.FaceColor=areaCols(mod(nt,3)+1,:);
    areas{nt}=cArea;
end
%legend({'b2x','b2xOFF','b2xON','b2v','b2vOFF','b2vON'})
legend({'b2xOFF','b2xON','b2vOFF','b2vON'})


 %%
synLists={b2x,b2xOFF,b2xON,b2v,b2vOFF,b2vON};
 %% Checking the depths of the amacrine inputs for the 
 % estimating the number of pacv inputs to VG3

synLists={a2v,o2v,x2v}; 
 
areaCols=[[1 0 1];[.6 .6 .6];[0 1 1]];
areaLinStyl={'--','--','--',':',':',':'};
divs=[99 99 99 499 499 499];
figure(); hold on
selektor=[1 2 3];
for t=1:length(selektor)
    nt=selektor(t);
    cSyn=synLists{nt};
    [~,~,cDeps]=getIPLdepth(curTis.syn.pos(cSyn,3),curTis.syn.pos(cSyn,1),curTis.syn.pos(cSyn,2),[],[]);
    cHist=histcounts(cDeps,[0:.025:1]);
    %cHist=cHist./divs(t);
    cArea=area(cHist);
    cArea.FaceAlpha=.1;
    cArea.LineStyle=areaLinStyl{nt};
    cArea.FaceColor=areaCols(mod(nt,3)+1,:);
    areas{nt}=cArea;
end
%legend({'b2x','b2xOFF','b2xON','b2v','b2vOFF','b2vON'})
legend({'named AC in','unk in','PACv in'})



%% figure 0
compFig=figure();
hold on
lit=camlight;


%clf
%cf1=mightyMorph({[1026,1180,3008,1196,2002,1110],[1026,1180,3008,1196,2002,1110]},{fvDir,medFVdir}, ...
%    offsets=[[0 230 295];[0 233 290];[0 233 290];[0 230 295];[0 233 290];[0 233 290]]);


cf2=mightyMorph({[1026;3007;2003;3049;3223;3008],[1026;3007;2003;3049;3223;3008]},{fvDir,medFVdir});
line10=plot3([5 5],[150 160],[-10 -10]);
line50=plot3([0 0],[150 200],[-15 -15]);
citer=1;
cList=[2 7 9 4 5 12];
cList=[2 4 12 7 9 14];
for j=1:size(cf2.patches,2)
    for k=1:size(cf2.patches,1)
        cf2.patches{k,j}.FaceColor=pCols(cList(citer),:);
        cf2.patches{k,j}.FaceAlpha=1;
    end
    citer=citer+1;
end

lit.Position=[250 150 100];
lit.Style='ambient';

cf3=compareMorph([],[3119 1070 6051 3051 3102] ,fvDir,'NON');
view(90,88);

an=85;
view(90,an);
an=an+.1;



%% figure 1
pacCid=1026;
inPAC=find(curTis.syn.edges(:,1)==pacCid);
outPAC=find(curTis.syn.edges(:,2)==pacCid);
prePACtypeDat=alPre{1}(inPAC);
prePACsubtypeDat=alPre{3}(inPAC);
prePACtypeDatStr=alPre{2}(inPAC)';
prePACsubtypeDatStr=alPre{4}(inPAC)';
ppsdCat=horzcat(prePACtypeDatStr,prePACsubtypeDatStr);

inPAClocs=curTis.syn.pos(inPAC,:);
[~,~,inPACdepths]=getIPLdepth(inPAClocs(:,3),inPAClocs(:,1),inPAClocs(:,2),[],[]);
[~,~,inPACdepths2]=getIPLdepth(inPAClocs(:,3),inPAClocs(:,2),inPAClocs(:,1),[],[]);
%[~,inPACdepths3]=getIPLdepth(inPAClocs(:,3),inPAClocs(:,1),inPAClocs(:,2),[],[]);
%figure(); scatter(inPAClocs(:,3),inPACdepths,50,HCcolMap(mod(prePACsubtypeDat,15)+1,:)./255,'filled');

cents=[[40 115 150]];
rads=[[6 12 6]];
pacCid=1026;
curFV=getFV(pacCid,fvDir);
curFV=curFV{1};
curFV=flipFV(curFV);
curFVtrm=trimFV(curFV,cents,rads);

[~,~,vertDep]=getIPLdepth(curFVtrm.vertices(:,1),curFVtrm.vertices(:,2),curFVtrm.vertices(:,3),[],[]);
[~,~,synDep]=getIPLdepth(curTis.syn.pos(:,3),curTis.syn.pos(:,1),curTis.syn.pos(:,2),[],[]);
vertDep=vertDep.*40;
synDep=synDep.*40;
curFVtrm.vertices(:,1)=vertDep;
%[19.1,170.1,57.15] vertex
%[112.784,132.856,24.92] synapse
pf1=figure('Position',[1100 100 800 800]); 
p1=patch(curFVtrm);
view(-90,0);
p1.FaceAlpha=0.2;
p1.FaceColor=[.4 .4 .4];
p1.EdgeColor='none';

b2x=getClassConn({'bpc'},{'all'},{'amc'},{'xac'},curTis);
v2x=getClassConn({'amc'},{'vgc'},{'amc'},{'xac'},curTis);
a2x=getClassConn({'amc'},{'all'},{'amc'},{'xac'},curTis);
x2v=getClassConn({'amc'},{'xac'},{'amc'},{'vgc'},curTis);
x2o=getClassConn({'amc'},{'xac'},{'all'},{'all'},curTis);
a2x=a2x(~ismember(a2x,v2x));
x2o=x2o(~ismember(x2o,x2v));

curSynsOut=find(curTis.syn.edges(:,2)==pacCid);
curSynsIn=find(curTis.syn.edges(:,1)==pacCid);

b2x=intersect(b2x,curSynsIn);
v2x=intersect(v2x,curSynsIn);
a2x=intersect(a2x,curSynsIn);
x2v=intersect(x2v,curSynsOut);
x2o=intersect(x2o,curSynsOut);

vCidIn=curTis.syn.edges(v2x,2);
vCidOut=curTis.syn.edges(x2v,1);
tabIn=tabulate(vCidIn);
tabOut=tabulate(vCidOut);
tabIn=sortrows(tabIn,2,'descend');
tabOut=sortrows(tabOut,2,'descend');

vCidIn2=curTis.syn.edges(v2x,1);
vCidOut2=curTis.syn.edges(x2v,2);
tabIn2=tabulate(vCidIn2);
tabOut2=tabulate(vCidOut2);
tabIn2=sortrows(tabIn2,2,'descend');
tabOut2=sortrows(tabOut2,2,'descend');


synLists={b2x,v2x,a2x,x2v,x2o};
synCounts=cellfun(@length,synLists);
colLists=[[0 1 1];[1 0 1];[0 1 0];[0 1 1];[1 1 0]];
%colLists=pCols([3 4 7 8 5],:);
colLists=pCols([4 3 10 14 2],:);
symList={'v','s','d','o','o'};
sizeList=[100,100,100,75,75];
%legs={'','','','','','','','','','','','','','','','','', ...
%    '','','','','','','','','','','','','','','','','','', ...
%    '','','','','','','','','','','','','', ...
legs={'','BPC in','VG3 in','amacrine in','PACv to VG3','PACv to other'};
set(gcf,'color',[1 1 1]);
hold on;
scatz={};
for j=1:5%length(synLists)
    curList=synLists{j};
    curCol=colLists(j,:);
    curSym=symList{j};
    scatz{j}=scatter3(synDep(curList),curTis.syn.pos(curList,1), ...
        curTis.syn.pos(curList,2),sizeList(j),curCol,curSym,'filled');
end
legend(legs);
set(gca,'visible','off');
ax=gca;
ax.Position=[0 0 1 1];
ax.ZLim=[50 200];
ax.XLim=ax.ZLim-100;
ax.YLim=ax.ZLim;
view(-90,0)
%exportgraphics(ax,'front.eps','ContentType','vector')
view(0,0)
%exportgraphics(ax,'side.eps','ContentType','vector')
%exportgraphics(gcf,'labs.pdf','ContentType','vector')
%% time for pie
pf2=figure();
pi1=pie(synCounts(1:3));
pi1.Visible
ax=gca;
ax.Colormap=colLists([1 2 3],:);
%exportgraphics(gcf,'pie.pdf','ContentType','vector')


%% Do the same thing for a VG
while 0
    
%x2v    
    
pacCid=2;
inPAC=find(curTis.syn.edges(:,1)==pacCid);
outPAC=find(curTis.syn.edges(:,2)==pacCid);
prePACtypeDat=alPre{1}(inPAC);
prePACsubtypeDat=alPre{3}(inPAC);
prePACtypeDatStr=alPre{2}(inPAC)';
prePACsubtypeDatStr=alPre{4}(inPAC)';
ppsdCat=horzcat(prePACtypeDatStr,prePACsubtypeDatStr);

inPAClocs=curTis.syn.pos(inPAC,:);
[~,~,inPACdepths]=getIPLdepth(inPAClocs(:,3),inPAClocs(:,1),inPAClocs(:,2),[],[]);
[~,~,inPACdepths2]=getIPLdepth(inPAClocs(:,3),inPAClocs(:,2),inPAClocs(:,1),[],[]);
%[~,inPACdepths3]=getIPLdepth(inPAClocs(:,3),inPAClocs(:,1),inPAClocs(:,2),[],[]);
%figure(); scatter(inPAClocs(:,3),inPACdepths,50,HCcolMap(mod(prePACsubtypeDat,15)+1,:)./255,'filled');

%cents=[[40 115 150]];
%rads=[[6 12 6]];
%pacCid=1026;
curFV=getFV(pacCid,fvDir);
curFV=curFV{1};
curFV=flipFV(curFV);
curFVtrm=curFV;
%curFVtrm=trimFV(curFV,cents,rads);

[~,~,vertDep]=getIPLdepth(curFVtrm.vertices(:,1),curFVtrm.vertices(:,2),curFVtrm.vertices(:,3),[],[]);
[~,~,synDep]=getIPLdepth(curTis.syn.pos(:,3),curTis.syn.pos(:,1),curTis.syn.pos(:,2),[],[]);
vertDep=vertDep.*40;
synDep=synDep.*40;
curFVtrm.vertices(:,1)=vertDep;
%[19.1,170.1,57.15] vertex
%[112.784,132.856,24.92] synapse
pf1=figure('Position',[1100 100 800 800]); 
p1=patch(curFVtrm);
view(-90,0);
p1.FaceAlpha=0.2;
p1.FaceColor=[.4 .4 .4];
p1.EdgeColor='none';

b2x=getClassConn({'bpc'},{'all'},{'amc'},{'xac'},curTis);
v2x=getClassConn({'amc'},{'vgc'},{'amc'},{'xac'},curTis);
a2x=getClassConn({'amc'},{'all'},{'amc'},{'xac'},curTis);
x2v=getClassConn({'amc'},{'xac'},{'amc'},{'vgc'},curTis);
x2o=getClassConn({'amc'},{'xac'},{'all'},{'all'},curTis);
u2v=getClassConn({'amc'},{'unk'},{'amc'},{'vgc'},curTis);
o2v=getClassConn({'all'},{'all'},{'amc'},{'vgc'},curTis);
a2v=getClassConn({'amc'},{'all'},{'amc'},{'vgc'},curTis);
o2v=o2v(~ismember(o2v,a2v));
a2v=a2v(~ismember(a2v,x2v));
a2x=a2x(~ismember(a2x,v2x));
x2o=x2o(~ismember(x2o,x2v));

curSynsOut=find(curTis.syn.edges(:,2)==pacCid);
curSynsIn=find(curTis.syn.edges(:,1)==pacCid);

b2x=intersect(b2x,curSynsIn);
v2x=intersect(v2x,curSynsIn);
a2x=intersect(a2x,curSynsIn);
x2v=intersect(x2v,curSynsOut);
x2o=intersect(x2o,curSynsOut);

vCidIn=curTis.syn.edges(v2x,2);
vCidOut=curTis.syn.edges(x2v,1);
tabIn=tabulate(vCidIn);
tabOut=tabulate(vCidOut);
tabIn=sortrows(tabIn,2,'descend');
tabOut=sortrows(tabOut,2,'descend');

x2vTable=sortrows(tabulate(curTis.syn.edges(x2v,1)),2,'descend');


synLists={b2x,v2x,a2x,x2v,x2o};
synCounts=cellfun(@length,synLists);
colLists=[[0 1 1];[1 0 1];[0 1 0];[0 1 1];[1 1 0]];
%colLists=pCols([3 4 7 8 5],:);
colLists=pCols([4 3 10 14 2],:);
symList={'v','s','d','o','o'};
sizeList=[100,100,100,75,75];
%legs={'','','','','','','','','','','','','','','','','', ...
%    '','','','','','','','','','','','','','','','','','', ...
%    '','','','','','','','','','','','','', ...
legs={'','BPC in','VG3 in','amacrine in','PACv to VG3','PACv to other'};
set(gcf,'color',[1 1 1]);
hold on;
scatz={};
for j=1:5%length(synLists)
    curList=synLists{j};
    curCol=colLists(j,:);
    curSym=symList{j};
    scatz{j}=scatter3(synDep(curList),curTis.syn.pos(curList,1), ...
        curTis.syn.pos(curList,2),sizeList(j),curCol,curSym,'filled');
end
legend(legs);
set(gca,'visible','off');
ax=gca;
ax.Position=[0 0 1 1];
ax.ZLim=[50 200];
ax.XLim=ax.ZLim-100;
ax.YLim=ax.ZLim;
view(-90,0)
exportgraphics(ax,'front.eps','ContentType','vector')
view(0,0)
exportgraphics(ax,'side.eps','ContentType','vector')
exportgraphics(gcf,'labs.pdf','ContentType','vector')
%% time for pie
pf2=figure();
pi1=pie(synCounts(1:3));
pi1.Visible
ax=gca;
ax.Colormap=colLists([1 2 3],:);
exportgraphics(gcf,'pie.pdf','ContentType','vector')


end




%% Chapter4 figures
fig42=figure('Position',[1000 100 800 800]);
mm1=mightyMorph({xaCidsDend,xaCidsDend},{fvDir,medFVdir});
set(gcf,'color',[1 1 1]);
hold on
%lit=camlight;
%cf2=mightyMorph({[1026;3007;2003;3049;3223;3008],[1026;3007;2003;3049;3223;3008]},{fvDir,medFVdir});
line10=plot3([5 5],[150 160],[-10 -10]);
line50=plot3([0 0],[150 200],[-15 -15]);
line10b=plot3([0 0],[180 190],[75 75]);
citer=1;
cList=[2 7 9 4 5 12];
cList=[2 6 7 12 4 9 14 3 5 8 13 11 10 1];
for j=1:size(mm1.patches,2)
    for k=1:size(mm1.patches,1)
        mm1.patches{k,j}.FaceColor=pCols(cList(citer),:);
        mm1.patches{k,j}.FaceAlpha=0;
    end
    citer=citer+1;
end
tf=gcf;
mm1.axes.Position=[0 0 1 1];
exportgraphics(gca,'pac_wide_top.png','Resolution',2400);
mm1.patches{1,1}.FaceAlpha=0;
exportgraphics(gca,'pac_wide_top_no1206.png','Resolution',2400);
mm1.patches{1,1}.FaceAlpha=0.8;
mm1.axes.XLimMode='manual';
mm1.axes.YLimMode='manual';
mm1.axes.ZLimMode='manual';
mm1.axes.CameraTargetMode='manual';
mm1.axes.CameraUpVectorMode='manual';
mm1.axes.CameraPositionMode='manual';
xlim=mm1.axes.XLim; ylim=mm1.axes.YLim; zlim=mm1.axes.ZLim;
%mm1.axes.DataAspectRatio=[1,1,1];
mm1.axes.Clipping='on';
mm1.axes.XLim=[0 30];
mm1.axes.YLim=ylim;
mm1.axes.ZLim=zlim;
exportgraphics(gca,'pac_wide_top_nosoma.png','Resolution',2400);
vIN=find(ismember(curTis.syn.edges(:,1),xaCidsDend)&alPre{1}'==8&alPre{3}'==1);
ps1=plotSyns(vIN,curTis,depth=0);
ps1.Marker='d';
ps1.MarkerFaceColor=[0 .5 1];
ps1.MarkerEdgeColor=[0 0 1];
ps1.SizeData=30;
exportgraphics(gca,'pac_wide_VGinputs.png','Resolution',2400);





%%

%lit.Position=[250 150 100];
%lit.Style='ambient';


fig43=figure();
cmA=compareMorph(xaCidsDend,fvDir,depth=1);
citer=1;
cList=[2 7 9 4 5 12];
cList=[2 6 7 12 4 9 14 3 5 8 13 11 10 1];
for k=1:size(cmA.patches,2)
    %for k=1:size(cmA.patches,1)
        cmA.patches(k).FaceColor=pCols(cList(citer),:);
        cmA.patches(k).FaceAlpha=.25;
    %end
    citer=citer+1;
end
vIN=find(ismember(curTis.syn.edges(:,1),xaCidsDend)&alPre{1}'==8&alPre{3}'==1);
ps1=plotSyns(vIN,curTis,depth=1);
ps1.Marker='d';
ps1.MarkerFaceColor=[0 .5 1];
ps1.MarkerEdgeColor=[0 0 1];
line10b=plot3([0 0],[180 190],[75 75]);



xaCidsAxon=curTis.syn.edges(find(ismember(curTis.syn.edges(:,2),xaCids)),2);
xaCidsAxon=unique(xaCidsAxon);

cmA=compareMorph(xaCidsAxon,fvDir,depth=0);
set(gcf,'color',[1 1 1]);

citer=1;
cList=[2 7 9 4 5 12];
cList=[2 6 7 12 4 9 14 3 5 8 13 11 1 2 10 6 7 12 4 9 14 3 5 8 13 11 10 1 15 2 6 7 12 4 9];
for k=1:size(cmA.patches,2)
    %for k=1:size(cmA.patches,1)
        cmA.patches(k).FaceColor=pCols(cList(citer),:);
        cmA.patches(k).FaceAlpha=.8;
    %end
    citer=citer+1;
end
vOUT=find(ismember(curTis.syn.edges(:,2),xaCidsAxon)&alPost{1}'==8&alPost{3}'==1);
ps1=plotSyns(vOUT,curTis,depth=1);
ps1.Marker='d';
ps1.MarkerFaceColor=[1 .5 0];
ps1.MarkerEdgeColor=[1 0 0];
ps1.MarkerEdgeAlpha=1;
ps1.MarkerFaceAlpha=1;
line10b=plot3([0 0],[180 190],[75 75]);
exportgraphics(gca,'pac_close_syns.png','Resolution',2400);
