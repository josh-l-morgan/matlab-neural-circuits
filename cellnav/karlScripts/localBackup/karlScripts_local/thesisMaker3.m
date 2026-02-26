%% updaste the database if needed
updateThesisData(['fvLibrary','dsObj','scripts']);

localDir='C:\work\localDat\';

%% load the standard vars
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

%% Basic Preprocs

HCcolMap=load('HCcolmap.mat');
HCcolMap=HCcolMap.HCcolmap;
pCols=round(HCcolMap./255,3);
pCols=pCols([2:2:end 1:2:end],:);
%figure() ;scatter(1:15,repmat(1,[15 1]),repmat(100,[15 1]),pCols,'filled');%colormap(pCols);
alPre=cid2type(curTis.syn.edges(:,2),curTis);
alPost=cid2type(curTis.syn.edges(:,1),curTis);

%% find the synapses to survey for the PACv questions
ac2vg=find(alPre{1}==8&alPost{1}==8&alPost{3}==1)';
na2vg=find(alPre{1}==0&alPost{1}==8&alPost{3}==1)';
pv2vg=find(alPre{1}==8&alPost{1}==8&alPost{3}==1&alPre{3}==3)';
ac2vg=vertcat(ac2vg,na2vg);
ac2vg=ac2vg(~ismember(ac2vg,pv2vg));

allVG3cids=type2cid({'amc'},{'vgc'},tis);
allVG3cids=allVG3cids{1};

fA=figure();
axA=compareMorph([],allVG3cids,fvDir,'INL');
hold on;
params=struct();
params.color=[1 0 1];
[synPlotAC2VG,params]=plotSyns(ac2vg,tis,params);
synPlotAC2VG.MarkerFaceColor='flat';
synPlotAC2VG.SizeData=10;
params2=params;
params2.color=[0 1 1];
[synPlotPV2VG,params2]=plotSyns(pv2vg,tis,params2);
synPlotPV2VG.MarkerFaceColor='flat';

%Okay, the 25 micron box is 125,125 to 150 150
boxCutter=[125,125,150,150];
plotBox=plot3([20,20,20,20,20],[125,125,150,150,125],[125,150,150,125,125]);
plotBox.LineWidth=2;
inBox=find(tis.syn.pos(:,1)>boxCutter(1)& ...
    tis.syn.pos(:,2)>boxCutter(2)& ...
    tis.syn.pos(:,1)<boxCutter(3)& ...
    tis.syn.pos(:,2)<boxCutter(4));
acInside=intersect(ac2vg,inBox);
pvInside=intersect(pv2vg,inBox);

checkLocs=tis.syn.pos(acInside,:);
checkLocsVast=uint16(checkLocs(:,[2 1 3]).*[250 250 25]);
checkSegids=tis.syn.edges(acInside,3);
checkVtargs=tis.syn.edges(acInside,1);
checkACs=tis.syn.edges(acInside,2);

%% Figure XX - XAC cells
XACproc=1;
if XACproc
    xacCidList=type2cid({'amc'},{'xac'},curTis);
    xacCidList=xacCidList{1};
    xacCidList=xacCidList(~ismember(xacCidList,[6866]));%,6859,6860,6890,6891,6892,6897,6904,6925,6933,7012,7013,7017,7018]));
    xacCidList=[1026,3008,1110,2002];
    xacCidListMed=[1026,3008,1110,2002];
    refFig=figure();
    refAx=gca;
    xacMrph=compareMorph(refAx,xacCidList,fvDir);
    
    xacFig=figure();
    %xacAx=gca();
    %xacMrph=compareMorph(xacAx,xacCidList,fvDir);
    hold on;
    highFVs=getFV(xacCidList,fvDir);
    medFVs=getFV(xacCidListMed,medFVdir);
    adjustMedFV=1;
    if adjustMedFV
        for cFV=1:length(medFVs)
            verts=medFVs{cFV}.vertices;
            newVerts=verts-[230 295 0];
            medFVs{cFV}.vertices=newVerts;
        end
    end
    mFVs={highFVs,medFVs};
    patz={};
    for i=1:length(highFVs)
        for j=1:2
            hFV=mFVs{j}{i};
            hPat=patch(hFV);
            hPat.EdgeColor='none';
            hPat.FaceColor=pCols(i,:);
            patz{i,j}=hPat;
        end
    end
    
    v2x=find(alPre{1}==8&alPre{3}==3&alPost{1}==8&alPost{3}==1)';
    x2v=find(alPre{1}==8&alPre{3}==1&alPost{1}==8&alPost{3}==3)';
    synLists={v2x,x2v};
    synCid=1026;
    synLists{1}=find(curTis.syn.edges(:,2)==synCid);
    synLists{2}=find(curTis.syn.edges(:,1)==synCid);
    cols=[[1 0 1];[0 1 1]];
    scatSize=50;
    scatz={};
    for i=1:length(synLists)
        curList=synLists{i};
        scatz{i}=scatter3(curTis.syn.pos(curList,1),curTis.syn.pos(curList,2),curTis.syn.pos(curList,3), ...
            scatSize,cols(i,:),'o','filled');
    end
    
    
 
    
end




%%
allXacCids=type2cid({'amc'},{'xac'},curTis);
allXacCids=allXacCids{1};
allXacVox=getCidVox(allXacCids,1,dsObj,curTis);
allXacSize=cellfun(@size,allXacVox,'UniformOutput',false);
allXacSize=cell2mat(allXacSize);
allXacSize=allXacSize(1:2:end);
allXacSynNums=zeros(length(allXacCids),2);
for i=1:length(allXacCids)
    ins=find(curTis.syn.edges(:,1)==allXacCids(i));
    outs=find(curTis.syn.edges(:,2)==allXacCids(i));
    allXacSynNums(i,:)=[length(ins) length(outs)];
end
xacDat=horzcat(allXacCids,allXacSize',allXacSynNums);
%%
%% 
fA=figure();
hold on;
%%% trimtest
clf
%cents=[[40 115 150];[19 125 149];[19 148 125];[20 157 85];[18 154 96]; ...
%    [22 178 104];[18 154 96];[22 132 96]];
%rads=[[6 12 6];[5 25 25];[15 25 25];[1 10 10];[4 20 20]; ...
%    [10 30 30];[10 15 5];[10 10 10]];
%cents=loxTest;
%rads=repmat([3 3 3],size(cents,1),1);

cents=[[40 115 150]];
rads=[[6 12 6]];
tCid=1026;
curFV=getFV(tCid,fvDir);
curFV=curFV{1};
curFV=flipFV(curFV);
curFVtrm=trimFV(curFV,cents,rads);

pB=patch(curFV);
pA=patch(curFVtrm);
pA.LineStyle='none';
pB.LineStyle='none';
pA.FaceColor=[0 1 1];
pB.FaceColor=[1 0 1];
view(-90,0)
allLox=[0 0 0];
loxNlox={};
roundNum=1;
%% get clicks
while 0
rotate3d off
datacursormode off

set(gcf, 'KeyPressFcn', @keyPressFcn);
while true
    oldTestVar=testVar;
    testVar = getappdata(0, 'Location');
    exitBool = getappdata(0, 'Exit');
    if testVar~=oldTestVar
        allLox=[allLox;testVar];
    end
    if ~exitBool
        break
    end
    pause(0.02);
end
loxNlox{roundNum}=allLox;
roundNum=roundNum+1;
end
%% 

%rotate3d off
%datacursormode off
%cursorInfo = getCursorInfo(dcm);
%set(gcf,'ButtonDownFcn',@ButtFn)
dcm=datacursormode(gcf);
dcm.DisplayStyle='window';
totalSum=0;
testVar=[0 0 0];
while totalSum<500
    oldTestVar=testVar;
    %testVar = getappdata(a, 'Location');
    curse = getCursorInfo(dcm);
    if ~isempty(curse)
    testVar = curse(1).Position;
    if testVar~=oldTestVar
        allLox=[allLox;testVar];
        testVar;
    end
    end
    pause(0.1);
    totalSum=totalSum+1
end
%% 
hold on
scatter3(allLox(2:end,1),allLox(2:end,2),allLox(2:end,3));

%pts=drawpolyline
%curFVtrm=curFV;
%%
bpcCidsOFF=type2cid({'bpc','bpc','bpc'},{'bc5o','bc5i','bc5t'},curTis);
bpcCidsON=type2cid({'bpc','bpc','bpc'},{'bc3a','bc3b','bc4'},curTis);
bpcCidsOFF=bpcCidsOFF{:};
bpcCidsON=bpcCidsON{:};
bpcPatchesOFF=getFV(bpcCidsOFF,fvDir);
bpcPatchesON=getFV(bpcCidsON,fvDir);

medFVs=getFV(tCid,medFVdir);

adjustMedFV=1;
if adjustMedFV
    for cFV=1:length(medFVs)
        verts=medFVs{cFV}.vertices;
        newVerts=verts-[230 295 0];
        medFVs{cFV}.vertices=newVerts;
    end
end

mFVs={highFVs,medFVs};
medFV=medFVs{1};

%curFVtrm.vertices=curFVtrm.vertices(:,[2 1 3]);
%curFV.vertices=curFV.vertices(:,[3 1 2]);
medFV.vertices=medFV.vertices(:,[3 1 2]);


f1=figure('Position',[1100 100 800 800]); 
hold on; 
bpcsOFF=compareMorph([],bpcCidsOFF,fvDir,'NON');
bpcsON=compareMorph([],bpcCidsON,fvDir,'NON');
for k=1:length(bpcsOFF.patches)
    bpcsOFF.patches(k).FaceColor=[1 0 0];
    bpcsOFF.patches(k).FaceAlpha=0.06;
end
for k=1:length(bpcsON.patches)
    bpcsON.patches(k).FaceColor=[0 1 0];
    bpcsON.patches(k).FaceAlpha=0.06;
end

%p1=patch(curFV,'EdgeColor','none','FaceColor',[1 0 1],'FaceAlpha',0.2);
p2=patch(curFVtrm,'EdgeColor','none','FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.2);
p3=patch(medFV,'EdgeColor','none','FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.1);


b2x=getClassConn({'bpc'},{'all'},{'amc'},{'xac'},curTis);
v2x=getClassConn({'amc'},{'vgc'},{'amc'},{'xac'},curTis);
a2x=getClassConn({'amc'},{'all'},{'amc'},{'xac'},curTis);
x2v=getClassConn({'amc'},{'xac'},{'amc'},{'vgc'},curTis);
x2o=getClassConn({'amc'},{'xac'},{'all'},{'all'},curTis);
a2x=a2x(~ismember(a2x,v2x));
x2o=x2o(~ismember(x2o,x2v));

curSynsOut=find(curTis.syn.edges(:,2)==tCid);
curSynsIn=find(curTis.syn.edges(:,1)==tCid);

b2x=intersect(b2x,curSynsIn);
v2x=intersect(v2x,curSynsIn);
a2x=intersect(a2x,curSynsIn);
x2v=intersect(x2v,curSynsOut);
x2o=intersect(x2o,curSynsOut);

synLists={b2x,v2x,a2x,x2v,x2o};
colLists=[[0 0 1];[1 0 1];[0 1 0];[0 1 1];[1 1 0]];
symList={'v','s','d','o','o'};
sizeList=[25,50,25,50,25];
legs={'','','','','','','','','','','','','','','','','', ...
    '','','','','','','','','','','','','','','','','','', ...
    '','','','','','','','','','','','','', ...
    'BPC in','VG3 in','amacrine in','PACv to VG3','PACv to other'};

set(gcf,'color',[1 1 1]);

%scatter3(curTis.syn.pos(curSynsOut,1),curTis.syn.pos(curSynsOut,2), ...
%    curTis.syn.pos(curSynsOut,3),50,'mo','filled');

scatz={};
for j=1:length(synLists)
    curList=synLists{j};
    curCol=colLists(j,:);
    curSym=symList{j};
    scatz{j}=scatter3(curTis.syn.pos(curList,3),curTis.syn.pos(curList,1), ...
        curTis.syn.pos(curList,2),sizeList(j),curCol,curSym,'filled');
end

legend(legs);


%%

if 0
p=0;
clickedLocation=[0 0 0];
%setappdata(0, 'Exit', false);
%clickedLocations=[];
if strcmp(event.Key, 'space')
    % Get the current cursor info
    dcm = datacursormode(gcf);
    cursorInfo = getCursorInfo(dcm);
    if ~isempty(cursorInfo)
        % Extract the clicked location and add it to the array
        clickedLocation = cursorInfo(1).Position;
        setappdata(0, 'Location', clickedLocation);
        %clickedLocations = [clickedLocations; clickedLocation];
        %fprintf('Clicked location added: [X: %.2f, Y: %.2f, Z: %.2f]\n', ...
        %    clickedLocation(1), clickedLocation(2), clickedLocation(3));
    end
elseif strcmp(event.Key, 'escape')
    setappdata(0, 'Exit', true);
end
end


%%

compFig=figure();
hold on
cf1=mightyMorph({[1026,3008,1196,2002],[1026,3008,1196,2002]},{fvDir,medFVdir}, ...
    offsets=[[0 230 295];[0 233 290];[0 230 295];[0 233 290]]);


%% functions
function ButtFn(~, event)
    clickedPt = get(gca,'CurrentPoint');
    setappdata(a, 'Location', clickedPt);
end



function keyPressFcn(~, event)
setappdata(0, 'Exit', true);
key = event.Key;
ptCt=0;
switch key
    case 'space'
        'pressed space'
        dcm = datacursormode(gcf);
        %dcm.Enable='off';
        cursorInfo = getCursorInfo(dcm);
        if ~isempty(cursorInfo)
            clickedLocation = cursorInfo(1).Position;
            setappdata(0, 'Location', clickedLocation);
        end
        ptCt=ptCt+1
    case 'q'
        'pressed q'
        setappdata(0, 'Exit', false);
end
end



function successList=updateThesisData(updateList,updateVolumes,rootDir,localDir)
if ~exist('updateList','var')
    updateList=["tis"];
end

if ~exist('updateVolumes','var')
    updateVolumes=["Final","medRes"];
end


if ~exist('libDir','var')
    rootDir='Y:\karlsRetina\CellNavLibrary_IxQ\';
end

if ~exist('localDir','var')
    localDir='C:\work\localDat\';
end

if ~isfolder(localDir)
    mkdir('localDir');
end



if ~isfolder([localDir '\karlScripts\'])
    mkdir([localDir '\karlScripts\']);
end

for j=1:length(updateVolumes)
    curVol=updateVolumes(j);
    if ~isfolder([localDir 'Volumes\' char(curVol) '\Analysis\fvLibrary\'])
        mkdir([localDir '\Volumes\' char(curVol) '\Analysis\fvLibrary\']);
    end
    
    if ~isfolder([localDir '\Volumes\' char(curVol) '\Merge\'])
        mkdir([localDir '\Volumes\' char(curVol) '\Merge\']);
    end
    
    if ismember('dsObj',updateList)
        copyfile([rootDir 'Volumes\' char(curVol) '\Merge\dsObj.mat'],[localDir 'Volumes\' char(curVol) '\Merge\dsObj.mat'],'f');
    end
end
if ismember('scripts',updateList)
    remScripts=dir('Y:\MATLAB\cellNav\karlScripts\localBackup\karlScripts_local');
    locScripts=dir([localDir '\karlScripts\']);
    for k=1:length(remScripts)
        if remScripts(k).isdir==0
            if sum(ismember({locScripts(:).name},remScripts(k).name))
                if remScripts(k).datenum>locScripts(ismember({locScripts(:).name},remScripts(k).name)).datenum
                    copyfile([remScripts(k).folder '\' remScripts(k).name],[locScripts(1).folder '\' remScripts(k).name],'f');
                end
            else
                copyfile([remScripts(k).folder '\' remScripts(k).name],[locScripts(1).folder '\' remScripts(k).name],'f');
            end
        end
    end
    for k=1:length(locScripts)
        if locScripts(k).isdir==0
            if sum(ismember({remScripts(:).name},locScripts(k).name))
                if locScripts(k).datenum>remScripts(ismember({remScripts(:).name},locScripts(k).name)).datenum
                    copyfile([locScripts(k).folder '\' locScripts(k).name],[remScripts(1).folder '\' locScripts(k).name],'f');
                end
            else
                copyfile([locScripts(k).folder '\' locScripts(k).name],[remScripts(1).folder '\' locScripts(k).name],'f');
            end
        end
    end
end

if ismember('fvLibrary',updateList)
    
    for j=1:length(updateVolumes)
        curVol=updateVolumes(j);
        
        fvDir=[rootDir 'Volumes\' char(curVol) '\Analysis\fvLibrary\'];
        fvDirList=dir(fvDir);
        fvDestList=dir([localDir 'Volumes\' char(curVol) '\Analysis\fvLibrary\']);
        fileList={fvDirList(:).name}';
        destList={fvDestList(:).name}';
        fprintf('\nupdating fvLibrary for volume %s\n',curVol)
        for curFilit=1:numel(fileList)
            curFile=fvDirList(curFilit);
            %curFile.name
            if ~curFile.isdir
                if ~ismember(curFile.name,destList)
                    copyfile([fvDir curFile.name],[localDir 'Volumes\' char(curVol) '\Analysis\fvLibrary\' curFile.name],'f');
                else
                    destFile=fvDestList(find(ismember(destList,curFile.name)));
                    if datetime(destFile.date)<datetime(curFile.date)
                        curFile.name
                        copyfile([fvDir curFile.name],[localDir 'Volumes\' char(curVol) '\Analysis\fvLibrary\' curFile.name],'f');
                        %elseif datenum(destFile.date)==datenum(curFile.date)
                        %    same=1;
                        %elseif datenum(destFile.date)>datenum(curFile.date)
                        %    newer=1;
                    end
                end
            end
        end
        fprintf('\nfinished updating fvLibrary for volume %s\n',curVol)
    end
end

successList=1;
end

