
global glob
load([glob.useFvDir 'tis.mat'])

postT=8;
postS={1};
preT=7;
preS={0};

bpcSyns={};
for i=1:14
    bpcSyns{i}=getClassConn(postT,postS,preT,{i});
    
end


%% clean the syn list for duplicates
for i=1:14
    curSyns=bpcSyns{i};
    for synIt=1:length(curSyns.pos(:,1))
        curPos=curSyns.pos(synIt,:);
        %sameSynIDs=find(curSyns.edges(:,1)==curEdges(1)&curSyns.edges(:,2)==curEdges(2));
        nearbySyns=find(sum(abs(curSyns.pos(:,:)-curPos),2)<.25);
        if ~isempty(nearbySyns)
            for j=1:length(nearbySyns)
                if nearbySyns(j)~=synIt
                    if curSyns.edges(synIt,1:2)==curSyns.edges(nearbySyns(j),1:2)
                        nearbySyns(j)
                        bpcSyns{i}.edges(nearbySyns(j),:)
                        bpcSyns{i}.pos(nearbySyns(j),:)
                        bpcSyns{i}.edges(nearbySyns(j),:)=zeros(size(bpcSyns{i}.edges(nearbySyns(j),:)));
                        bpcSyns{i}.pos(nearbySyns(j),:)=zeros(size(bpcSyns{i}.pos(nearbySyns(j),:)));
                    end
                end
            end
        end
    end
end

%% fig for josh

makefig=1;
if makefig==1
    bpcfig=figure();
    hold on
    for i=1:14
        bpcSubList=unique(bpcSyns{i}.edges(:,2));
        if ~isempty(bpcSubList)
            for j=1:length(bpcSubList)
                curXjitter=j/length(bpcSubList)*0.4-0.2;
                curBpcCid=bpcSubList(j);
                if curBpcCid~=0
                    curBpcTargs=bpcSyns{i}.edges(bpcSyns{i}.edges(:,2)==curBpcCid,1);
                    curBpcTargList=unique(curBpcTargs);
                    for k=1:length(curBpcTargList)
                        curXmjitt=k/length(curBpcTargList)*.08-0.05;
                        plot(i+curXjitter+curXmjitt,length(curBpcTargs(curBpcTargs==curBpcTargList(k))),'o');
                    end
                end
            end
        end
    end
    yticks(1:16);
    ylabel('bpc->vgc syn#');
    xticks(1:14);
    xlabel('bpc subtype');
    xticklabels(tis.cells.type.subTypeNames{7});
end

%% We got the Func
function outDat=getSynHistDat(cidList,typeStr)
%typeStr (on, off; )
outDat=struct;
for i=1:length(cidList)
    curCid=cidList(i);
    curSynDat=getSynDat(curCid);
    inSynIDs=find(curSynDat.edges(:,1)==curCid);
    outSynIDs=find(curSynDat.edges(:,2)==curCid);
    histIn=histogram(curSynDat.zdepth(inSynIDs));
    histOut=histogram(curSynDat.zdepth(outSynIDs));
    outDat.histIn(i)=histIn;
    outDat.histOut(i)=histOut;
end
end

function outputDat=getTypeDat(cidList)
for i=1:length(cidList)
    cellID=find(tis.cells.cids==cidList(i));
end
end

function synDat=getSynDat(cid)
%gets the synapses corresponding to a cid. Returns the edges,
%positions, and z depth data in a structure.
global tis;
global Gparams;
global Iparams;
synDat=struct;
inputSynIDs=find(tis.syn.edges(:,1)==cid);
outputSynIDs=find(tis.syn.edges(:,2)==cid);
synDat.edges=tis.syn.edges([inputSynIDs;outputSynIDs],:);
synDat.pos=tis.syn.pos([inputSynIDs;outputSynIDs],:);
%cidList=[tis.syn.edges(inputSynIDs,2);tis.syn.edges(inputSynIDs,2)];
%outType,subType, inType,subTyps
%typeMat=zeros(length(cidList),4);
%cellTypeIDs=arrayfun( @(x)( find(tis.cells.cids==x) ), cidList , 'UniformOutput', false);

for i=1:length(synDat.edges)
    curEdge=synDat.edges(i,:);
    curPostCid=curEdge(1);
    curPreCid=curEdge(2);
    curPostType=tis.cells.type.typeID(find(tis.cells.cids==curPostCid));
    curPreType=tis.cells.type.typeID(find(tis.cells.cids==curPreCid));
    curPostSubtype=tis.cells.type.subTypeID(find(tis.cells.cids==curPostCid));
    curPreSubtype=tis.cells.type.subTypeID(find(tis.cells.cids==curPreCid));
    synDat.types(i,1:4)=[0 0 0 0];
    if ~isempty(curPostType)
        synDat.types(i,1)=curPostType;
    end
    if ~isempty(curPreType)
        synDat.types(i,3)=curPreType;
    end
    if ~isempty(curPostSubtype)
        synDat.types(i,2)=curPostSubtype;
    end
    if ~isempty(curPreSubtype)
        synDat.types(i,4)=curPreSubtype;
    end
    curPos=synDat.pos(i,:);
    [zG,zI,zSyn] = getIPLdepth2(curPos(3),curPos(2),curPos(1));
    %sprintf('%d',zSyn)
    synDat.zdepth(i)=zSyn;
end
end

%gets a distribution with the synapses, their targets, and the z-depths at
%those points.
function synDepthDat=getSynDatDepthDist(synDat,INLplane,GCLplane)
synDepthDat=struct;
global tis
end

%useful later
%     if plotBool==1
%         plotFilename=sprintf('%i.jpg',cid);
%         plotFilePath=[plotLoc plotFilename];
%         figure();
%         plot(plotDatXpercs,plotDatYpercs);
%         xlim([0 1]);
%         camroll(-90);
%         title(plotFilename);
%         saveas(gcf,plotFilePath);
%         close;
%     end


%get the depth of a point in the IPL from the fitted planes
%test with 160,90,30 and 160,190,30
% x     y	z	IPLdepth
% 160	90	20	0.0914 IPL
% 160	190	20	0.2317  INL
% 100	150	10	1.047   GCL
% 150	150	10	0.9071  IPL

function [zGCL,zINL,IPLdepth] = getIPLdepth2(z,x,y)
global B
zGCL=B(1,1)* x + B(2,1) * y + B(3,1);
zINL=B(1,2)* x + B(2,2) * y + B(3,2);
IPLdepth=abs(zINL-z)/abs(zINL-zGCL);
end

function [zGCL,zINL,IPLdepth] = getIPLdepth(z,x,y,GCparams,INparams)

zGCL=(-x*GCparams(2)-y*GCparams(3)-GCparams(4))/GCparams(1);

zINL=(-x*INparams(2)-y*INparams(3)-INparams(4))/INparams(1);

IPLdepth=abs(z-zINL)/abs(zINL-zGCL);

end

%get the downsampled voxel list for a cid in Y X Z
function voxList=getCidVox(cid)
global tis
global dsObj 
allCidList=tis.obI.nameProps.cids;
X = cellfun(@(m)isequal(m,cid),allCidList(1,:));
obIdList=find(X);
for i=1:length(obIdList)
    curObId=obIdList(i);
    if ~exist('curVox')
        curVox=dsObj(curObId).subs;
    else
        curVox=[curVox;dsObj(curObId).subs];
    end
end
voxList=curVox;
end

%get class-wise connectivity between two types/subtypes
function ourSyns=getClassConn(postType,postSub,preType,preSub)
ourSyns=struct;
global tis
postCidList=getTypeCids(postType,postSub);
preCidList=getTypeCids(preType,preSub);
postAll=postCidList{:};
preAll=preCidList{:};
allEdges=tis.syn.edges;
allPos=tis.syn.pos;
synList=find(ismember(allEdges(:,1),postAll) & ismember(allEdges(:,2), preAll));
ourSyns.edges=tis.syn.edges(synList,:);
ourSyns.pos=tis.syn.pos(synList,:);
end

%get all cids for a certain cell type
function cidLists=getTypeCids(neurType,neurSubtype)
global tis
%this needs at least one neurType, and then a 2D subType cell array with
%the subtypes wanted for each type
cidLists={1:length(neurType)};
for i=1:length(neurType)
    curType=neurType(i);
    Subtypes=neurSubtype{i};
    if Subtypes==0
        curCidList=tis.cids(tis.cells.type.typeID==curType);
    else
        curCidList=tis.cids(tis.cells.type.typeID==curType & ismember(tis.cells.type.subTypeID,Subtypes));
    end
    cidLists{i}=curCidList;
end
end