% Headings for table describing VG3 connectome:
% Total synapse count
% Arbor microns
% Input and output counts (bpc in / amc in / unk in / rgc out / amc out / unk out)
% 0.4 (43)	(syn/um) (total)

%% total synapse counts from BPCs and AMCs onto the 6 important VG3s.
% pulled from the tis.mat
curTis=highTis;
vgcSkelList=[2 3 4 5 13 14];
totalInputs=find(ismember(curTis.syn.edges(:,1),vgcSkelList));
%get the ribbon/conv % for each of the VG3
inputStrs=repmat("_",size(totalInputs));
for t=1:length(inputStrs)
    curEdge=curTis.syn.edges(totalInputs(t),:);
    curStr=strjoin([num2str(curEdge(1)), "_", num2str(curEdge(4))]);
    inputStrs(t)=curStr;
end
inputSynTab=tabulate(inputStrs);


%% get the numbers of ID'd BPCs and the %s of each type onto VG3
allPre=cid2type(curTis.syn.edges(:,2),curTis);
allVG3Cids=type2cid({'amc'},{'vgc'},curTis);
allVG3Cids=allVG3Cids{1};
preV=find(ismember(curTis.syn.edges(:,1),allVG3Cids));
bPre=find(allPre{1}'==7);
b2v=intersect(preV,bPre);
b2vTarg=curTis.syn.edges(b2v,1);
tabStr=repmat("",size(b2v));
for u=1:length(b2v)
    tabStr(u)=strjoin([allPre{4}(b2v(u)),"_",b2vTarg(u)]);
end
bpcTypeTab=tabulate(tabStr);
doubleCheck=allPre{4}(b2v)';
% 
% 
% b2vTargStr=num2str(b2vTarg);
% testCell=horzcat(allPre{4}(b2v)',b2vTargCell);
% testStr=strjoin(testCell);
% bpcTypeTab=tabulate(allPre{4}(b2v));

bpcTypeList={'bc3a','bc3b','bc4','bc5o','bc5i','bc5t','xbc','bc6'};
bpcTypeList2={'bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc'};
bpcIDdCidList=type2cid(bpcTypeList2,bpcTypeList,curTis);
[a,b,c,d,e,f,g,h]=deal(bpcIDdCidList{:});
bpcList=vertcat(a,b,c,d,e,f,g,h);
doubleCheck=cid2type(bpcList,curTis);
test1=doubleCheck{4}'

%% 
synPerBpc={};
pairDat={};
results=cell(8,7);
for i=1:length(bpcIDdCidList)
    curCidList=bpcIDdCidList{i};
    typBpcSyns={};
    pairDat{i}=zeros(length(curCidList),length(vgcSkelList));
    for j=1:length(curCidList)
        
        curCid=curCidList(j);
        typBpcSyns{j}=length(find(ismember(curTis.syn.edges(:,1),vgcSkelList)&curTis.syn.edges(:,2)==curCid));
        BpcPair=zeros(length(vgcSkelList),1);
        for k=1:length(vgcSkelList)
            curTarg=vgcSkelList(k);
            BpcPair(k)=length(find(curTis.syn.edges(:,1)==curTarg&curTis.syn.edges(:,2)==curCid));
            pairDat{i}(j,k)=length(find(curTis.syn.edges(:,1)==curTarg&curTis.syn.edges(:,2)==curCid));
        end
    end
    curDat=pairDat{i};
    pp=sum(curDat,2);
    pp2=sum(curDat>0,2);
    pp3=sum(curDat>0,1);
    pp4=sum(curDat(:)>0);
    perPair=pp(pp2>0)./pp2(pp2>0);
    ppsem=mean(perPair)/sqrt(sum(pp2>0)*sum(pp3>0));
    
    fprintf('\n%s, %.0f, %.0f, %.2f (%.3f), %.2f (%.3f)\n', ...
        bpcTypeList{i},length(curCidList),sum(curDat(:)),mean(sum(curDat,2)), ...
        std(sum(curDat,2))/sqrt(length(curCidList)), ...
        mean(perPair),ppsem);
    results(i,:)={bpcTypeList{i},length(curCidList),sum(curDat(:)),mean(sum(curDat,2)), ...
        std(sum(curDat,2))/sqrt(length(curCidList)), ...
        mean(perPair),ppsem};
end


%% in text RGC numbers
% total output synapses
vPre=find(ismember(curTis.syn.edges(:,2),vgcSkelList));
allvPre=find(ismember(curTis.syn.edges(:,2),allVG3Cids));
allPost=cid2type(curTis.syn.edges(:,1),curTis);
aPost=find((allPost{1}==8&allPost{3}~=1)|allPost{1}==0);
rPost=find(allPost{1}==1);
v2r=intersect(vPre,rPost);
v2a=intersect(vPre,aPost);

683/1279
519/1279

% doing the rgc table
rCids=unique(curTis.syn.edges(v2r,1));
rTypes=cid2type(rCids,curTis);
rTypeTab=tabulate(rTypes{4});
%tabulate(allPost{4}(v2r));
rTypesUnique=unique(rTypes{4});
%rTypesUnique{1}=[];
rTypesUnique(ismember(rTypesUnique,{' ','none','unk','off','on'})) = {''};
rTypesUnique=rTypesUnique(~cellfun('isempty',rTypesUnique));
rtm=repmat("rgc",size(rTypesUnique));

rgcIDdCidList=type2cid(rtm,rTypesUnique,curTis);


synPerBpc={};
pairDat={};
rgcresults=cell(length(rTypesUnique),5);
for i=1:length(rTypesUnique)
    curCidList=rgcIDdCidList{i};
    typRgcSyns={};
    pairDat{i}=zeros(length(curCidList),length(vgcSkelList));
    for j=1:length(curCidList)
        
        curCid=curCidList(j);
        typRgcSyns{j}=length(find(ismember(curTis.syn.edges(:,2),vgcSkelList)&curTis.syn.edges(:,1)==curCid));
        BpcPair=zeros(length(vgcSkelList),1);
        for k=1:length(vgcSkelList)
            curTarg=vgcSkelList(k);
            BpcPair(k)=length(find(curTis.syn.edges(:,2)==curTarg&curTis.syn.edges(:,1)==curCid));
            pairDat{i}(j,k)=length(find(curTis.syn.edges(:,2)==curTarg&curTis.syn.edges(:,1)==curCid));
        end
    end
    curDat=pairDat{i};
    pp=sum(curDat,2);
    pp2=sum(curDat>0,2);
    pp3=sum(curDat>0,1);
    perPair=pp(pp2>0)./pp2(pp2>0);
    ppsem=mean(perPair)/sqrt(sum(pp2>0)*sum(pp3>0));
    newPerPairMean=sum(curDat(:))/sum(curDat(:)>0);
    newPerPairSEM=std(curDat(curDat(:)>0))/sqrt(sum(curDat(:)>0));
    
%     fprintf('\n%s, %.0f, %.0f, %.2f (%.3f), %.2f (%.3f)', ...
%         rTypesUnique{i},length(curCidList),sum(curDat(:)),mean(sum(curDat,2)), ...
%         std(sum(curDat,2))/sqrt(length(curCidList)), ...
%         mean(perPair),ppsem);
     fprintf('\n%s, %.0f, %.0f, %.2f (%.3f)', ...
        rTypesUnique{i},length(curCidList),sum(curDat(:)), ...
        newPerPairMean,newPerPairSEM);
%     rgcresults(i,:)={rTypesUnique{i},length(curCidList),sum(curDat(:)), ...
%         mean(sum(curDat,2)), ...
%         std(sum(curDat,2))/sqrt(length(curCidList)), ...
%         mean(perPair),ppsem};
    rgcresults(i,:)={rTypesUnique{i},length(curCidList),sum(curDat(:)), ...
        newPerPairMean,newPerPairSEM};
end

% we were able to identify XX RGCs
sum(cell2mat(rgcresults(:,2)))

% onto which VG3 formed XX synapses
sum(cell2mat(rgcresults(:,3)))

% comprising XX subtypes
size(rgcresults,1)

% with two-thirds? being represented by one or two cells
sum(cell2mat(rgcresults(:,2))<3)

% targets were 4i/4on/etc with XX% (+-xx%)
VGpairDat=cell(length(vgcSkelList),1);
percs=zeros(size(rgcresults,1)+1,length(vgcSkelList));
for u=1:length(vgcSkelList)
    curVdat={};
for k=1:length(pairDat)
   curDat=pairDat{k};
   curVdat{k}=curDat(:,u);
   percs(k,u)=sum(curDat(:,u));
end
    VGpairDat{u}=curVdat;
end
percs(21,:)=sum(percs,1);
percsDiv=percs./percs(21,:);
means=mean(percsDiv,2);
sems=zeros(size(rgcresults,1),1);
for p=1:size(rgcresults,1)
    sems(p)=std(percsDiv(p,:))/sqrt(6);
end

finalResults=horzcat(percsDiv(1:20,:),means(1:20,:),sems);

%% loading stuff
% load allSkels if it's not already here

%% go through the skeletons and get total number of arbor microns
vgcSkelList=[2 3 4 5 13 14];
arborRadCut=1.5;
resultStruct=struct();
for skelIt=1:length(vgcSkelList)
    curSM=allSkels{skelIt}.sm;
    curCid=curSM.cid;
    %need to find the total length of arbor that is reconstructed
    allEdgeLength=curSM.nep.props.edgeLength;
    %figure(); histogram(allEdgeLength);
    arborTot=sum(curSM.nep.props.edgeLength);
    arborNodePos=curSM.nep.pos;
    %figure(); histogram(arborNodePos(:,3));
    arborEdgeRad=curSM.nep.edgeRad;    
    arborNodeRad=curSM.nep.nodeRad;
    %figure(); histogram(arborNodeRad);
    %
    arborTotNoSoma=sum(curSM.nep.props.edgeLength(arborEdgeRad<arborRadCut));
    arborTotNoSoma=arborTotNoSoma;
    %
    testFig=1
    if testFig
        radMult=50;
        nodeCols=repmat([.5 .5 .5],length(arborNodePos(:,1)),1);
   nodeCols(arborNodeRad>arborRadCut,:)=repmat([1 0 1],sum(arborNodeRad>arborRadCut),1);
    figure(); hold on; scatter3(arborNodePos(:,1),arborNodePos(:,2),arborNodePos(:,3),arborNodeRad.*radMult,nodeCols,'.')
    curFV=getFV(curCid,[]);
   plotFV=curFV{1};
   plotFV.vertices=plotFV.vertices(:,[1 2 3]);
   p1=patch(plotFV);
   p1.FaceAlpha=0.2;
   p1.EdgeColor='none';
    end
    resultStruct(skelIt).arborTotNoSoma=arborTotNoSoma;
end
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);

noTypeSyns=[0 0 0 0];
for synIt=1:length(vgcSkelList)
    curCid=vgcSkelList(synIt);
    curSyns=find(any(curTis.syn.edges(:,[1 2])==curCid,2));
    resultStruct(synIt).totalSyn=length(curSyns);
    inputs=find(curTis.syn.edges(:,1)==curCid);
    outputs=find(curTis.syn.edges(:,2)==curCid);
    amcIn=find(curTis.syn.edges(:,1)==curCid & (curTis.syn.edges(:,2)==0 ... 
        | allPreTypes{1}'==8));
    bpcIn=find(curTis.syn.edges(:,1)==curCid & allPreTypes{1}'==7);
    amcOut=find(curTis.syn.edges(:,2)==curCid & (curTis.syn.edges(:,1)==0 ... 
        | allPostTypes{1}'==8));
    rgcOut=find(curTis.syn.edges(:,2)==curCid & allPostTypes{1}'==1);
    resultStruct(synIt).bpcIn=length(bpcIn);
    resultStruct(synIt).amcIn=length(amcIn);
    resultStruct(synIt).amcOut=length(amcOut);
    resultStruct(synIt).rgcOut=length(rgcOut);
    allClassdSyns=vertcat(bpcIn,amcIn,amcOut,rgcOut);
    test2=curTis.syn.edges(curSyns(~ismember(curSyns,allClassdSyns)),:);
    noTypeSyns=vertcat(noTypeSyns,test2);
end