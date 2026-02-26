%% final scripts for the paper

medTis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\medRes\Analysis\fvLibrary\tis.mat');
highTis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\tis.mat');
medTis=medTis.tis;
highTis=highTis.tis;
medFV='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\medRes\Analysis\fvLibrary\';
highFV='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
mainFV='Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\fvLibrary\';
medObj=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\medRes\Merge\dsObj.mat');
highObj=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\dsObj.mat');
medObj=medObj.dsObj;
highObj=highObj.dsObj;


%% get the bpc and vg3 innervation of the important RGCs

curTis=highTis;
alPre=cid2type(highTis.syn.edges(:,2),highTis);
alPost=cid2type(highTis.syn.edges(:,1),highTis);
allVG3cids=highTis.syn.edges(find(alPre{1}==8&alPre{3}==1),2);
allVG3cids=unique(allVG3cids);
preVG3bpcCids=highTis.syn.edges(find(ismember(highTis.syn.edges(:,1),allVG3cids)&alPre{1}'==7),2);
preVG3bpcCids=unique(preVG3bpcCids);
postVGCrgcCids=highTis.syn.edges(find(ismember(highTis.syn.edges(:,2),allVG3cids)&alPost{1}'==1),1);
postVGCrgcCids=unique(postVGCrgcCids);

exRGCs=[2002 3051 3119];
for i=1:length(exRGCs)
    curRGC=exRGCs(i);
    sBPCind=find(ismember(curTis.syn.edges(:,2),preVG3bpcCids)&curTis.syn.edges(:,1)==curRGC&alPre{1}'==7);
    nsBPCind=find(~ismember(curTis.syn.edges(:,2),preVG3bpcCids)&curTis.syn.edges(:,1)==curRGC&alPre{1}'==7);
    [curRGC length(sBPCind) length(nsBPCind)]
    length(sBPCind)/(length(nsBPCind)+length(sBPCind))
end

%% find the number of ribbons with X postsynaptic cells
k=1;
ribsults={};
ribsultDebug={};
ribPartTypes={};
ribPartSubtypes={};
ribPartTypeMat={};
preBpcSubtypes={};
preBpcCid={};
for j=1:length(curTis.syn.synProp)
    curSynProp=curTis.syn.synProp{j};
    if curSynProp.preID & curSynProp.postID
    preType=cid2type(curSynProp.preID,curTis);
    postTypes=cid2type(curSynProp.postID,curTis);
%     if sum(curSynProp.postID==3048)>0
%         curSynProp.name
%     end
    if sum(ismember(curSynProp.postID,allVG3Cids))>0
    if preType{1}==7
        preBpcCid{k}=curSynProp.preID;
        ribsults{k}=length(curSynProp.postID);
        ribPartCids{k}=curSynProp.postID;
        ribsultDebug{k}=curSynProp.name;
        %fprintf('\n%s',curSynProp.name);
        ribPartTypes{k}=postTypes{1};
        ribPartSubtypes{k}=postTypes{3};
        preBpcSubtypes{k}=preType{3};
        ribPartTypeMat{k}=horzcat(postTypes{1}',postTypes{3}');
        k=k+1;
    end
    end
    end
end

%badSubs=[0 1 2 3 4 57 58 59];
badSubs=[0 1 2 3 4 52 53 54 55 58 59 ];

RGCpart=0;
nonRGCpart=0;
typedRGCpart=0;
typedRVrib=zeros(length(ribPartTypes),1);
postSub=repmat("na",length(ribPartTypes),1);
preSub=repmat("na",length(ribPartTypes),1);
for m=1:length(ribPartTypes)
    curPostIDs=ribPartCids{m};
    curPartTypeDat=ribPartTypes(m);
    curPartTypeDat=curPartTypeDat{1};
    curPartTypeMat=ribPartTypeMat{m};
    curPreSub=preBpcSubtypes{m};
    for n=1:size(curPartTypeMat,1)
        curPartCid=curPostIDs(n);
        if ismember(curPartCid,postVGCrgcCids)
            typedRGCpart=typedRGCpart+1;
        end
        curRow=curPartTypeMat(n,:);
        if sum(curRow==[8 1])<2
            if curRow(1)==1
                RGCpart=RGCpart+1;
                if sum(ismember(curRow(2),badSubs))==0
                    typedRVrib(m)=1;
                    postSub(m)=string(curTis.cells.type.subTypeNames{1}(curRow(2)));
                    if ~curPreSub==0
                        preSub(m)=string(curTis.cells.type.subTypeNames{7}(curPreSub));
                    else
                        preSub(m)="zero";
                    end
                    %fprintf('\n%s',string(curTis.cells.type.subTypeNames{1}(curRow(2))));
                end
            else
                nonRGCpart=nonRGCpart+1;
            end
        %fprintf('\ncid%d pvr%d r%d nr%d',curPartCid,typedRGCpart,RGCpart,nonRGCpart);
        %fprintf('\n %d',RGCpart-typedRGCpart);
        end
    end
end

bpcTabs=tabulate(preSub);
rgcTabs=tabulate(postSub);
BRtabs=[preSub+"_"+postSub];
combTabs=tabulate(BRtabs);

preCids=cell2mat(preBpcCid);
preTypeMat=cid2type(preCids,curTis);
preTypeCellStr=preTypeMat{4};
preTypeStrings=repmat("na",length(preTypeCellStr),1);
for n=1:length(preTypeCellStr)
preTypeStrings(n)=preTypeCellStr{n};
end

debug=horzcat(preTypeCellStr',ribsultDebug');

%% checking the state of the 4ow results
v2offa=getClassConn({'amc'},{'vgc'},{'rgc'},{'4ow'},curTis);

v2offaEdges=curTis.syn.edges(v2offa,:);


%% filling in XXs

%We identified 1208 output synapses in the reconstructed VG3 plexus
vpre=find(ismember(curTis.syn.edges(:,2),allVG3cids));
vpreSum=length(vpre);
fprintf('We identified %d output synapses in the reconstructed VG3 plexus',vpreSum);

%The majority of VG3 output synapses targeted other amacrine cells
vpre=find(ismember(curTis.syn.edges(:,2),allVG3cids));
apost1=find(allPostTypes{1}'==8);
apost2=find(allPostTypes{1}'==0|curTis.syn.edges(:,1)==0);
apost=[apost1;apost2];

v2a=intersect(vpre,apost);

vpreSum=length(vpre);
v2aSum=length(v2a);
v2aPerc=v2aSum*100/vpreSum;
fprintf('\nThe majority of VG3 output synapses targeted other amacrine cells (n=%d (%3.1f%%))',v2aSum,v2aPerc);

%We identified 547 (43.6%) probable VG3 to RGC synapses
vpre=find(ismember(curTis.syn.edges(:,2),allVG3cids));
rpost=find(allPostTypes{1}'==1);

v2r=intersect(vpre,rpost);

rpostSum=length(rpost);
v2rSum=length(v2r);
v2rPerc=v2rSum*100/vpreSum;
fprintf('\nWe identified %d (%3.1f%%) probable VG3 to RGC synapses',v2rSum,v2rPerc);

%which synapses are to neither RGC nor AMC?
sanity1=intersect(v2a,v2r);
mysterySyns=setdiff(vpre,[v2a;v2r]);
mystSynTargs=curTis.syn.edges(mysterySyns,1);
unclassCids=tabulate(mystSynTargs);
[g UCord]=sort(unclassCids(:,2),'descend');
unclassCids=unclassCids(UCord,:);
%looks like they are mostly bpc mistakes.
if 0
for u=1:length(mysterySyns)
    curSynID=mysterySyns(u);
    curSynEdges=curTis.syn.edges(curSynID,:);
    curSynPos=curTis.syn.pos(curSynID,:);
    curSynPosVast=curSynPos([2 1 3]).*[250 250 25];
    clipboard('copy',curSynPosVast);
    fprintf('%d to %d\n',curSynEdges(2),curSynEdges(1));
    pause();
end
end
% We were able to reconstruct and identify the subtypes of 43 retinal ganglion
% cells innervated by 289 synapses from the VG3 plexus (Fig XX).
rgcCids=type2cid({'rgc'},{'all'},curTis);
rgcCids=rgcCids{1};
rgcType=cid2type(rgcCids,curTis);
typRind=find(rgcType{1}==1&~ismember(rgcType{3},badSubs));
rgcTypTab=tabulate(rgcType{4}(typRind));

typedRGCsum=sum(cell2mat(rgcTypTab(:,2)));

typedRGCcids=rgcCids(typRind);
v2rtypd=find(ismember(curTis.syn.edges(:,1),typedRGCcids) ...
    & ismember(curTis.syn.edges(:,2),allVG3cids));
v2rtypdSum=length(v2rtypd);
fprintf('We were able to reconstruct and identify the subtypes of %d retinal ganglion cells innervated by %d synapses from the VG3 plexus\n', ...
    typedRGCsum,v2rtypdSum);

%These retinal ganglion cells included 12 subtypes with half of the subtypes
[g, rgcOrd]=sort(cell2mat(rgcTypTab(:,2)),'descend');
rgcTypTabSrtd=rgcTypTab(rgcOrd,:);
numOfRGCtypes=size(rgcTypTabSrtd,1);
fprintf('These retinal ganglion cells included %d subtypes with \n',numOfRGCtypes);

%The type 4 RGCs can be divided into 4i (17.2%,+-1.3% of VG3 to RGC synapses,
%mean and standard error for 6 VG3s), 4on (no synapses detected),
%and 4ow (transient OFF Alpha RGC, 18.1%, +-3%) RGCs. 
%The second most common RGC subtype synapse was to type 5ti RGCs 
%(W3, UHD, small-field transient ON/OFF, 24.7%, +-2.7%) RGCs.

targCids=curTis.syn.edges(v2rtypd,1);
preV=curTis.syn.edges(v2rtypd,2);
targTyp=cid2type(targCids,curTis);
targSubTypeNum=targTyp{3};
targSubType=targTyp{4};
targTypTab=tabulate(targSubType);
rgc4iSum=sum(contains(targSubType,{'4i'}));
rgc4iPerc=rgc4iSum/v2rtypdSum;
mainVGcids=[2 3 4 5 13 14];

rtList={'4i','4ow','5ti','63','37','6sw','8w'};
for t=1:length(rtList)
    curRT=rtList(t);
percs=zeros(length(mainVGcids),1);
for r=1:length(mainVGcids)
    curV=mainVGcids(r);
    curTab=tabulate(targSubType(preV==curV));
    curRgc4iSum=sum(contains(targSubType(preV==curV),curRT));
    curV4iPerc=curRgc4iSum/sum(preV==curV);
    percs(r)=curV4iPerc;
end
fprintf('%s (%3.1f%%,+-%3.1f%% of VG3 to RGC synapses\n',curRT{:},mean(percs)*100,std(percs)*100/sqrt(length(mainVGcids)));
end

%%
testCids=[2002,3051,3119];

for testCidIt=1:length(testCids)
    testCid=testCids(testCidIt);
vPost=find(ismember(highTis.syn.edges(:,1),allVG));
vPre=find(ismember(highTis.syn.edges(:,2),allVG));
toaPost=find(highTis.syn.edges(:,1)==testCid);
bPre=find(allPreInfo{1}==7);
bInput=intersect(bPre,toaPost);
bInCids=highTis.syn.edges(bInput,2);
bInCids=unique(bInCids);
sharedBins=intersect(bInCids,preVG3bpcCids);
length(sharedBins)
length(bInCids)
length(sharedBins)/length(bInCids)
vInput=intersect(vPre,toaPost);
b2vSegIds=highTis.syn.edges(intersect(bPre,vPost),3);
b2rSegIds=highTis.syn.edges(bInput,3);
shared=intersect(b2vSegIds,b2rSegIds);
outp=[ length(shared) ...
    length(vInput) length(bInput) ...
    length(vInput)/(length(bInput)+length(vInput))];

m=0;
for j=1:length(ribPartCids)
    curParts=ribPartCids{j};
    if ismember(testCid,curParts)
        m=m+1;
    end
end
n=0;
for k=1:length(preBpcCid)
    curCid=preBpcCid{k};
    if ismember(curCid,preVG3bpcCids)
        n=n+1;
    end
    
end
%m
%n
end
%%
pairs=repmat("none",length(curTis.syn.edges(:,1)),1);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
postSubs=allPostTypes{4};
for k=1:length(curTis.syn.edges(:,1))
    curPostType=string(postSubs(k));
    curPreCid=string(curTis.syn.edges(k,2));
    curString=curPostType+"_"+curPreCid;
    pairs(k)=curString;
end
pairTabCln=tabulate(pairs(vpre));
pairTab=tabulate(pairs);
[x,pairTabClnSrtd]=sort(pairTabCln(:,1));
pairTabFinal=pairTabCln(pairTabClnSrtd,:);

%%
%%
vPost=find(ismember(curTis.syn.edges(:,1),allVG3cids));
b2vInds=intersect(bPre,vPost);
b2vPreSubs=alPre{4}(b2vInds)';
b2vPreSubTab=tabulate(b2vPreSubs);

pairs=repmat("none",length(curTis.syn.edges(:,1)),1);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
postSubs=allPostTypes{4};
for k=1:length(curTis.syn.edges(:,1))
    curPostType=string(postSubs(k));
    curPreCid=string(curTis.syn.edges(k,2));
    curString=curPostType+"_"+curPreCid;
    pairs(k)=curString;
end
pairTabCln=tabulate(pairs(vpre));
pairTab=tabulate(pairs);
[x,pairTabClnSrtd]=sort(pairTabCln(:,1));
pairTabFinal=pairTabCln(pairTabClnSrtd,:);

%%
%bcstringList=["bc3a","bc3b","bc4","bc5o","bc5i","bc5t","xbc","bc6"];
bpcCidLists=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc',}, ...
        {'bc3a','bc3b','bc4','bc5o','bc5i','bc5t','xbc','bc6'},curTis);

results={};
for k=1:length(bpcCidLists)
    curCidList=bpcCidLists{k};
    synNum=zeros(length(curCidList),1);
    for cidIt=1:length(curCidList)
        curCid=curCidList(cidIt);
        curPre=find(curTis.syn.edges(:,2)==curCid);
        curSyns=intersect(vPost,curPre);
        synNum(cidIt)=length(curSyns);
    end
    results{k,1}=synNum;
    results{k,2}=mean(synNum);
    results{k,3}=std(synNum)/sqrt(length(synNum));
end






