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
vgcSkelList=[2 3 4 5 13 14];

%% get the bpc and vg3 innervation of the important RGCs

allVG=type2cid({'amc'},{'vgc'},curTis);
allVG=allVG{1};
curTis=highTis;
alPre=cid2type(highTis.syn.edges(:,2),highTis);
alPost=cid2type(highTis.syn.edges(:,1),highTis);
allVG3cids=highTis.syn.edges(find(alPre{1}==8&alPre{3}==1),2);
allVG3cids=unique(allVG3cids);
preVG3bpcCids=highTis.syn.edges(find(ismember(highTis.syn.edges(:,1),allVG3cids)&alPre{1}'==7),2);
preVG3bpcCids=unique(preVG3bpcCids);
postVGCrgcCids=highTis.syn.edges(find(ismember(highTis.syn.edges(:,2),allVG3cids)&alPost{1}'==1),1);
postVGCrgcCids=unique(postVGCrgcCids);

% Altogether, we identified 427 synapses between RGCs innervated by the VG3 plexus and bipolar cells that innervate the VG3 plexus
preVB2postVR=find(ismember(curTis.syn.edges(:,2),preVG3bpcCids)&ismember(curTis.syn.edges(:,1),postVGCrgcCids));

% That is way off what it was before
oldObi=load('C:\work\liteLib\obI.mat');
oldObi=oldObi.obI;

highObi=load([highFV 'obI.mat']);
highObi=highObi.obI;

curObi=highObi;
%oldObi;



synCounter=0;
synCounter2=0;
synCounter3=0;
synCounter4=0;
for i=1:length(curObi.nameProps.synProp)
    curStruct=curObi.nameProps.synProp{i};
    curName=curStruct.name;
    curPre=curStruct.preID;
    curPost=curStruct.postID;
    if length(intersect(curPost,allVG))>0 & ismember(curPre,preVG3bpcCids)
        if ismember(curPre,bpcTypdCidList)
            synCounter4=synCounter4+1;
        end
        synCounter3=synCounter3+1;
    end
    if ismember(curPre,preVG3bpcCids) & length(intersect(curPost,postVGCrgcCids))>0
        synCounter=synCounter+1;
        synCounter2=synCounter2+length(intersect(curPost,postVGCrgcCids));
    end
end

bpcTypd=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc'}, ...
    {'bc2','bc3a','bc3b','bc4','bc5o','bc5i','bc5t','xbc','bc6','bc7'},highTis);
bpcTypdCidList=[];
for q=1:length(bpcTypd)
    bpcTypdCidList=[bpcTypdCidList;bpcTypd{q}];
end
testList=intersect(preVG3bpcCids,bpcTypdCidList);
fprintf('We were able to identify the subtype of %.0f bipolar cells\n forming %.0f of the %.0f (%.1f%%) ribbon synapses\n', ...
    length(testList),synCounter4,synCounter3,100*synCounter4/synCounter3);



%%
oldTis=load('C:\work\liteLib\tis.mat');
oldTis=oldTis.tis;
%curTis=oldTis;


% For most RGC subtypes, bipolar-VG3-RGC drive is parallel to direct input from the same bipolar cell types.
exRGCs=[2002 3051 3119];
for i=1:length(exRGCs)
    curRGC=exRGCs(i);
    sBPCind=find(ismember(curTis.syn.edges(:,2),preVG3bpcCids)&curTis.syn.edges(:,1)==curRGC);%&alPre{1}'==7);
    nsBPCind=find(~ismember(curTis.syn.edges(:,2),preVG3bpcCids)&alPre{1}'==7&curTis.syn.edges(:,1)==curRGC);%&alPre{1}'==7);
    vinputInd=find(ismember(curTis.syn.edges(:,2),allVG3cids)&curTis.syn.edges(:,1)==curRGC);
    fprintf('cid %.0f; %.0f from preVG3 BPCs; %.0f from VG3\n',curRGC,length(sBPCind),length(vinputInd))
    fprintf(' %.3f%% of syns, %.0f of %.0f\n',100*length(sBPCind)/(length(nsBPCind)+length(sBPCind)),length(sBPCind),(length(nsBPCind)+length(sBPCind)))
end

length(find(ismember(highTis.syn.edges(:,2),preVG3bpcCids)))
length(find(ismember(oldTis.syn.edges(:,2),preVG3bpcCids)))

length(find(oldTis.syn.edges(:,1)==2002))
length(find(highTis.syn.edges(:,1)==2002))

length(highTis.syn.edges(:,1))
length(oldTis.syn.edges(:,1))

oldFV='C:\work\liteLib\';
% f3=figure();
% compareMorph(f3,2002,highFV)
% newSynPos=highTis.syn.pos(highTis.syn.edges(:,1)==2002,:);
% oldSynPos=oldTis.syn.pos(oldTis.syn.edges(:,1)==2002,:);
% scatter3(newSynPos(:,3),newSynPos(:,1),newSynPos(:,2),50,'mo');
% scatter3(oldSynPos(:,3),oldSynPos(:,1),oldSynPos(:,2),25,'co','filled');

%% find the number of ribbons with X postsynaptic cells
k=1;
ribsults={};
ribsultDebug={};
ribPartTypes={};
ribPartSubtypes={};
ribPartTypeMat={};
preBpcSubtypes={};
preBpcCid={};
ribsWithPostVRgcs=0;
partsOfIDdRGCs=repmat("empty",length(curTis.syn.synProp),1);
partsRGCind=repmat("empty",1,1);
for j=1:length(curTis.syn.synProp)
    curSynProp=curTis.syn.synProp{j};
    if curSynProp.preID & curSynProp.postID
        preType=cid2type(curSynProp.preID,curTis);
        postTypes=cid2type(curSynProp.postID,curTis);
        %     if sum(curSynProp.postID==3048)>0
        %         curSynProp.name
        %     end
        if contains(curSynProp.tag,'rib')
            if sum(ismember(curSynProp.postID,allVG3cids))>0
                %if preType{1}==7
                preBpcCid{k}=curSynProp.preID;
                ribsults{k}=length(curSynProp.postID);
                if length(curSynProp.postID)==1
                    tx1=curSynProp.name;
                    tx2=curTis.obI.colStruc.anchors(curSynProp.segID,:);
                    tx3=curTis.obI.fuse.exportDir{curTis.obI.fuse.obSource(curSynProp.segID)};
                    fprintf('%s ; loc=%.0f %.0f %.0f ; src=%s\n',tx1,tx2(1),tx2(2),tx2(3),tx3)
                    % pause();
                end
                ribPartCids{k}=curSynProp.postID;
                
                ribsultDebug{k}=curSynProp.name;
                %fprintf('\n%s',curSynProp.name);
                ribPartTypes{k}=postTypes{1};
                ribPartSubtypes{k}=postTypes{3};
                preBpcSubtypes{k}=preType{3};
                ribPartTypeMat{k}=horzcat(postTypes{1}',postTypes{3}');
                strs=string([preType{2} preType{4} postTypes{2} postTypes{4}]);
                %strjoin(strs)
                if sum(ismember(curSynProp.postID,postVGCrgcCids))>0
                    ribsWithPostVRgcs=ribsWithPostVRgcs+1;
                    for o=1:length(postTypes{1})
                        if postTypes{1}(o)==1
                            minstrs=string([preType{2} preType{4} postTypes{2}(o) postTypes{4}(o)]);
                            sinstr=strjoin(minstrs);
                            partsRGCind=[partsRGCind(:);sinstr];
                            
                        end
                    end
                    partsOfIDdRGCs(k)=strjoin(strs);
                end
                k=k+1;
            end
            %end
        end
    end
end

identifiableRpart=find(~contains(partsOfIDdRGCs,'empty'));
b3aToVG4i=find(contains(partsOfIDdRGCs,'3a')&contains(partsOfIDdRGCs,'4i'));
b3aToVG4on=find(contains(partsOfIDdRGCs,'3a')&contains(partsOfIDdRGCs,'4on'));

fprintf('\n1.  of %.0f ribbons innervating VG3 and an IDd RGC, \n%.0f were from a 3a to a VG3 and a 4i\nand %.0f were from 3a to VG3 and 4on\n',length(identifiableRpart),length(b3aToVG4i),length(b3aToVG4on));

rgcTab=tabulate(partsRGCind);

ribMat=cell2mat(ribsults);
ribTab=tabulate(ribMat);

fprintf('\n2.  of the ribbons that had a VG3 involved,\n%.0f ribbons had two partners, and %.0f had three\n',ribTab(2,2),ribTab(3,2))




%RGC vs AMC postSyn partners

%number of postSynpartners that are AMC and RGC
test=vertcat(ribPartTypeMat{:});
amPart=sum(test(:,1)==8&test(:,2)~=1);
rgcPart=sum(test(:,1)==1);

%number of postRib RGCs that are postVG3
test2=vertcat(ribPartCids{:});
postRibRGCcids=test2(test(:,1)==1);
triadRGCcids=ismember(postRibRGCcids,postVGCrgcCids);
%mean(triadRGCcids)
%sum(triadRGCcids)
%length(triadRGCcids)
fprintf('\n3.  When an RGC shared a ribbon with a VG3, \n it was also innervated by VG3 %.1f%% (%.0f of %.0f) of the time\n',mean(triadRGCcids)*100,sum(triadRGCcids),length(triadRGCcids))

%% find out where the ribbons went

ribCnt=0;
for i=1:length(curTis.syn.synProp)
    curSynProp=curTis.syn.synProp{i};
    %curSynProp.tag
    if contains(curSynProp.tag,'rib')
        ribCnt=ribCnt+1;
        curSynProp.postID
    end
end

%% debug

badSubs=[0 1 2 3 4 52 53 54 55 58 59 ];
rgcBadSubNams=curTis.cells.type.subTypeNames{1}(badSubs(2:end))

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
apost1=find(alPost{1}'==8);
apost2=find(alPost{1}'==0|curTis.syn.edges(:,1)==0);
apost=[apost1;apost2];

v2a=intersect(vpre,apost);

vpreSum=length(vpre);
v2aSum=length(v2a);
v2aPerc=v2aSum*100/vpreSum;
fprintf('\nThe majority of VG3 output synapses targeted other amacrine cells (n=%d (%3.1f%%))',v2aSum,v2aPerc);

%We identified 547 (43.6%) probable VG3 to RGC synapses
vpre=find(ismember(curTis.syn.edges(:,2),allVG3cids));
rpost=find(alPost{1}'==1);

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

%% We counted bipolar inputs on subregions surrounding t
testCids=[2002,3051,3119];

bboxes{1}=[[0 500];[0 500]];
bboxes{2}=[[0 500];[0 500]];%[[130 172];[110 165]];
bboxes{3}=[[0 500];[0 500]];

if 0
    testFig1=figure();
    compareMorph(testFig1,3119,highFV);
    hold on
    bIn=find(curTis.syn.edges(:,1)==3119&alPre{1}'==7);
    vIn=find(curTis.syn.edges(:,1)==3119&alPre{1}'==8&alPre{3}'==1);
    vIn=intersect(vIn,vPre);
    vInPos=curTis.syn.pos(vIn,:);
    bInPos=curTis.syn.pos(bIn,:);
    %goodVin=vIn(vInPos(:,1)
    scatter3(vInPos(:,3),vInPos(:,1),vInPos(:,2),40,'mo','filled');
    scatter3(bInPos(:,3),bInPos(:,1),bInPos(:,2),40,'co','filled');
    
    
end
testFigs={};
for testCidIt=1:length(testCids)
    testCid=testCids(testCidIt);
    %vPost=find(ismember(highTis.syn.edges(:,1),allVG));
    bPre=find(alPre{1}==7);
    vPre=find(ismember(highTis.syn.edges(:,2),allVG));
    toaPost=find(highTis.syn.edges(:,1)==testCid);
    v2curr=intersect(vPre,toaPost);
    b2curr=intersect(bPre,toaPost);
    v2currPos=curTis.syn.pos(v2curr,:);
    b2currPos=curTis.syn.pos(b2curr,:);
    curBox=bboxes{testCidIt};
    goodV2curr=v2curr(v2currPos(:,1)>curBox(1,1)&v2currPos(:,1)<curBox(1,2)&v2currPos(:,2)>curBox(2,1)&v2currPos(:,2)<curBox(2,2));
    goodB2curr=b2curr(b2currPos(:,1)>curBox(1,1)&b2currPos(:,1)<curBox(1,2)&b2currPos(:,2)>curBox(2,1)&b2currPos(:,2)<curBox(2,2));
    if 1
        testFigs{testCidIt}=figure();
        compareMorph(testFigs{testCidIt},testCid,highFV);
        hold on;
        scatter3(curTis.syn.pos(goodV2curr,3),curTis.syn.pos(goodV2curr,1),curTis.syn.pos(goodV2curr,2),40,'mo','filled');
        scatter3(curTis.syn.pos(goodB2curr,3),curTis.syn.pos(goodB2curr,1),curTis.syn.pos(goodB2curr,2),40,'co','filled');
    end
    fprintf('%.0f, vg3 inputs=%.0f, bpc inputs=%.0f, %.1f%%\n',testCid,length(goodV2curr),length(goodB2curr),length(goodV2curr)*100/(length(goodV2curr)+length(goodB2curr)));
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
alPost=cid2type(curTis.syn.edges(:,1),curTis);
postSubs=alPost{4};
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
alPost=cid2type(curTis.syn.edges(:,1),curTis);
postSubs=alPost{4};
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
bpcCidLists=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc'}, ...
    {'bc2','bc3a','bc3b','bc4','bc5o','bc5i','bc5t','xbc','bc6'},curTis);

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

%% Altogether, we identified 587 synapses between RGCs innervated by the VG3 plexus and bipolar cells that innervate the VG3 plexus (Figure 9A).
% preVBs
% postVRs
% find syns from a to b



%% testing if the skels look okay.

for i=1:length(allSkels)
    curSkel=allSkels{i};
    curSM=curSkel.sm;
    curSWC=curSkel.swc;
    curPred=curSWC.pred(curSWC.arbor2swcID)+1;
    curPredInv=zeros(size(curPred));
    curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
    curSM.pred=curPred;
    allEdges=curSM.arbor.edges;
    uniqueNodes=unique(allEdges(:));
    nodeCounts=histcounts(allEdges,uniqueNodes);
    tipIDs=find(nodeCounts==1);
    forkIDs=find(nodeCounts>2);
    distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
    [srtd srtIdx] = sort(distFromRoot,'descend');
    tipIDsSrtd=tipIDs(srtIdx');
    usedNodes=[];
    usedEdges=[];
    branchNodes=cell(length(tipIDsSrtd),1);
    branchSize=zeros(length(tipIDsSrtd),1);
    for j=1:length(tipIDsSrtd)
        curTip=tipIDsSrtd(j);
        curBranch=curTip;
        curNode=curTip;
        nextEdge=0;
        stopBool=0;
        counter=0;
        while counter<5000 & curNode>0
            prevEdge=nextEdge(1);
            prevNode=curNode;
            usedNodes=[usedNodes;curNode];
            nextEdge=find(any(allEdges==curNode,2));
            nextEdge=nextEdge(nextEdge~=prevEdge);
            if isempty(nextEdge)
                break
            end
            edgeNodes=allEdges(nextEdge(1),:);
            curNode=edgeNodes(edgeNodes~=prevNode);
            curBranch=[curBranch;curNode];
            if ismember(curNode,usedNodes)
                break
            end
            counter=counter+1;
        end
        branchNodes{j}=curBranch;
        branchSize(j)=length(curBranch);
    end
    
    
    
    if 1
        figures{i}=figure();
        figStruct=compareMorph(figures{i},curSM.cid,highFV);
        title(curSM.cid)
        hold on
        for k=1:length(branchNodes)
            curBranch=branchNodes{k};
            branchPos=curSM.nep.pos(curBranch,:);
            figStruct.skel(k)=plot3(branchPos(:,3),branchPos(:,1),branchPos(:,2),'c-','LineWidth',1);
            %credge=allEdges(k,:);
            %plot3(curSM.nep.pos(credge,1),curSM.nep.pos(credge,2),curSM.nep.pos(credge,3));
            %drawnow();
        end
        figStruct.patches(1).FaceAlpha=0.1;
        figStruct.patches(1).FaceColor=[.7 0 .7];
        %figStruct.skel.LineWidth=3;
    end
end


%% make a skelgraph
for i=1:length(allSkels)
    curSM=allSkels{i}.sm;
    allEdges=curSM.nep.edges;
    uniqueNodes=curSM.nep.nodes;
    nodeCounts=tabulate(allEdges(:));
    tipIDs=nodeCounts(nodeCounts(:,2)==1);
    forkIDs=nodeCounts(nodeCounts(:,2)>2);
    
    rootNode=curSM.nep.seedNode;
    
    
end

%%

for i=1:length(allSkels)
    curSM=allSkels{i}.sm;
    curSMpreTypes=cid2type(curSM.syn.edges(:,2),highTis);
    curSMpostTypes=cid2type(curSM.syn.edges(:,1),highTis);
    allEdges=curSM.nep.edges;
    uniqueNodes=curSM.nep.nodes; %unique(allEdges(:));
    nodeCounts=tabulate(allEdges(:));
    tipIDs=nodeCounts(nodeCounts(:,2)==1);
    forkIDs=nodeCounts(nodeCounts(:,2)>2);
    rootID=intersect(tipIDs,find(curSM.nep.pos(:,1)==min(curSM.nep.pos(:,1))));
    distFromRoot=curSM.skel2skel.linDist(tipIDs,rootID);
    [srtd srtIdx] = sort(distFromRoot,'descend');
    tipIDsSrtd=tipIDs(srtIdx');
    [predList,nodeList]=getSkelGraph(allEdges,rootID,[],[rootID 0]);
    [srtd,idx]=sort(predList(:,1));
    preds=predList(idx,2);
    
    %     usedNodes=[];
    %     usedEdges=[];
    %     branchNodes=cell(length(tipIDsSrtd),1);
    %     branchSize=zeros(length(tipIDsSrtd),1);
    %     for j=1:length(tipIDsSrtd)
    %         curTip=tipIDsSrtd(j);
    %         curBranch=curTip;
    %         curNode=curTip;
    %         nextEdge=0;
    %         stopBool=0;
    %         counter=0;
    %         while counter<5000 & curNode>0
    %             prevEdge=nextEdge(1);
    %             prevNode=curNode;
    %             usedNodes=[usedNodes;curNode];
    %             nextEdge=find(any(allEdges==curNode,2));
    %             nextEdge=nextEdge(nextEdge~=prevEdge);
    %             if isempty(nextEdge)
    %                 break
    %             end
    %             edgeNodes=allEdges(nextEdge(1),:);
    %             curNode=edgeNodes(edgeNodes~=prevNode);
    %             curBranch=[curBranch;curNode];
    %             if ismember(curNode,usedNodes)
    %                 break
    %             end
    %             counter=counter+1;
    %         end
    %         branchNodes{j}=curBranch;
    %         branchSize(j)=length(curBranch);
    %     end
    %
    if 0
        figures{i}=figure();
        figStruct=compareMorph(figures{i},curSM.cid,highFV);
        figStruct.patches(1).FaceAlpha=0.1;
        figStruct.patches(1).FaceColor=[.7 0 .7];
        title(curSM.cid)
        hold on
        for j=1:length(tipIDsSrtd)
            nde=tipIDsSrtd(j);
            branchIDs=[];
            while nde~=rootID
                branchIDs=[branchIDs;nde];
                nde=preds(nde);
            end
            curp=plot3(curSM.nep.pos(branchIDs,3),curSM.nep.pos(branchIDs,1), ...
                curSM.nep.pos(branchIDs,2));
            curp.Color=rand(1,3);
            drawnow();
        end
    end
    if 0
        figures{i}=figure();
        figStruct=compareMorph(figures{i},curSM.cid,highFV);
        title(curSM.cid)
        hold on
        plotdNodes=[];
        for k=1:length(branchNodes)
            curBranch=branchNodes{k};
            branchPos=curSM.nep.pos(curBranch,:);
            figStruct.skel(k)=plot3(branchPos(:,3),branchPos(:,1),branchPos(:,2),'c-','LineWidth',1);
            nodes2plot=curBranch(~ismember(curBranch,plotdNodes));
            
            %credge=allEdges(k,:);
            scatter3(curSM.nep.pos(nodes2plot,3),curSM.nep.pos(nodes2plot,1),curSM.nep.pos(nodes2plot,2),10,'yo','filled');
            plotdNodes=[plotdNodes;nodes2plot];
            drawnow();
        end
        figStruct.patches(1).FaceAlpha=0.1;
        figStruct.patches(1).FaceColor=[.7 0 .7];
        %figStruct.skel.LineWidth=3;
    end
    
    %%
    % loop through the tips and get the synaptic environments at the last
    % micron
    synDists=[0 0 0];
    for y=1:length(tipIDs)
        curNode=tipIDs(y);
        curNode
        if curNode>0
            curTipEnv=[tipIDs(y) 0 0];
            nextNode=preds(curNode);
            curNodeLength=curSM.nep.props.nodeLength(curNode);
            curTipLength=curNodeLength;
            %curTipHistDat=[curNode,9,curNodeLength,curTipLength];
            beginningPos=curSM.arbor.nodes.pos(curNode,:);
            curNode=nextNode;
            while curTipLength<5
                if curNode==0
                    break
                end
                curNodeLength=curSM.nep.props.nodeLength(curNode);
                closeSyns=find(curSM.syn2Skel.closest==curNode);
                if ~isempty(closeSyns)
                    for h=1:length(closeSyns)
                        curSyn=closeSyns(h);
                        synDist=curSM.syn2Skel.syn2SkelDist(curSyn,curNode);
                        synDists=[synDists;[tipIDs(y),curSyn, ...
                            curTipLength+(curNodeLength/2)+synDist]];
                    end
                end
                curTipLength=curTipLength+curNodeLength;
                nextNode=preds(curNode);
                curNode=nextNode;
            end
            
            
        end
    end
    
    
    %% graph to check
    synDists=synDists(2:end,:);
    testFig3=figure();
    title(curSM.cid)
    hold on
    morphObj=compareMorph(testFig3,curSM.cid,highFV);
    morphObj.patches.FaceAlpha=0.05;
    morphObj.patches.FaceColor=[0.5 0.5 0];
    for curSynIt=1:size(synDists,1)
        [curTip,curSyn,curDist]=deal(synDists(curSynIt,1),synDists(curSynIt,2),synDists(curSynIt,3));
        curTipPos=curSM.nep.pos(curTip,:);
        curSynPos=curSM.syn.pos(curSyn,:);
        if curSMpreTypes{1}(curSyn)==7
            synCol=[0 1 1];
        else
            synCol=[1 0 1];
        end
        plot3([curTipPos(3) curSynPos(3)],[curTipPos(1) curSynPos(1)], ...
            [curTipPos(2) curSynPos(2)],'k:');
        scatter3(curSynPos(3),curSynPos(1),curSynPos(2),10,synCol,'v');
        
    end
end


%% The median distance between the VG3 input and output was 8.8 um with only 9% under 1 um (Figure 9E).
b2rMatrix=zeros(length(curTis.cells.type.subTypeNames{1})+1,length(curTis.cells.type.subTypeNames{7})+1);
allB2Rsyns=getClassConn({'bpc'},{'all'},{'rgc'},{'all'},curTis);
goodBPCtypes=[1:12,20,23];
for k=1:length(allB2Rsyns)
    curSyn=allB2Rsyns(k);
    b2rMatrix(alPost{3}(curSyn)+1,alPre{3}(curSyn)+1)=b2rMatrix(alPost{3}(curSyn)+1,alPre{3}(curSyn)+1)+1;
end
%figure(); image(b2rMatrix*5);
validFFLists={};
validFFListsText={};
validList=repmat("",1,1);
padRGCnams=["",curTis.cells.type.subTypeNames{1}];
padBPCnams=["",curTis.cells.type.subTypeNames{7}];
figure(); image(b2rMatrix*50);
xticks([1:length(padBPCnams)]); yticks([1:length(padRGCnams)]);
xticklabels(padBPCnams); yticklabels(padRGCnams);
b2rClipMat=b2rMatrix(2:end,2:end);
for w=1:length(padRGCnams)
    validFFLists{w}=intersect(find(b2rMatrix(w,:)>0)-1,goodBPCtypes);
    preRbpcTypes=find(b2rMatrix(w,:)>0)-1;
    validFFListsText{w}=padBPCnams(intersect(preRbpcTypes+1,goodBPCtypes+1));
    for x=1:length(intersect(find(b2rMatrix(w,:)>0),goodBPCtypes+1))
        validList=[validList;strjoin([padRGCnams(w) "_"  ...
            padBPCnams(preRbpcTypes(x)+1)])];
    end
end
validList=unique(validList);

bpcTypd=type2cid({'bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc','bpc'}, ...
    {'bc2','bc3a','bc3b','bc4','bc5o','bc5i','bc5t','xbc','bc6','bc7'},highTis);
bpcTypdCidList=[];
for q=1:length(bpcTypd)
    bpcTypdCidList=[bpcTypdCidList;bpcTypd{q}];
end
preVbpcTypdCidList=intersect(preVG3bpcCids,bpcTypdCidList);

nearestFFarray={};
allFFdists=[];

for p=1:length(allSkels)
    nearestFFdists=[];
    nearestFFdistsEuc=[];
    curSM=allSkels{p}.sm;
    curSMpreTypes=cid2type(curSM.syn.edges(:,2),highTis);
    curSMpostTypes=cid2type(curSM.syn.edges(:,1),highTis);
    fromKnown=ismember(curSM.syn.edges(:,2),preVbpcTypdCidList);
    toKnown=ismember(curSM.syn.edges(:,1),typedRGCcids)& ...
        ismember(curSM.syn.edges(:,2),allVG);
    for curSynIt=1:length(curSM.syn.edges(:,1))
        curEdges=curSM.syn.edges(curSynIt,:);
        if toKnown(curSynIt) %| fromKnown(curSynIt)
            connBtypes=validFFLists{curSMpostTypes{3}(curSynIt)+1};
            allowedBpcStrs=strjoin(padBPCnams(connBtypes+1));
            fprintf('%.0f: %s || %s\n',curEdges(1),curSMpostTypes{4}{curSynIt},allowedBpcStrs);
            validBinputs=fromKnown&ismember(curSMpreTypes{3},connBtypes)';
            validBindices=find(validBinputs);
            [a,nearestBpcInput]=min(curSM.syn2Skel.syn2SynDist(curSynIt,validBindices));
            nearestSyndex=validBindices(nearestBpcInput);
            %             fprintf('%.0f ; %.0f %.0f %.0f ; %.0f %.0f %.0f ; %.2f\n', curSynIt, ...
            %                 curSM.syn.pos(nearestSyndex,1),curSM.syn.pos(nearestSyndex,2),curSM.syn.pos(nearestSyndex,3), ...
            %                 curSM.syn.pos(curSynIt,1),curSM.syn.pos(curSynIt,2), ...
            %                 curSM.syn.pos(curSynIt,3),pdist2(curSM.syn.pos(nearestSyndex,:), ...
            %                 curSM.syn.pos(curSynIt,:)));
            nearestFFdistsEuc=[nearestFFdistsEuc;pdist2(curSM.syn.pos(nearestSyndex,:), ...
                curSM.syn.pos(curSynIt,:))];
            nearestFFdists=[nearestFFdists;curSM.syn2Skel.syn2SynDist(curSynIt,nearestSyndex)];
            
        end
    end
    nearestFFarray{p}=nearestFFdists;
    allFFdists=[allFFdists;nearestFFdists];
end

figure(); histogram(allFFdists,[0:1:30]);
median(allFFdists)

%% check with the looser definition

nearestFFarray={};
allFFdists=[];
for p=1:length(allSkels)
    nearestFFdists=[];
    curSM=allSkels{p}.sm;
    curSMpreTypes=cid2type(curSM.syn.edges(:,2),highTis);
    curSMpostTypes=cid2type(curSM.syn.edges(:,1),highTis);
    fromKnown=ismember(curSM.syn.edges(:,2),preVbpcTypdCidList);
    toKnown=ismember(curSM.syn.edges(:,1),postVGCrgcCids)& ...
        ismember(curSM.syn.edges(:,2),allVG);
    for curSynIt=1:length(curSM.syn.edges(:,1))
        curEdges=curSM.syn.edges(curSynIt,:);
        if toKnown(curSynIt) %| fromKnown(curSynIt)
            %connBtypes=validFFLists{curSMpostTypes{3}(curSynIt)+1};
            validBinputs=fromKnown;%&ismember(curSMpreTypes{3},connBtypes)';
            validBindices=find(validBinputs);
            [a,nearestBpcInput]=min(curSM.syn2Skel.syn2SynDist(curSynIt,validBinputs));
            [b,nearestBpcInputEuc]=min(pdist2(curSM.syn.pos(curSynIt,:), ...
                curSM.syn.pos(validBinputs,:)));
            if a<b
                nearestSyndex=validBindices(nearestBpcInput);
            elseif b<a
                nearestSyndex=validBindices(nearestBpcInputEuc);
            else
                nearestSyndex=validBindices(nearestBpcInput);
            end
            %             fprintf('%.0f ; %.0f %.0f %.0f ; %.0f %.0f %.0f ; %.2f\n', curSynIt, ...
            %                 curSM.syn.pos(nearestSyndex,1),curSM.syn.pos(nearestSyndex,2),curSM.syn.pos(nearestSyndex,3), ...
            %                 curSM.syn.pos(curSynIt,1),curSM.syn.pos(curSynIt,2), ...
            %                 curSM.syn.pos(curSynIt,3),pdist2(curSM.syn.pos(nearestSyndex,:), ...
            %                 curSM.syn.pos(curSynIt,:)));
            nearestFFdists=[nearestFFdists;pdist2(curSM.syn.pos(nearestSyndex,:), ...
                curSM.syn.pos(curSynIt,:))];
        end
    end
    nearestFFarray{p}=nearestFFdists;
    allFFdists=[allFFdists;nearestFFdists];
end

figure(); histogram(allFFdists,[0:1:30]);

%% For one micrometer-long samples of neurite tips, only 24% of 1853
% neurites had any synapses and only 0.3% of the neurites had bipolar
% cell inputs,RGC outputs, and no amacrine cell inputs. Examining longer
% lengths of neurite (5, 10 um) revealed more combinations of multiple
% synapse types, but amacrine cell input-free lengths that included
% bipolar cells and RGCs remained rare (Figure 4D).

%first, find the closest tip to each of the synapses.
for i=1:length(allSkels)
    curSM=allSkels{i}.sm;
    curSMpreTypes=cid2type(curSM.syn.edges(:,2),highTis);
    curSMpostTypes=cid2type(curSM.syn.edges(:,1),highTis);
    allSynEdges=curSM.syn.edges;
    allEdges=curSM.nep.edges;
    uniqueNodes=curSM.nep.nodes; %unique(allEdges(:));
    nodeCounts=tabulate(allEdges(:));
    tipIDs=nodeCounts(nodeCounts(:,2)==1);
    forkIDs=nodeCounts(nodeCounts(:,2)>2);
    closestTip=zeros(length(curSM.syn.synID),1);
    closestTipDist=zeros(length(curSM.syn.synID),1);
    for j=1:length(curSM.syn.synID)
        %curEdges=allEdges(j,:);
        %curPreType=[curSMpreTypes{1}(j),curSMpreTypes{3}(j)];
        %curPostType=[curSMpostTypes{1}(j),curSMpostTypes{3}(j)];
        [a,closestTipInd]=min(curSM.syn2Skel.syn2SkelDist(j,tipIDs));
        closestTip(j)=tipIDs(closestTipInd);
        closestTipDist(j)=a;
    end
    curSM.syn.closestTip=closestTip;
    curSM.syn.closestTipDist=closestTipDist;
    allSkels{i}.sm.syn.closestTip=closestTip;
    allSkels{i}.sm.syn.closestTipDist=closestTipDist;
    if 0
        testFig4=figure();
    title(curSM.cid)
    hold on
    morphObj=compareMorph(testFig4,curSM.cid,highFV);
    morphObj.patches.FaceAlpha=0.05;
    morphObj.patches.FaceColor=[0.5 0.5 0];
    for curSynIt=1:size(closestTip,1)
        curTipPos=curSM.nep.pos(closestTip(curSynIt,:),:);
        curSynPos=curSM.syn.pos(curSynIt,:);
        if allSynEdges(curSynIt,1)==curSM.cid
            synCol=[0 1 1];
            synMarker='v';
        elseif allSynEdges(curSynIt,2)==curSM.cid
            synCol=[1 0 1];
            synMarker='^';
        else
            synCol=[1 1 0];
            synMarker='o';
        end
        plot3([curTipPos(3) curSynPos(3)],[curTipPos(1) curSynPos(1)], ...
            [curTipPos(2) curSynPos(2)],'k:');
        scatter3(curSynPos(3),curSynPos(1),curSynPos(2),10,synCol,synMarker);
        
    end
    end
    
end

% continuing, get the closest synapses for each tip.
tipResults=cell(length(allSkels),1);
fullResults={};
totalLocalSynTips=[];
totalLocalSynTipsFF=[];
for i=1:length(allSkels)
    curTipResults=[0 0 0];
    curSM=allSkels{i}.sm;
    curSMpreTypes=cid2type(curSM.syn.edges(:,2),highTis);
    curSMpostTypes=cid2type(curSM.syn.edges(:,1),highTis);
    allSynEdges=curSM.syn.edges;
    inputs=find(allSynEdges(:,1)==curSM.cid);
    outputs=find(allSynEdges(:,2)==curSM.cid);
    binputs=intersect(inputs,find(curSMpreTypes{1}==7));
    routputs=intersect(outputs,find(curSMpostTypes{1}==1));
    synCode=zeros(size(curSM.syn.edges(:,1)));
    synCode(inputs)=3;
    synCode(outputs)=4;
    synCode(binputs)=1;
    synCode(routputs)=2;
    allEdges=curSM.nep.edges;
    uniqueNodes=curSM.nep.nodes; %unique(allEdges(:));
    nodeCounts=tabulate(allEdges(:));
    tipIDs=nodeCounts(nodeCounts(:,2)==1);
    nearbySyns=cell(size(tipIDs));
    nearSynDists=cell(size(tipIDs));
    nearSynCode=cell(size(tipIDs));
    forkIDs=nodeCounts(nodeCounts(:,2)>2);
    synLocalMicron=zeros(size(tipIDs));
    synLocalMicronFF=zeros(size(tipIDs));
    for k=1:length(tipIDs)
        curTip=tipIDs(k);
        closeSynIDs=find(curSM.syn.closestTip==curTip);
        if ~isempty(closeSynIDs)
            nearbySyns{k}=closeSynIDs;
            nearSynDists{k}=curSM.syn.closestTipDist(closeSynIDs);
            nearSynCode{k}=synCode(closeSynIDs);
            if min(curSM.syn.closestTipDist(closeSynIDs))<1
                synLocalMicron(k)=1;
                if ismember(1,synCode(closeSynIDs))&ismember(2,synCode(closeSynIDs))& ...
                        ~ismember(3,synCode(closeSynIDs))
                    synLocalMicronFF(k)=1;
                end
            end
            curt=horzcat(closeSynIDs,curSM.syn.closestTipDist(closeSynIDs),synCode(closeSynIDs));
            curTipResults=[curTipResults;curt];
        end
    end
    tipResults{i}=curTipResults;
    fullResults{i}={nearbySyns,nearSynDists,nearSynCode,curTipResults};
    totalLocalSynTips=[totalLocalSynTips;synLocalMicron];
    totalLocalSynTipsFF=[totalLocalSynTipsFF;synLocalMicronFF];
end

fprintf('For one micrometer-long samples of neurite tips, only %.1f%% of %.0f neurites had any synapses\n', ...
    mean(totalLocalSynTips)*100,length(totalLocalSynTips));

fprintf('and only %.2f%% had bpc input, rgc output, and no amacrine input\n', ...
    mean(totalLocalSynTipsFF)*100);

%% test distance calculation

curSM=allSkels{1}.sm;
distTestIDs=[357 350 349 348];
curSM.syn2Skel.syn2SkelDist(11,distTestIDs);
posMat=curSM.nep.pos(distTestIDs,:);
synPos=curSM.syn.pos(11,:);
posMat=[synPos;posMat];
testFig4=figure();
title(curSM.cid)
hold on
morphObj=compareMorph(testFig5,curSM.cid,highFV);
morphObj.patches.FaceAlpha=0.05;
morphObj.patches.FaceColor=[0.5 0.5 0];

scatter3(posMat(1,3),posMat(1,1),posMat(1,2),50,'cv');
scatter3(posMat(2:end,3),posMat(2:end,1),posMat(2:end,2),50,'mo','filled');

plot3(posMat(:,3),posMat(:,1),posMat(:,2),'k:');

avgLoc=mean(posMat,1);

vertBox=[avgLoc-0.2;avgLoc+0.2];
goodVerts=find(curSM.nep.fv.vertices(:,1)>vertBox(1,1)&curSM.nep.fv.vertices(:,1)<vertBox(2,1)& ...
    curSM.nep.fv.vertices(:,2)>vertBox(1,2)&curSM.nep.fv.vertices(:,2)<vertBox(2,2)& ...
    curSM.nep.fv.vertices(:,3)>vertBox(1,3)&curSM.nep.fv.vertices(:,3)<vertBox(2,3));

goodVerts2=find(pdist2(avgLoc,curSM.nep.fv.vertices)'<0.3);

scatter3(curSM.nep.fv.vertices(goodVerts,3),curSM.nep.fv.vertices(goodVerts,1), ...
    curSM.nep.fv.vertices(goodVerts,2),50,'y.');

xlim([avgLoc(3)-0.2 avgLoc(3)+0.2]);
ylim([avgLoc(1)-0.2 avgLoc(1)+0.2]);
zlim([avgLoc(2)-0.2 avgLoc(2)+0.2]);

%% trying with another synapse
synID=26; nodeList=curSM.syn.closestTip(synID);
closeTip=curSM.syn.closestTip(synID);
closeTipDist=curSM.syn.closestTipDist(synID);
curNode=nodeList;
closestNode=curSM.syn2Skel.closest(synID);
nextNode=0;
while nextNode~=2340
    nextNode=curSM.nep.pred(curNode);
    nodeList=[nodeList;nextNode];
    curNode=nextNode;
end

nodeListInv=nodeList(length(nodeList):-1:1);

curSM=allSkels{1}.sm;
distTestIDs=nodeListInv;
curSM.syn2Skel.syn2SkelDist(synID,distTestIDs);
posMat=curSM.nep.pos(distTestIDs,:);
synPos=curSM.syn.pos(synID,:);
posMat=[synPos;posMat];

testDist=pdist2(posMat(1,:),posMat(end,:))

ran=range(posMat);
ranRad=(ran/2)+0.2;
testFig6=figure();
title(curSM.cid)
hold on
morphObj=compareMorph(testFig6,curSM.cid,highFV);
morphObj.patches.FaceAlpha=0.05;
morphObj.patches.FaceColor=[0.5 0.5 0];

scatter3(posMat(1,3),posMat(1,1),posMat(1,2),100,'cv','filled');
scatter3(posMat(2:end,3),posMat(2:end,1),posMat(2:end,2),50,'mo','filled');

plot3(posMat(:,3),posMat(:,1),posMat(:,2),'k:');

avgLoc=mean(posMat,1);

vertBox=[avgLoc-ranRad;avgLoc+ranRad];
goodVerts=find(curSM.nep.fv.vertices(:,1)>vertBox(1,1)&curSM.nep.fv.vertices(:,1)<vertBox(2,1)& ...
    curSM.nep.fv.vertices(:,2)>vertBox(1,2)&curSM.nep.fv.vertices(:,2)<vertBox(2,2)& ...
    curSM.nep.fv.vertices(:,3)>vertBox(1,3)&curSM.nep.fv.vertices(:,3)<vertBox(2,3));

goodVerts2=find(pdist2(avgLoc,curSM.nep.fv.vertices)'<max(ranRad));

scatter3(curSM.nep.fv.vertices(goodVerts,3),curSM.nep.fv.vertices(goodVerts,1), ...
    curSM.nep.fv.vertices(goodVerts,2),50,'y.');

xlim([avgLoc(3)-ranRad(3) avgLoc(3)+ranRad(3)]);
ylim([avgLoc(1)-ranRad(1) avgLoc(1)+ranRad(1)]);
zlim([avgLoc(2)-ranRad(2) avgLoc(2)+ranRad(2)]);
