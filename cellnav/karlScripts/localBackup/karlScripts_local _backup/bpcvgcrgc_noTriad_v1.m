%% Is the shared ribbon thing real?

% get all the synapses where there are bpc inputs to vgc and to rgc

% get the ones that are less than 0.1um that have the same upstream bpc
% Check the distances and see if there is a pattern there

if ~exist('curTis')
    global glob tis
    curTis=tis;
end

cleanTis=findBadSynapse(curTis);

curTis=cleanTis;

% get the type information for all the synapses
synTypePre=cid2type(curTis.syn.edges(:,2),curTis);
synTypePost=cid2type(curTis.syn.edges(:,1),curTis);
allSynPos=curTis.syn.pos;
%get the distances from pos to pos
synDists=pdist2(curTis.syn.pos,curTis.syn.pos);
distMat=pdist2(allSynPos,allSynPos);
noDistSyns=find(distMat<0.1);
self=sub2ind(size(distMat),[1:size(distMat,1)],[1:size(distMat,2)]);
noDistSynsClean=noDistSyns(~ismember(noDistSyns,self));
[badx,bady]=ind2sub(size(distMat),noDistSynsClean);
coLocInd=horzcat(badx,bady);
[idx,idx]=unique(sort(coLocInd')','rows','stable');
coLocInd=coLocInd(idx,:);

rib2vgcCount=0; rib2VGCOFFa=0; rib2VGCW3=0; rib2VGCRGC=0;
rgcTargSubList=[];
rgcCoTargCidList=[];
%go through the 
for coLocIt=1:length(coLocInd(:,1))
    curPair=coLocInd(coLocIt,:);
    preCids=curTis.syn.edges(curPair,2);
    postCids=curTis.syn.edges(curPair,1);
    if preCids(1)==preCids(2) && curTis.syn.synType(curPair(1))==2 ...
        && synTypePre{1}(curPair(1))==7
        %pair with same loc and preSyn BPC + syntype==2
        postAtype=[synTypePost{1}(curPair(1)) synTypePost{3}(curPair(1))];
        postBtype=[synTypePost{1}(curPair(2)) synTypePost{3}(curPair(2))];
        postTypes=vertcat(postAtype,postBtype);
        if sum(postBtype==[8 1])==2 | sum(postAtype==[8 1])==2
            %one ribTarg is a VGC
            if postAtype(1)==1 | postBtype(1)==1
                rib2VGCRGC=rib2VGCRGC+1;
                rgcTargSub=postTypes(postTypes(:,1)==1,2);
                rgcTargSubList=[rgcTargSubList;rgcTargSub];
                rgcCoTargCidList=[rgcCoTargCidList;postCids(postTypes(:,1)==1)];
            end
            if sum(postBtype==[1 23])==2 | sum(postAtype==[1 23])==2
                rib2VGCOFFa=rib2VGCOFFa+1;
                curTis.syn.pos(curPair(1),:)
                curTis.syn.pos(curPair(2),:)
            end
            if sum(postBtype==[1 24])==2 | sum(postAtype==[1 24])==2
                rib2VGCW3=rib2VGCW3+1;
                
            end
            %[postAtype postBtype]
            rib2vgcCount=rib2vgcCount+1;
        end
    end
end

rgcSubHistDat=histcounts(rgcTargSubList,-0.5:1:59.5);
rgcNameLabels=curTis.cells.type.subTypeNames{1};
rgcNameLabels={'none',rgcNameLabels{:}};
h1=figure();
bar(rgcSubHistDat);
xticks([1:length(rgcNameLabels)]);
xticklabels(rgcNameLabels);

rgcTypeColMap=turbo(24);
rgcRandColMap=rgcTypeColMap(randperm(length(rgcTypeColMap)),:);

rgcCoTargs=unique(rgcCoTargCidList);
rgcCoTargTypeDat=cid2type(rgcCoTargs,curTis);
rgcCoTargSubs=rgcCoTargTypeDat{3};
uniqueTargSubList=unique(rgcCoTargSubs,'stable');
idx = arrayfun(@(x)find(uniqueTargSubList==x,1),rgcCoTargSubs);
rgcBarCol=rgcRandColMap(idx,:);
rgcBarColUn=unique(rgcBarCol,'rows','stable');
legendSubNames=curTis.cells.type.subTypeNames{1}(uniqueTargSubList(uniqueTargSubList~=0));

% Change to a stacked bar or something so that we can see the VGC->RGC
% BPC->RGC and BPC->RGC all at the same time.
%VGC - W3 synapses should be able to share BPC inputs (they are in the
%right places), but they do not. What is driving the avoidance of dyadic
%inputs from OFFBPC->VGC/W3?

% go back through the loop and make it so that you go through each of the
% subtypes and plot out those bars and do the colors on the spot and add
% the name of the subtype to the list.

rgcCidHistDat=histcounts(rgcCoTargCidList,[0;rgcCoTargs-0.5]);
rgcCidLabels=cellstr(num2str(unique(rgcCoTargCidList)));
h2=figure();
for i=1:length(rgcBarColUn)
    scatter(0,0,1,rgcBarColUn(i,:),'filled');
    hold on
end
barg=bar(rgcCidHistDat,'facecolor','flat');
hold on
barg.CData=rgcBarCol;
xticks([2:length(rgcCidLabels)]);
xticklabels(rgcCidLabels);
legend(legendSubNames);
