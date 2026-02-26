
global tis
vgcCidList=[2 3 4 5 10 11 13 14 20];
UScids=tis.syn.edges(ismember(tis.syn.edges(:,1),vgcCidList),2);
UStypes=cid2type(UScids,tis);
%make the first figure
f=figure();
hold on;
histogram(UStypes{1},length(tis.cells.type.typeNames)); 
xticks([0:length(tis.cells.type.typeNames)]+0.5);
typeNames=tis.cells.type.typeNames(:);
typeNames2=convertCharsToStrings(typeNames);
xticklabels({'none' typeNames{:}})
title('# of input synapses by cell type');



BPCids=find(UStypes{1}==7);
BPCcids=UScids(BPCids);
BPCtypes{1}=UStypes{1}(BPCids);
BPCtypes{2}=UStypes{2}(BPCids);
BPCtypes{3}=UStypes{3}(BPCids);
BPCtypes{4}=UStypes{4}(BPCids);

f=figure();
hold on;
histogram(BPCtypes{3},length(tis.cells.type.subTypeNames{7})); 
xticks([0:length(tis.cells.type.subTypeNames{7})]+0.5);
typeNames=tis.cells.type.subTypeNames{7};
typeNames2=convertCharsToStrings(typeNames);
xticklabels({'none' typeNames{:}})
title('# of VGC inputs by BPC subtype');

% Get the cids for these different things
TOIlist={'none', 'bcunk', 'bcon', 'bcoff'};
for k=1:length(TOIlist)
    curSubType=TOIlist{k};
    if strcmp(curSubType,'none')
        curSubTypeID=0;
    else
        curSubTypeID=find(contains(tis.cells.type.subTypeNames{7},curSubType)==1);
    end
    curSubType
    curSubCidList=BPCcids(BPCtypes{3}==curSubTypeID)
    
end


