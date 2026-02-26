function [outNum outIDs]=getVGinputs(cids,tis)
allVGCcids=[2 3 4 5 10 11 12 13 14 20];
outIDs=cell(length(cids),1);
outNum=zeros(length(cids),1);
for i=1:length(cids)
    cid=cids(i);
    synIDs=find(tis.syn.edges(:,1)==cid & ismember(tis.syn.edges(:,2),allVGCcids));
    outNum(i)=length(synIDs);
    outIDs{i}=synIDs;


end