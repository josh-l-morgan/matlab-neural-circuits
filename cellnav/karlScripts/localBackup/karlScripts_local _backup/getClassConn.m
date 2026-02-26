function output=getClassConn(preType,preSub,postType,postSub,curTis)
%get the synapses connecting the cell types.
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);


if isnumeric(preType)
    if preSub==99
        preSub=unique(allPreTypes{3}(allPreTypes{1}==preType));
    end
    preMatch=find(ismember(allPreTypes{1},preType)&ismember(allPreTypes{3},preSub));
else
    if ismember(preSub,{'all'})
        preSub=unique(allPreTypes{4}(ismember(allPreTypes{2},preType)));
    end
    preMatch=find(ismember(allPreTypes{2},preType)&ismember(allPreTypes{4},preSub));
end

if isnumeric(postType)
    if postSub==99
        postSub=unique(allPostTypes{3}(allPostTypes{1}==postType));
    end
    postMatch=find(ismember(allPostTypes{1},postType)&ismember(allPostTypes{3},postSub));
else
    if ismember(postSub,{'all'})
        postSub=unique(allPostTypes{4}(ismember(allPostTypes{2},postType)));
    end
    postMatch=find(ismember(allPostTypes{2},postType)&ismember(allPostTypes{4},postSub));
end

matchIDs=intersect(preMatch,postMatch);
output=matchIDs;

end