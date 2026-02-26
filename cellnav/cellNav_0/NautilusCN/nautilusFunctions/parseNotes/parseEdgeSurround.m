function[conProp] = parseEdgeSurround(nam,edgeTag);

postID = [];
preID = [];
edge1 = [];
preGroup = [];
postGroup = [];

s = regexp(nam,edgeTag);

if ~isempty(s)
    
    namPre = nam(1:s(1)-1);
    pF = regexp(namPre,'(');
    pB = regexp(namPre,')');
    if ~isempty(pF) & ~isempty(pB)
        pF = pF(end);
        pB = pB(end);
    end
    
    [cids cidPos] = getCids(namPre);
    if pB>pF
        inCid = cids((cidPos>pF) & (cidPos<pB));
    else
        inCid = [];
    end
    lastCid = cids(cidPos==max(cidPos));
    if sum(pB>max(cidPos)) & ~isempty(inCid)
        preID = inCid;
    else
        preID = lastCid;
    end
    
    namPost = nam(s(1)+2:end);
    pF = regexp(namPost,'(');
    pB = regexp(namPost,')');
    if ~isempty(pF) & ~isempty(pB)
        pF = pF(end);
        pB = pB(end);
    end
    
    [cids cidPos] = getCids(namPost);
    if pB>pF
        inCid = cids((cidPos>pF) & (cidPos<pB));
    else
        inCid = [];
    end
    firstCid = cids(cidPos==min(cidPos));
    if sum(pF<max(cidPos)) & ~isempty(inCid)
        postID = inCid;
    else
        postID = firstCid;
    end
    
    %% make edge
    c = 0;
    edge1 = [];
    for y = 1:length(preID)
        for x = 1:length(postID)
            c = c+1;
            edge1(c,:) = [preID(y) postID(x)];
        end
    end
    
    if length(preID)>1
        preGroup = nchoosek(preID,2);
    end
    
    if length(postID)>1
        postGroup = nchoosek(postID,2);
    end
    
else
    
    %disp(sprintf('The tag "%s" was not found in "%s".',edgeTag,nam));
end


conProp.name = nam;
conProp.tag = edgeTag;
conProp.preID = preID;
conProp.postID = postID;
conProp.preGroup = preGroup;
conProp.postGroup = postGroup;
conProp.edges = edge1;




