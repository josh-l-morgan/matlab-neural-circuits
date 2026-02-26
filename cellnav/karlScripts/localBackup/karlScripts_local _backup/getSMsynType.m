function typeInfo=getSMsynType(synID,curCid,curSM,tisDat)
synType=5;
curEdges=curSM.syn.edges(synID,:);
edgeTypes=cid2type(curEdges(1:2),tisDat);
if curEdges(1)==curCid %input
    if edgeTypes{1}(2)==7
        synType=1; %BPC IN
    elseif curEdges(2)==0 | edgeTypes{1}(2)==8
        synType=2; %AMC IN
    end
elseif curEdges(2)==curCid %output
    if edgeTypes{1}(1)==1
        synType=3; %RGC OUT
    elseif curEdges(1)==0 | edgeTypes{1}(1)==8
        synType=4; %AMC OUT
    end
end
typeInfo=synType;
end
