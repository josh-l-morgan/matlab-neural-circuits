function [outputDat stats]=getDRG(cidCellArrayA,cidCellArrayB,SMcellArray,binEdges,binWidth)
%% Just brainstorming for a function to handle this stuff.
% take in a cell array of SMs fr the relevant vgcs
% take in two cid lists to form the axes of the matrix for the connectivity
% so, like, you could have one axes be lists for all the ON and OFF bpcs
% and then the other axis could be some different types of RGCs
% The function will:
% Go through each of the supplied SM files,
% Make a matrix of the size m*n for the input lists
% Find the IDs of the synapses (OUTPUT)
% Find the distances for each of the matrix boxes (OUTPUT)
% generate the relative influences (OUTPUT)
% generate the density recovery profile (OUTPUT)
synRad=binEdges(end);
% copied stuff from the DRG script
distBin=binWidth;
g2g = zeros(29,length(cidCellArrayA),length(binEdges));
nearLength = g2g;
countG = zeros(2,1);
%outputDat=zeros(length(cidCellArrayA),length(cidCellArrayB),length(SMcellArray));
stats=zeros(length(cidCellArrayB),1);
testBins=binEdges;
outputDat=struct();
%% 
for curSMit=1:length(SMcellArray)
    curSM=SMcellArray{curSMit};
    edgeCids=curSM.syn.edges(:,[1 2]);
    %get the max differences of the distances between synapses
    nodeLengths=curSM.nep.props.nodeLength;
    pos = curSM.syn.pos;
    numSyn = size(pos,1);
    dif1 = pos(:,1)-pos(:,1)';
    dif2 = pos(:,2)-pos(:,2)';
    dif3 = pos(:,3)-pos(:,3)';
    syn2synEuc = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
    syn2synLin = curSM.syn2Skel.syn2SynDist;
    syn2synMax = max(cat(3,syn2synEuc,syn2synLin),[],3);
    curPreTypes=cid2type(curSM.syn.edges(:,2),curTis);
    curPostTypes=cid2type(curSM.syn.edges(:,1),curTis);
    curInputs=find(curSM.syn.edges(:,1)==curCid);
    curOutputs=find(curSM.syn.edges(:,2)==curCid);
    closeSyns=find(dist2synsComb<synRad);
    nearInputs=intersect(closeSyns,curInputs);
    nearOutputs=intersect(closeSyns,curOutputs);
    %go through the close synapses and make the density of the
    %different types
    bpcInIDs=intersect(nearInputs,find(curPreTypes{1}==7));
    amcOutIDs=intersect(nearOutputs,find(curPostTypes{1}==0|curPostTypes{1}==8));
    rgcOutSynIDs=intersect(nearOutputs,find(curPostTypes{1}==1));
    rgcFFSynIDs=intersect(rgcOutSynIDs,find(ismember(curPostCids,rgcCid)));
    rgcOutSynIDs=rgcOutSynIDs(~ismember(rgcOutSynIDs,rgcFFSynIDs));
    amcInIDs=intersect(nearInputs,find(curPreTypes{1}==0|curPreTypes{1}==8));
    synIDcell={bpcInIDs,amcInIDs,rgcOutSynIDs,rgcFFSynIDs,amcOutIDs};
    
    tots=tots+[length(bpcInIDs),length(amcInIDs),length(rgcOutSynIDs),length(rgcFFSynIDs),length(amcOutIDs)];
    for itA=1:length(cidCellArrayA)
        curCidsA=cidCellArrayA{itA};
        listAIDs=find(ismember(edgeCids(:,1),curCidsA)|ismember(edgeCids(:,2),curCidsA));
        outputDat(curSMit).locsA{itA}=curSM.syn.pos(listAIDs,:);
        outputDat(curSMit).edgesA{itA}=curSM.syn.edges(listAIDs,:);
        outputDat(curSMit).idsA{itA}=listAIDs;
        for itB=1:length(cidCellArrayB)
            curCidsB=cidCellArrayB{itB};
            listBIDs=find(ismember(edgeCids(:,1),curCidsB)|ismember(edgeCids(:,2),curCidsB));
            [g,i,listBdepths]=getIPLdepth(curSM.syn.pos(listBIDs,3),curSM.syn.pos(listBIDs,1),curSM.syn.pos(listBIDs,2),[],[]);
            %listBIDs=listBIDs(listBdepths>.47);
            outputDat(curSMit).locsB{itB}=curSM.syn.pos(listBIDs,:);
            outputDat(curSMit).edgesB{itB}=curSM.syn.edges(listBIDs,:);
            outputDat(curSMit).idsB{itB}=listBIDs;
            for curSynIt=1:length(listBIDs)
                curSyn=listBIDs(curSynIt);
                curSynPos=curSM.syn.pos(curSyn,:);
                %curCloseNode=curSM.syn2Skel.closest(curSyn);
                dist2syns=syn2synMax(curSyn,:);
                dist2syns(curSyn) = inf;
                dist2skels=curSM.syn2Skel.syn2SkelDist(curSyn,:);
                for it2A=1:length(cidCellArrayA)
                    curCids2A=cidCellArrayA{it2A};
                    listAIDs2=find(ismember(edgeCids(:,1),curCids2A)|ismember(edgeCids(:,2),curCids2A));
                    countG(it2A)=countG(it2A)+1;
                    %graphDenom(curSynIt,it2A)=length(synIDcell{it2A});
                    synDG2 = dist2syns(listAIDs2);
                    for b = 1:(length(binEdges))
                        %get the indices and items for everything within a binwidth of that range
                        [ y x ] = find((synDG2>=(binEdges(b)-distBin/2)) & (synDG2<(binEdges(b)+distBin/2)));
                        %add the number of things at that range to the data
                        g2g(curSynIt,it2A,b) =  g2g(curSynIt,it2A,b) + length(x);
                        %g2g2(curSynIt,it2A,b) = g2g2(curSynIt,it2A,b) + length(x); %/length(synIDcell{hT}) ;
                    end
                end
                for b = 1:(length(binEdges)) %get length for density
                    %find all the skeleton nodes at that range from the point
                    [ y x] = find((dist2skels>=(binEdges(b)-distBin/2)) & (dist2skels<(binEdges(b)+distBin/2)));
                    %that is the denominator for the density!
                    %how many synapses at that range / sum of the nodes at that range
                    nearLength(it2A,:,b) = nearLength(it2A,:,b) + sum(nodeLengths(x));
                end
            end
            distMat=curSM.syn2Skel.syn2SynDist(listAIDs,listBIDs);
            closeDistMat=distMat(distMat<binEdges(end));
            histDat=histcounts(closeDistMat,testBins);
            outputDat(curSMit).histDat(:,itA,itB)=histDat;
            stats(itB)=stats(itB)+length(listBIDs);
            
        end
    end
end
