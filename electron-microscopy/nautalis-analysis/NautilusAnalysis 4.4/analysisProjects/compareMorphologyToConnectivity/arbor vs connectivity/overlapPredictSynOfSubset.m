function[predSyn nodeFrac] = overlapPredictSynOfSubset(varargin);

%%make overlap based prediction using mask, requires overlap structure from skelOverlap or filteredSkelOverlap




skelOverlapPred = varargin{1};
binRes = varargin{2};
testMask = varargin{3};
useShift = 0;

for v = 1:length(varargin); 
    
   varin = varargin{v};
   if ischar(varin)
   switch varin
       case 'shift';
           useShift = varargin{v+1};
   end
   end
end


axCellHists = skelOverlapPred.axCellHists; % Get real overlaps


%Get shifted overlaps if shift is used
if ~useShift
    useOverlapAxCellHists = axCellHists;
else
    for y = 1:size(skelOverlapPred.shiftCellHists,1)
        for x = 1:size(skelOverlapPred.shiftCellHists,2)
            useOverlapAxCellHists{y,x} = skelOverlapPred.shiftCellHists{y,x,useShift};
        end
    end
end



conDists = skelOverlapPred.conDists;
con = skelOverlapPred.con;
histBin = skelOverlapPred.histBin;

useOverlap = testMask;
sumHists = zeros(1,length(histBin));
synDists = [];
for i = 1:size(conDists,1)  %% Get synapse distances
    for p = 1:size(conDists,2)
        if  useOverlap(i,p);  % only use traced post synaptic cells
            synDists = cat(1,synDists,  conDists{i,p});
            sumHists = sumHists + axCellHists{i,p}';
        end
    end
end

synHistDist = hist(synDists,histBin);

newBins = [0:binRes:100];
clear nodeCount synCount
for i =1 :length(newBins)-1
    groupBin = find((histBin>=newBins(i)) & (histBin<newBins(i+1)));
    nodeCount = sum(sumHists(groupBin));
    synCount = sum((synDists>=newBins(i)) & (synDists<newBins(i+1)));
    hFrac = synCount/nodeCount;
    nodeFrac(groupBin) = hFrac;
end


nodeFrac(isnan(nodeFrac)) = 0;

predSyn = con * 0;
for i = 1:size(conDists,1)
    for p = 1:size(conDists,2)
        pairCount= useOverlapAxCellHists{i,p} ;
        predSyn(i,p) = sum(pairCount(1:length(nodeFrac))'.*nodeFrac);
    end
end

