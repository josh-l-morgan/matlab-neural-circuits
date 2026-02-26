function[colCon] = graphCon(syn,typeLab)



[sortType typePos] = sort(typeLab);

synTypes = typeLab(syn);
useCells = unique(syn(:));
cellTypes = typeLab(useCells);

[sortCells cellIdx] = sort(cellTypes); %use fopr types
sortIDs = useCells(cellIdx);

lookUp = zeros(1,max(syn(:)));
lookUp(sortIDs) = 1:length(sortIDs);

sortSyn = lookUp(syn); %new sorted synapse matrix
numSyn = length(sortSyn);
numCell = length(useCells);
cells = 1:numCell;
typeKey = [' RGC,       Inhibitory,       thalamocortical,       glia,       other'];


%%  Show connections

conG = zeros(numCell,numCell);
conG = conG + repmat(mod(sortCells+1,2),[length(sortCells) 1]);
conG = conG + repmat(mod(sortCells+1,2),[length(sortCells) 1]);
conG = xor(conG,repmat(mod(sortCells'+1,2),[1 length(sortCells)]));

conG = conG ;

conS = conG*0;

colCon = zeros(numCell,numCell,3,'uint8');

for i = 1:size(sortSyn,1)

    if ~sum(sortSyn(i,:)==0) & abs(diff(sortSyn(i,:)))
        prePos = find(typePos == sortSyn(i,1));
        postPos = find(typePos == sortSyn(i,2));
        conS(sortSyn(i,1),sortSyn(i,2)) = conG(sortSyn(i,1),sortSyn(i,2)) +1;
        %image(conS)
    end
end

conG = repmat(conG,[1 1 3]);

%conS = repmat((conS)*100,[1 1 3]);

colCon = uint8(conG*20);
colCon(:,:,3) = colCon(:,:,3) + uint8(conS*300);
background = colCon;
 
image(colCon)
title(typeKey)