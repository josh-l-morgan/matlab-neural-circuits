function[col] = driftNetCol(nep,chng);

if ~exist('chng','var') 
    chng = 1;
end

seedNode = nep.seedNode;
edges = nep.edges;
pos = nep.pos;
coli = nep.nodes*0;
coli(seedNode) = 1;

checkNodes = seedNode;
hold on
set(gca,'clipping','off')
while sum(coli==0)

    nextCheck = [];
    for n = 1:length(checkNodes)
        cn = checkNodes(n);
        cCol = coli(cn);
        hitN = cat(1,edges(edges(:,1) == cn,2), edges(edges(:,2) == cn,1));
        openN = hitN(coli(hitN)==0); %not colored
        %openN = hitN;
        newCol = cCol + rand(length(openN),1)*chng-chng/2;
        coli(openN) = newCol;
        nextCheck = cat(1,nextCheck,openN);
    end
    checkNodes = nextCheck;

end
hold off


cNum = 1000;
cmap = hsv(cNum);
coli = mod(round(coli*cNum),cNum)+1;

col = cmap(coli,:);







