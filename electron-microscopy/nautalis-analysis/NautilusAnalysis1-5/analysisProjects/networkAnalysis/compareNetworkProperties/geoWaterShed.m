


loadData = 0;
if loadData
    clear all
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    
    seedList = [ 108  201 109 907 903];
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
end


%%
con = seedPref.covMat;
cellList = seedPref.cellList;
[sortLinks idx] = sort(con(:),'ascend');
[iy ix] = ind2sub(size(con),idx);
num = length(cellList);

cutLinks = sortLinks(sortLinks>0);
for i = 1:length(cutLinks)
    
    
    con2  = con;
    con2(con2 < cutLinks(i))=0;
    subplot(2,1,1)
    image(con2*100),pause(.0001)
    
    cellGroup = zeros(1,num);
    for g = 1:num;
            

        nextY = find(cellGroup ==0,1);
        while 1
            Y = nextY;
            cellGroup(Y) = g;
            con2(:,Y) = 0;
            [oldY nextY] = find(con2(Y,:)>0);
            if isempty(nextY),break,end
        end
        if ~sum(cellGroup==0),break,end
    end
    
        subplot(2,1,2)

                image(cellGroup*1); pause(.001)

    
    
end

