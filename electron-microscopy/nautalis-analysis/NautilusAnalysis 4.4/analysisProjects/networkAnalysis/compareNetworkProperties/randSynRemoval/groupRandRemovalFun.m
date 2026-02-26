function[sumGroup] = groupRandRemovalFun(con,forTrans);
    



nn = size(con,1);
sn = sum(con(:));
 cellGroups = zeros(sum(con(:))+1,nn);
    symCon = con + con';
    [cellGroup] = segmentCon(symCon);
    cellGroups(1,:)  = cellGroup;
    sumGroup = symCon * 0;
    for s = 1:sn;
        [y x] = find(con>0);
        pick = ceil(length(y)*rand);
        con(y(pick),x(pick)) = con(y(pick),x(pick)) - 1;
        
        symCon = con + con';
        [cellGroup] = segmentCon(symCon);
        cellGroups(s+1,:)  = cellGroup;
        
        compGroup = cellGroup';
        shiftGroup = compGroup;
        sameGroup = zeros(length(cellGroup));
       
        for i = 1:length(compGroup)
            isSame = (shiftGroup-compGroup) == 0;
            sameGroup(i,:) = isSame;
            shiftGroup = circshift(shiftGroup,1);
        end
        newGroup = sameGroup(forTrans);
        sumGroup = sumGroup+newGroup;
        
    end