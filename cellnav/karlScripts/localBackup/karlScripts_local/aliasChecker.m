function outList = aliasChecker(inputCidList,aliasList)
if isempty(aliasList)
    load('Y:\MATLAB\cellNav\karlScripts\localBackup\karlScripts_local\aliasList.mat');
end
outList=inputCidList;
iters=numel(inputCidList);
for i=1:iters
    curCid=inputCidList(i);
    if curCid~=0
        aliasRow=find(any(aliasList==curCid,2));
        if isempty(aliasRow)
            outCid=curCid;
        else
            outCid=aliasList(aliasRow(1),1);
        end
    else
        outCid=0;
    end
    outList(i)=outCid;
end



end