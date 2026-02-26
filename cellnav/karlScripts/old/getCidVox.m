function cidVox=getCidVox(cid,n,DS,curtis)
cidVox={};
for i=1:length(cid)
    curCid=cid(i);
    %cidObjIDs=cellfun(@isequal, tis.obI.nameProps.cids, {curCid}, size(tis.obI.nameProps.cids));
    CurCidID=find(curtis.cells.cids==curCid);
    if ~isempty(CurCidID)
        curCidObIDs=curtis.obI.cell.obIDs{CurCidID};
        curCidVoxCell=struct2cell(DS(curCidObIDs));
        bigCoords=zeros(1,3);
        for k=1:length(curCidVoxCell)
            bigCoords= [bigCoords; curCidVoxCell{:,:,k}];
        end
        cidVox{i}=bigCoords(1:n:end,:);
    end
end

end