function cidVox=getCidVox(cid,density,dsObj,tis)
cidVox={};
for i=1:length(cid)
    curCid=cid(i);
    %cidObjIDs=cellfun(@isequal, tis.obI.nameProps.cids, {curCid}, size(tis.obI.nameProps.cids));
    CurCidID=find(tis.cells.cids==curCid);
    if ~isempty(CurCidID)
        curCidObIDs=tis.obI.cell.obIDs{CurCidID};
        curCidVoxCell=struct2cell(dsObj(curCidObIDs));
        bigCoords=zeros(1,3);
        for k=1:length(curCidVoxCell)
            bigCoords= [bigCoords; curCidVoxCell{:,:,k}];
        end
        cidVox{i}=bigCoords(1:density:end,:);
    end
end

end