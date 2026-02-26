
MPN = GetMyDir;
TPN = 'D:\LGNs1\Analysis\exports\';
fileName = 'cellVectors.mat';


load([MPN 'dsObj.mat']);
load([MPN 'obI.mat']);
[cellList] = obI2cellList_tracedCells(obI);

uCellId = [cellList.preList cellList.postList];

clear vecCell
for i = 1:length(uCellId)
    
    subs = [];
    term = uCellId(i);
    
    targ = find(obI.cell.name==term);
    if isempty(targ)
        obTarg = [];
        disp(sprintf('Could not find %d',term))
    else
        obTarg = obI.cell.obIDs{targ};
    end
    
    for o = 1:length(obTarg)
        if obTarg(o)<=length(dsObj)
            sub = double(dsObj(obTarg(o)).subs);
            if ~isempty(sub)
                subs = cat(1,subs,sub);
            end
        end
        
    end
    
    
    fsize=(max(subs,[],1));
     inds = sub2ind(fsize,subs(:,1),subs(:,2), subs(:,3));
                uinds = unique(inds);
                [y x z]  = ind2sub(fsize,(uinds));
                usubs = [y x z];
    
              
                testInds = sub2ind(fsize(1:2),subs(:,1),subs(:,2));
                field = zeros(fsize(1:2));
                field(testInds) = 1000;
                image(field),pause(.01)
                
                
        vecCell(i).subs = subs;
        vecCell(i).ID = term;
                
end

save([TPN fileName],'vecCell')
                
       