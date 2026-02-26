function[Ic] = showCell(obI,dsObj,cellId, col,dim);

%%
if ~exist('dim','var')
    dim = 1;
end

if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end
%%
fsize = [1700 1700 1300];
I = zeros(fsize(dims),'uint8');
Ic = cat(3,I,I,I);dim 


uCellId = unique(cellId(cellId>0));

for i = 1:length(uCellId)
    targ = find(obI.cell.name==uCellId(i));
    obTarg = obI.cell.obIDs{targ};
    for o = 1:length(obTarg)
        if obTarg(o)<=length(dsObj)
        sub = double(dsObj(obTarg(o)).subs);
        if ~isempty(sub)
        inds = sub2ind(fsize(dims),sub(:,dims(1)),sub(:,dims(2)));
        
        I(inds) = 1000;
        
        Ic(:,:,1) = Ic(:,:,1) + I*col(i,1);
        Ic(:,:,2) = Ic(:,:,2) + I*col(i,2);
        Ic(:,:,3) = Ic(:,:,3) + I*col(i,3);
        
        image(Ic)
        pause(.1)
        I = I * 0;
        end
        end
    end
       
end