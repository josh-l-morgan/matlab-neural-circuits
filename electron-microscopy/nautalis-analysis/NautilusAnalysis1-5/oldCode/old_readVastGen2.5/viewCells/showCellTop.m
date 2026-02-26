function[Ic,maxHeight] = showCellTop(obI,dsObj,cellId, col,dim,fsize);

 minInt = 20;
contrastFactor = 5


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
%fsize = [1700 1700 1300];
I = zeros(fsize(dims));
newHeight = I;
maxHeight = I;
Ic = cat(3,I,I,I);


uCellId = unique(cellId(cellId>0));

for i = 1:length(uCellId)
    targ = find(obI.cell.name==uCellId(i));
    obTarg = obI.cell.obIDs{targ};
    for o = 1:length(obTarg)
        if obTarg(o)<=length(dsObj)
            sub = double(dsObj(obTarg(o)).subs);
            if ~isempty(sub)
                inds = sub2ind(fsize(dims),sub(:,dims(1)),sub(:,dims(2)));
                uinds = unique(inds);
                
                if length(uinds)>1
                    hinds = hist(inds,uinds);
                else
                    hinds = length(inds);
                end
                
                I(uinds) = I(uinds) + hinds';
                newHeight(inds) = max(newHeight(inds),sub(:,dim));
                
            end
            
           
            
        end
        
       
            
        
        
    end
            
            I = I * contrastFactor;
            I((I>0)) = I(I>0) + minInt;
            I(I>255) = 255;
            I((I>0) & (I<minInt)) = minInt;
            
            
            useI = newHeight>maxHeight;
            maxHeight(useI) = newHeight(useI);
            
            for c = 1:3
                Itemp = Ic(:,:,c);
                Itemp(useI) = I(useI)*col(i,c);
                Ic(:,:,c) = Itemp;
            end
            
            image(uint8(Ic))
            pause(.001)
            I = I * 0;
    
end

%Ic = uint8(Ic);



