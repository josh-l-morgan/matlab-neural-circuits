function[Ic] = showCellSum(obI,dsObj,cellId, col,dim,fsize);

 minInt = 50;
contrastFactor = 1;


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
                
            end
            
           
            
        end
        
       
            
        
        
    end
            I = I * contrastFactor;
            I((I>0)) = I(I>0) + minInt;
            I(I>255) = 255;
            I((I>0) & (I<minInt)) = minInt;
            
            
            Ic(:,:,1) = Ic(:,:,1) + I*col(i,1);
            Ic(:,:,2) = Ic(:,:,2) + I*col(i,2);
            Ic(:,:,3) = Ic(:,:,3) + I*col(i,3);
            
            image(uint8(Ic))
            pause(.001)
            I = I * 0;
    
end

%Ic = uint8(Ic);



