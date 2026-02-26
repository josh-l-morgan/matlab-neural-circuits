
%%

fsize = [2000 2000 2000]
ids = [1:650];
cols = [1 0 0; 0 1 0; 0 0 1];
cols = rand(length(ids),3);

dim = [3 2];

Itemp = zeros(fsize(dim(1)), fsize(dim(2)));
Ic = uint8(cat(3,Itemp,Itemp,Itemp)); 

for i = 1:length(ids)

    
indI = sub2ind([fsize(dim(1)) fsize(dim(2))],...
    dsObj(i).y,dsObj(i).x);
Itemp(indI) = 256;

Ic(:,:,1) = Itemp* cols(i,1);
Ic(:,:,2) = Itemp* cols(i,2);
Ic(:,:,3) = Itemp* cols(i,3);


image(Ic),pause(.1)

Ic = Ic * 0;
Itemp = Itemp*0;

end


