function[Ij] = juggleCh(I,j)
%%Juggle channels takes a 3 by 2 matrix that tells each channel where to
%%move to.  


if nargin == 1

    if sum(sum(sum(I(:,:,3,:))))==0;
        j=[2 1; 1 2; 1 3];
    else
        j=[1 3;2 2; 3 1];

    end
end
if size(I,3) >2
Ij(:,:,j(1,2),:)=I(:,:,j(1,1),:);
Ij(:,:,j(2,2),:)=I(:,:,j(2,1),:);
Ij(:,:,j(3,2),:)=I(:,:,j(3,1),:);
else
    Ij = I;
end