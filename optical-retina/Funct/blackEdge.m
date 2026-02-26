function[I] = blackEdge(I)

[ys xs zs] = size(I);

I(:,:,1) = 0;
I(:,:,zs) = 0;
I(:,1,:) = 0;
I(:,xs,:) = 0;
I(1,:,:) = 0;
I(ys,:,:)=0;