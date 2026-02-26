function[I] = fastCon(I,K);

% get size
[iys ixs izs] = size(I);
[kys kxs kzs] = size(K);

% pad with zeros

K = double(padarray(K,[iys ixs],'post'));
I = double(padarray(I,[kys kxs],'post'));

% convolve by fft
cI = ifftn(fftn(I).*fftn(K));

% recover non zeros
I = cI(round(kys/2):round(kys/2)+iys - 1,...
    round(kxs/2):round(kxs/2)+ixs - 1,...
    round(kzs/2):round(kzs/2)+izs - 1);
I = real(I);
% % image
% for i = 1:size(cutI,3)
%    image(fitH(cutI,i)),pause; 
% end



