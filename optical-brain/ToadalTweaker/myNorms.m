function[N] = myNorms(A,nType);

N = double(A);
if ~exist('nType','var')
    nType = 1;
end

if nType == 1

    N = N-min(N(:));
    N = N*256/max(N(:));
elseif nType ==2 %mode

    N = N-mode(N(:));
    N(N<0) = 0;
    N = N * 256/max(N(:));
end
