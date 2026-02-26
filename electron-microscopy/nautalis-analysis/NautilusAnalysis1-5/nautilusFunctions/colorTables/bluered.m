function[colTable] = bluered(N);


if ~exist('N','var')
    N = 100;
end

    colTable = zeros(N,3);
    colTable(:,1) = (1:N)/N;
    colTable(:,3) = (N:-1:1)/N;