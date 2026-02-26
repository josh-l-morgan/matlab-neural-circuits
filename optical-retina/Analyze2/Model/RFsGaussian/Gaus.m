function[G]=Gaus(SD,R,o)

%% Gaussian Equasion
N=1; % define the number of dimensions
r=1:R;
G(r) = (1/((2 * pi * SD^2)^(N/2)))* exp(-(r-o).^2/(2*SD^2));
G=G/max(G(:));