%% Seperate signals from two sources that are mixed into two channels. 

%{
VARIABLES
g = signial in green channel
G = real green signal
Gg = translatingin G into g 
Rg = translating R into g

OUTPUT
g = G * Gg + R * Rg,   r = G * Gr + R * Rr;

SOLUTION
G = (g - R * Rg)/Gg
r = (g - R * Rg)/Gg * Gr + R * Rr
r = g * Gr * 1/Gg - R * Rg * Gr * 1/Gg + R * Rr
R * Rg * Gr * 1/Gg + R * Rr = g * Gr * 1/Gg - r
R * (Rg * Gr * 1/Gg + Rr) = g * Gr * 1/Gg - r
R = (g * Gr * 1/Gg - r) / (Rg * Gr * 1/Gg + Rr)


FINAL
R = (g * Gr * 1/Gg - r) / (Rg * Gr * 1/Gg + Rr)
G = (g - R * Rg)/Gg

%}


%% Get image
clear all
colormap gray(256)
[DFN DPN] = uigetfile;
TargFile=[DPN DFN]

Iall = tiffread2(TargFile);
for i = 1: size(Iall,2)
I(:,:,i)=Iall(i).data;
end

r=I(:,:,1:size(I,3)/2);
g=I(:,:,size(I,3)/2+1:size(I,3));
g=double(g);
r=double(r);

% %% remove background
% V=sort(r(:));
% Mr=mean(V(1:fix(size(V,1)/10)+1));
% r=r-Mr;
% V=sort(g(:));
% Mg=mean(V(1:fix(size(V,1)/10)+1));
% g=g-Mg;
% 
% 
% 
% image(max(g,[],3)/10)
% image(max(r,[],3)/10)
% 
% %% Find ratio
% % rat=g./r;
% % 
% % IDX=kmeans(rat(:),2);
% % ratS=rat*0;
% % ratS(IDX==1)=1;
% % ratS(IDX==2)=2;
% % image(sum(ratS==2,3)*20)
% %% match medians and set Green minimallly brighter then blue with same
% 
% r2=r;
% %find backgrounds
% V=sort(r2(:));
% Mr=mean(V(1:fix(size(V,1)/10)+1));
% V=sort(g(:));
% Mg=mean(V(1:fix(size(V,1)/10)+1));
% 
% %find min difference
% Is=g-r2;
% Mindif=min(Is(:));
% [y x z] = find3(Is==Mindif);
% Br=r2(y(1), x(1), z(1));
% Bg=g(y(1), x(1), z(1));
% 
% %Fix contrast
% Dr=Br-Mr; %find differences
% Dg=Bg-Mg;
% 
% r2=r2*(Dg/Dr); % change blue contrast
% V=sort(r2(:));
% Mr2=mean(V(1:fix(size(V,1)/10)+1)); %find new backgrounds
% r2=r2+(Mg-Mr2); %Shift blue curve to match green
% r2(r2<1)=1;
% 
% %% Find ratio
% Rat=g./(r2.^3);
% 
% image(max(Rat,[],3)*200000)
% V=sort(Rat(:));
% Thresh=V(fix(length(V)-length(V)/100));
% 
% RatThresh=Rat>Thresh;
% for i = 1: size(RatThresh,3)
%    RatThresh(:,:,i)=filter2([1 1 1; 1 1 1; 1 1 1],RatThresh(:,:,i)); 
% end
% 
% GrTh=graythresh(g);
% image(sum(g>GrTh,3)*30)
% image(sum(RatThresh,3)*10)
% RatThresh=RatThresh*0;
% RatThresh(1:150,:,:)=1;
% RatThresh=RatThresh>0;
% R2s=r(RatThresh)./g(RatThresh);
% Gg=mean(g(RatThresh));
% Gr=mean(r(RatThresh));
% R2 = mean(R2s);
% 
% 
% %% use product to find double labled voxels
% P = r.*g;
% image(max(P,[],3)/20000)
% V=sort(P(:));
% Thresh=V(fix(length(V)-length(V)/1000));
% image(sum(P>Thresh,3)*20)
% R1s=r(P>Thresh)./g(P>Thresh);
% R1 = mean(R1s);
% Rg=mean(g(P>Thresh));
% Rr=mean(r(P>Thresh));
% 
% 

%% Unmix
% 
Gr = .1
Gg = .9
Rg = .8
Rr = .4

Runscaled = [Rg; Rr];
Rlength = sqrt(Runscaled'*Runscaled);
Ru = Runscaled / Rlength;

Gunscaled = [Gg; Gr];
Glength = sqrt(Gunscaled'*Gunscaled);
Gu = Gunscaled / Glength;

A = [Ru Gu];

[ys xs zs]=size(r);
nPixel = numel(r);

b = [r(:)'; g(:)'];


x = A\b;

xR=zeros(ys, xs, zs);
xR(:)=x(1,:);

xG=xR;
xG(:)=x(2,:);


S1 = xG;
S2 = xR;

image(max(S1,[],3)* 256/max(S1(:)))
image(max(S2,[],3) * 256/max(S2(:)))

C(:,:,1)=max(S1,[],3)*256/max(S1(:));
C(:,:,2)=max(S2,[],3)*256/max(S2(:));
C(:,:,3)=C(:,:,2);
image(uint8(C))



