%% Seperate signals from two sources that are mixed into two channels. 

%{
VARIABLES
g = signial in green channel
g1 = green signal from source 1
g2 = green signal from source 2

r = signal in red channel

R1 = ratio of g/r for source 1
R2 = ratio of g/r for source 2

OUTPUT
g = g1 + g2,   r = g1 * R1 + g2 * R2

SOLUTION
g1 = g - g2,
r = (g-g2) * R1 + g2 * R2
r = g * R1 - g2 * R1 + g2 * R2
g2 * R1 + g2 * R2 = g * R1 - r
g2 (R1 + R2) = g * R1 - r
g2 = ( g * R1 - r )/(R1 + R2)

FINAL
g2 = ( g * R1 - r )/(R1 + R2)
g1 = g - g2,

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
% image(sum(RatThresh,3)*10)
% R2s=r(RatThresh)./g(RatThresh);
% R2 = mean(R2s);
% 
% 
% %use product to find double labled voxels
% P = r.*g;
% image(max(P,[],3)/20000)
% V=sort(P(:));
% Thresh=V(fix(length(V)-length(V)/1000));
% image(sum(P>Thresh,3)*20)
% R1s=r(P>Thresh)./g(P>Thresh);
% R1 = mean(R1s);
% 
% 
% %% Manually set Rs
% % 
% % R1=1;
% % R2=10;
% 
% 
% %%Solve image
% 
% 
% g2 = ( g * R1 - r )/(R1 + R2);
% g1 = g - g2;
% 
% r1= g1 * R1;
% r2= g2 * R2;
% S1 = r1;
% S2 = g2 + r2;
% 
% image(max(S1,[],3)/7)
% image(max(S2,[],3)/10)
% 
% C(:,:,1)=max(r,[],3)*256/max(r(:));
% C(:,:,2)=max(S2,[],3)*256/max(S2(:));
% C(:,:,3)=C(:,:,2);
% image(uint8(C))

%% Unmix

RedSig=[80 116];
GreenSig=[33 87];
BackSig=[18 7];

RedSig2=RedSig-BackSig;
GreenSig2=GreenSig-BackSig;

RedS=RedSig2./(RedSig2 + GreenSig2);
GreenS=GreenSig2./(RedSig2+GreenSig2);


% Runscaled = [Rg; Rr];
% Rlength = sqrt(Runscaled'*Runscaled);
% Ru = Runscaled / Rlength;
% 
% Gunscaled = [Gg; Gr];
% Glength = sqrt(Gunscaled'*Gunscaled);
% Gu = Gunscaled / Glength;
% 
% A = [Ru Gu];

A = [RedS; GreenS];


[ys xs zs]=size(r);
nPixel = numel(r);

b = [g(:)'-BackSig(1); r(:)'-BackSig(2)];


x = A\b;

xR=zeros(ys, xs, zs);
xR(:)=x(1,:);

xG=xR;
xG(:)=x(2,:);


S1 = xG;
S2 = xR;

image(max(S1,[],3)/50)
image(max(S2,[],3)/50)

C(:,:,1)=max(S1,[],3)*300/max(S1(:))-20;
C(:,:,2)=max(S2,[],3)*256/max(S2(:));
C(:,:,3)=C(:,:,2);
image(uint8(C))



C(:,:,1)=S1(:,:,20)*256/max(S1(:));
C(:,:,2)=S2(:,:,20)*656/max(S2(:));
C(:,:,3)=C(:,:,2);
image(uint8(C))


