function[Is]=PreSyn(I)
%% Get Presynaptic field

'median filtering image'
[ys xs zs]=size(I);
for i = 1: zs     
    I(:,:,i)=medfilt2(I(:,:,i));
end
'median filtering done'


%% Iterative Erosion
'1pixel, 3d opening'

SE=strel('ball',3,3);
Con=ones(3,3,3);
frep=1;
Med=median(I(:));
It=I*0;
for i = 1:zs
    Samp=I(:,:,max(1,i-2):min(zs,i+2));
   It(:,:,i)=I(:,:,i)>median(Samp(:)); 
end

It=uint8(I>Med);
If=imfilter(It,Con);
Is=(If>=27) & (It>0);
Is=imfilter(Is,Con);
Is=Is>0;

Imax=max(Is,[],3);
Isum=sum(Is,3);

image(Isum*.04)

%imwriteNp(TPN,Is,'Is')
%imwriteNp(TPN,I,'I')
