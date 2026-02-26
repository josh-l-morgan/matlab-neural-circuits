function[If] = Gaus3D(I,siz,sig)

if nargin ==1
    siz=15;
end

if nargin <3
    sig=siz/4;
end

%%Count Channels
Dims=size(size(I),2);
if Dims>3
    Imax=max(I,[],Dims);
    ChUse=sum(sum(Imax,1),2)>0;
    [ys xs Chan zs]=size(I);
else
    ChUse=[1 0 0]>0;
    [ys xs zs] = size(I);
end


%% Median filter Laplas

If=I;
for ch = 1:3
    if ChUse(ch)
        for i = 1: size(I,4)
            Raw=I(:,:,ch,i);
            Raw=imfilter(Raw,fspecial('gaussian',siz,sig));
            If(:,:,ch,i)=Raw;
           
        end    
    end
end
