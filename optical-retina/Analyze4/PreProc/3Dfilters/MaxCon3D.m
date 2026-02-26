function[If] = MaxCon3D(I,Con)

%%Replaces outlyers

%%Set default rank
if nargin==1
    Con=.5; %% defines the rank of the pixel that replaces center
end

rankHigh=9-Out;
rankLow=1+Out; % define lower threshold

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
            MaxSurround=ordfilt2(Raw,rankHigh,[1 1 1; 1 1 1; 1 1 1]);
            MinSurround=ordfilt2(Raw,rankLow,[1 1 1; 1 1 1; 1 1 1]);
            Raw(Raw>MaxSurround)=MaxSurround(Raw>MaxSurround);
            Raw(Raw<MinSurround)=MinSurround(Raw<MinSurround);
            If(:,:,ch,i)=Raw;
           
        end    
    end
end
