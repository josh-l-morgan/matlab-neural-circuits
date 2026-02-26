function[If] = Shave(I,Out)

%%Replaces outlyers

%%Set default rank
if nargin==1
    Out=1; %% defines the rank of the pixel that replaces center
end

rankHigh=9-Out;
rankLow=1+Out; % define lower threshold


Dims=size(size(I),2);%%Count Channels
if Dims>3 %if multiple channels
    Imax=max(I,[],Dims);
    ChUse=sum(sum(Imax,1),2)>0;
    [ys xs Chan zs]=size(I);

    If=I;
    for ch = 1:Chan
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


else  %If one channel

    [ys xs zs] = size(I);

    for i = 1: size(I,4)
        Raw=I(:,:,i);
        MaxSurround=ordfilt2(Raw,rankHigh,[1 1 1; 1 1 1; 1 1 1]);
        MinSurround=ordfilt2(Raw,rankLow,[1 1 1; 1 1 1; 1 1 1]);
        Raw(Raw>MaxSurround)=MaxSurround(Raw>MaxSurround);
        Raw(Raw<MinSurround)=MinSurround(Raw<MinSurround);
        If(:,:,i)=Raw;

    end


end


%% Median filter Laplas

