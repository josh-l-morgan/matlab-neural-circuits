function[If] = MyMedian(I,Out)

%%Replaces outlyers

if nargin==1
    Out=3; %% 
end


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

                If(:,:,ch,i)=medfilt2(Raw,[Out Out]);

            end
        end
    end


else  %If one channel

    [ys xs zs] = size(I);

    for i = 1: size(I,4)
        Raw=I(:,:,i);
        If(:,:,i)=medfilt2(Raw,[Out Out]);

    end


end


%% Median filter Laplas

