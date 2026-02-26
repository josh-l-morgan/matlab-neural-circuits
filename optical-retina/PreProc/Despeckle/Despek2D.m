function[If] = Dspk(I,rank)

%%Replaces all pixels above a certain neighbor hood rank with the value at
%%that rank


rank=7; %% defines the rank of the pixel that replaces center

colormap gray(256)
Imax=max(I,[],4);
ChUse=sum(sum(Imax,1),2)>0;

% subplot(2,1,1)
image((Imax))

[ys xs Chan zs]=size(I);


%% Median filter Laplas
'median filtering',pause(.1)

If=I;
for ch = 1:3
    if ChUse(ch)
        for i = 1: size(I,4)
            Raw=I(:,:,ch,i);
            MaxSurround=ordfilt2(Raw,rank,[1 1 1; 1 0 1; 1 1 1]);
            Raw(Raw>MaxSurround)=MaxSurround(Raw>MaxSurround);
            If(:,:,ch,i)=Raw;
           
        end    
    end
end

IfMax=max(If,[],4);

subplot(2,1,1)
image(Imax(:,:,1)*5)
subplot(2,1,2)
image(IfMax(:,:,1)*5)