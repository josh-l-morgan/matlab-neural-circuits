%%Draw Mask

TPN = GetMyDir

%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'

if exist([TPN 'mask']) %% mask is image in folder 'mask' if exists
    d=dir([TPN 'mask']); %get number of files in directory
    d=d(3:size(d,1));
    
    clear I IM 
    for i=1:size(d,1)
        IM(:,:,i)=imread([TPN 'mask\' d(i).name]);
        PercentRead=i/size(d,1)*100
    end
else %default mask is D
    load([TPN 'data\D.mat']); %load thresholded dendrite arbor
    IM=D>0; %(make sure its in binary)
    clear D
end


imwriteNp(TPN,IM,'MaskUsed')