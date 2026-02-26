

SPN = 'D:\LGNs1\Segmentation\VAST\S8\otherTracers\exportCell148_14+06+25\';
if ~exist(SPN,'dir')
    SPN = GetMyDir;
end


%%(Start Matlab Pool)%%
tic
tif2point(SPN);
toc

tic
stackObs(SPN)
toc

