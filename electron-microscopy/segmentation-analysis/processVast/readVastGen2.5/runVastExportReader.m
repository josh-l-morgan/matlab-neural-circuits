

SPN = 'D:\LGNs1\Segmentation\VAST\S8\otherTracers\exportCell148_14+06+25\';

%%(Start Matlab Pool)%%
tic
tif2point(SPN);
toc

tic
stackObs(SPN)
toc

