
clear all
%SPN = 'D:\LGNs1\Segmentation\VAST\S8\otherTracers\exportCell107_14+06+27\';
%if ~exist('SPN','dir')
    SPN = GetMyDir;
%end


%%(Start Matlab Pool)%%
tic
tif2point(SPN);
toc

tic
stackObs(SPN);
toc

MPN = [SPN(1:end-1) '_mat\'];

tic
Ifull = easyShow(MPN,1);
subplot(1,1,1)
image(Ifull)
toc


