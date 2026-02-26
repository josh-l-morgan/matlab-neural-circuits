%% load data


global glob tis

figure
SPN = [glob.datDir 'Analysis\Data\preproc\'];
load([SPN 'ptDat.mat']);
load([SPN 'ROI.mat']);
%load([SPN 'ROI2.mat']);
load([SPN 'SOI.mat']);
load([SPN 'GOI.mat']);
load([SPN 'NOI.mat']);
load([SPN 'MOI.mat']);
load([SPN 'COI.mat']);



%%Load all sms
smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
clear sms
roiCid = ptDat(:,3);
vgcCids = COI.vgcCids;%unique(roiCid);%MOI.cids;
%runCids = unique(roiCid);%MOI.cids;
useSM = vgcCids * 0;
for i = 1:length(vgcCids);
    cid = vgcCids(i);
    %%Get distances between nodes
    disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(vgcCids)));
    fileName = sprintf('sm_cid%d.mat',cid);
    %sm = load([smDir fileName],'skel2skel','nep','syn2skel','syn');
    if exist([smDir fileName],'file')
        load([smDir fileName]);
        sm.nep.swcS = nep2swc(sm.nep);
        sms(i).sm = sm;
        useSM(i) = 1;    
    end
end

runCids = vgcCids(useSM>0);
useVgc = runCids;

s = makeSynapseClassifyer(COI); %make structure describing types of synapses







