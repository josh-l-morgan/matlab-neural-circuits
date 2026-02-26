%% load data


global glob tis

figure
SPN = [glob.datDir 'Analysis\Data\'];
% load([SPN 'SOI.mat']);
% load([SPN 'GOI.mat']);
% load([SPN 'NOI.mat']);
% load([SPN 'MOI.mat']);
load([SPN 'COI.mat']);

runCids = COI.targetTCs;
%runCids = 3204;

%%Load all sms
smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
clear sms
useSM = runCids * 0;
for i = 1:length(runCids);
    cid = runCids(i);
    %%Get distances between nodes
    disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
    fileName = sprintf('sm_cid%d.mat',cid);
    %sm = load([smDir fileName],'skel2skel','nep','syn2skel','syn');
    if exist([smDir fileName],'file')
        load([smDir fileName]);
        %sm.nep.swcS = nep2swc(sm.nep);
        sms(i).sm = sm;
        useSM(i) = 1;    
    end
end


%s = makeSynapseClassifyer(COI); %make structure describing types of synapses







