
%{
vasttools
parpool
%}

%% set scale variables
mipLev = 3;
%%%% standard resolutions
%LGN res : [6 4 30].*[4 4 1] 

res = [24 16 30]; 
vRes = [2^mipLev 2^mipLev 1].* res; %res times mip level
dsRes = [200 200 200];


%% Get directories

SPN = GetMyDir;
MPN = SPN;
TPN = SPN;


%% Run import
startTime = clock;

useSaved = 1;

% 
% mpSize = parpool('size');
% if ~mpSize
%     'opening matlab pool'
%     parpool close force
%     parpool 
% end

%}


%% Read in vast Colors
%Please connect to VAST with vasttools first!
if ~exist([MPN 'rleDat.mat'],'file')
    getRLEtoSubs(TPN,mipLev);
end

if ~exist([MPN 'vastSubs.mat'],'file')  
    rleTemp2Subs(TPN);
end

if useSaved & exist([MPN 'obI.mat'],'file')
    load([MPN 'obI.mat'])
else
    obI = makeOBI(SPN,MPN)
end


%% Down sample vastSubs


dsDim = dsRes./vRes;

useSaved = 0;
if useSaved & exist([MPN 'dsObj.mat'],'file')
    load([MPN 'dsObj.mat'])
else
    dsObj = downSampObj(MPN, dsDim);
end


obI.em.res = res; %or [4.6 4 30]?
obI.em.vRes =vRes;
obI.em.dsRes =dsRes/1000;
save([MPN 'obI.mat'],'obI')





