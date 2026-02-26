


%% load
load('MPN.mat')
load('WPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


if exists([MPN 'dat.m'],'file')
    load([MPN 'dat.m'],'file')
else
    dat = [];
end

%% 

