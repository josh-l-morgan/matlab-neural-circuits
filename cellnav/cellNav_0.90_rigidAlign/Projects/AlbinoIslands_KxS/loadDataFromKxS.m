

global glob tis
load('MPN.mat');
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])


% figure
% SPN = [glob.datDir 'Analysis\Data\preproc\'];
% load([SPN 'ptDat.mat']);
% load([SPN 'ROI.mat']);
% %load([SPN 'ROI2.mat']);
% load([SPN 'SOI.mat']);
% load([SPN 'GOI.mat']);
% load([SPN 'NOI.mat']);
% load([SPN 'MOI.mat']);
% load([SPN 'COI.mat']);

tcrs = tis.cids(find(tis.cells.type.typeID==2));
lins = tis.cids(find(tis.cells.type.typeID==3));

runCids = tcrs;

%%Load all sms
smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];
clear sms

for i = 1:length(runCids);
    cid = runCids(i);
    %%Get distances between nodes
    disp(sprintf('loading data for cell %d.  Cell %d of %d.',cid,i,length(runCids)));
    fileName = sprintf('sm_cid%d.mat',cid);
    %sm = load([smDir fileName],'skel2skel','nep','syn2skel','syn');
    if exist([smDir fileName],'file')
        load([smDir fileName]);
        if 0
        sm.nep.swcS = nep2swc(sm.nep);
        end
        sms(i).sm = sm;
        useSM(i) = 1;    
    end
end

%% Parse synaspse type
for i = 1:length(sms)
    sm = sms(i).sm;
    if ~isempty(sm)
    for o = 1:length(sm.syn.obID)
        nam = obI.nameProps.names{sm.syn.obID(o)};
        findRGC = regexp(lower(nam),'rgc');
        if ~isempty(findRGC)
           sm.syn.preClass(o) = 1;
        end
    end
    sms(i).sm = sm;
    end
end

%%Parse synaspse type

for o = 1:length(tis.syn.obID)
    nam = obI.nameProps.names{tis.syn.obID(o)};
    findRGC = regexp(lower(nam),'rgc');
    if ~isempty(findRGC)
        tis.syn.preClass(o) = 1;
        tis.syn.synType(o) = 3;
    else
        tis.syn.synType(o) = 4;
    end
end
if 1
    save([glob.fvDir 'tis.mat'],'tis')
end


%s = makeSynapseClassifyer(COI); %make structure describing types of synapses


