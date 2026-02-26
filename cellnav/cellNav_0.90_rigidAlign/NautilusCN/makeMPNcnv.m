function makeMPNcnv(shouldMake)

if ~exist('shouldMake','var')
    shouldMake = 0;
end

global glob

MPN = glob.NA.MPN;
slash = regexp(MPN,'\');
WPN = [MPN(1:slash(end-1)) 'Analysis\'];
glob.NA.WPN = WPN;
save('MPN.mat','MPN', 'WPN')

if shouldMake
    if ~exist(MPN,'dir'),mkdir(MPN);end
    if ~exist(WPN,'dir'),mkdir(WPN);end
end