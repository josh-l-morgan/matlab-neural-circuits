function makeVolMPNcnv(shouldMake)

if ~exist('shouldMake','var')
    shouldMake = 0;
end

global glob

MPN = [glob.dir.Volumes glob.vol.activeName '\Merge\'];
WPN =  [glob.dir.Volumes glob.vol.activeName '\Analysis\'];
glob.NA.WPN = WPN;
save('MPN.mat','MPN', 'WPN')
