function makeMPNcnv

% MPN = GetMyDir;
% WPN = GetMyDir;

% CPN = currentDirectory;


%
global glob

MPN = glob.NA.MPN;
slash = regexp(MPN,'\');
WPN = [MPN(1:slash(end-1)) 'Analysis\'];
glob.NA.WPN = WPN;
if ~exist(MPN,'dir'),mkdir(MPN);end
if ~exist(WPN,'dir'),mkdir(WPN);end
save('MPN.mat','MPN', 'WPN')
%end