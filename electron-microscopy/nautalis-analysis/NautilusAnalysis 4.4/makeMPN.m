%%Make MPN

% MPN = GetMyDir;
% WPN = GetMyDir;

% CPN = currentDirectory;

% 
%MPN =     'J:\Chas VAST Export Data P11 2\';
%MPN = 'D:\LGNs1\Export\export_joshm_LIN125_dendType_2018+11+08\';
%MPN = 'E:\LGNs1_analysis\mergeSeg_mat\';
%MPN = 'G:\IxQ\Matlab\Merge\'
MPN = 'C:\Users\morga\Downloads\cellNavData\Matlab\Volumes\HighResSCM\Merge\'
WPN = 'C:\Users\morga\Downloads\cellNavData\Matlab\Analysis\'

% 
% MPN = 'G:\IxQ\Matlab_lowRes\Merge\'
% WPN = 'G:\IxQ\Matlab_lowRes\Analysis\'

% WPN =     '..\..\..\joshm\LGNs1\Analysis\';
% MPN = GetMyDir
%WPN = 'J:\Chas VAST Export Data P11 2\';
%WPN = [ MPN 'Analysis\'];
%WPN = 'E:\LGNs1_analysis\Analysis\'
%MPN = GetMyDir;

save('C:\Users\morga\Downloads\cellNavData\Matlab\MPN.mat','MPN', 'WPN')