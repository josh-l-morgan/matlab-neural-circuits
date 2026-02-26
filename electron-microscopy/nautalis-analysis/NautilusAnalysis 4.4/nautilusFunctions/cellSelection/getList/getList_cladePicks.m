function[members] = getList_cladePicks(cladeNum);

load('.\data\clade\cladePick_six2.mat')

if exist('cladeNum','var')
    members = cladePick.members{cladeNum};
else
    
    members = cladePick.members;

end



