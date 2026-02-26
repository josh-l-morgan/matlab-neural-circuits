%% run all the mask folders in a directory 

clear all

DPN = GetMyDir;
dDPN = dir(DPN); dDPN = dDPN(3:end);
%% Run all folders
for d = 1:length(dDPN)
    display(['running ' num2str(d) ' of ' num2str(length(dDPN))])
    TPN = [DPN dDPN(d).name '\'];
    if isdir(TPN)
        if ~exist([TPN 'data\cellProps.mat'])
            cP = anaOneBip(TPN,0);
        end
    end
end




