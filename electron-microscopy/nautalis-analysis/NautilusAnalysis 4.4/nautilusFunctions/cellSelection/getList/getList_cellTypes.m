function[checkIDs checkProp]  = getList_cellTypes();

load('MPN.mat')
load([MPN 'obI.mat'])
% 
% checkIDs = unique(obI.nameProps.cellNum);
% 
% 
% tcrList = unique(obI.nameProps.cellNum(obI.nameProps.tcr));
% rgcList = unique(obI.nameProps.cellNum(obI.nameProps.rgc));
% linList = unique(obI.nameProps.cellNum(obI.nameProps.lin));
% ldmList = unique(obI.nameProps.cellNum(obI.nameProps.ldm));
% 
% checkIDs = setdiff(checkIDs,[0])';
% checkProp = checkIDs * 0;
% 
% for i = 1:length(checkIDs)
%     
%    if sum(find(rgcList == checkIDs(i)));
%        checkProp(i) = 1;
%    elseif sum(find(tcrList == checkIDs(i)));
%        checkProp(i) = 2;
%    elseif sum(find(linList == checkIDs(i)));
%        checkProp(i) = 3;
%    elseif sum(find(ldmList == checkIDs(i)));
%        checkProp(i) = 4;
%    end
% end


%%

checkIDs = unique(obI.cell.name);

checkIDs = setdiff(checkIDs,[0])';
checkProp = checkIDs * 0;

for i = 1:length(checkIDs);
    targ = obI.cell.mainObID(find(obI.cell.name == checkIDs(i)));
    if obI.nameProps.rgc(targ)
        checkProp(i) = 1;
    elseif obI.nameProps.tcr(targ)
        checkProp(i) = 2;
    elseif obI.nameProps.lin(targ)
        checkProp(i) = 3;
     elseif obI.nameProps.ldm(targ)
        checkProp(i) = 4;
    end
        
end


[checkIDs checkProp];



