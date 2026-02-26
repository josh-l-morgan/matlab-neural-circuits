function[types] = typeLists(obI,ids);

if ~exist('obI','var')
load('MPN.mat')
load([MPN 'obI.mat'])
end



%% cell types
if isfield(obI.cell, 'mainObID')
    rgcs = obI.cell.name(obI.nameProps.rgc(obI.cell.mainObID));
    tcrs = obI.cell.name(obI.nameProps.tcr(obI.cell.mainObID));
    lins = obI.cell.name(obI.nameProps.lin(obI.cell.mainObID));
end

rgcs = setdiff(rgcs,0);
tcrs = setdiff(tcrs,0);
lins = setdiff(lins,0);

types = ids * 0;
[a b] = intersect(ids,rgcs);
types(b) = 1;

[a b] = intersect(ids,tcrs);
types(b) = 2;
[a b] = intersect(ids,lins);
types(b) = 3;


