
clear all
load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])



%% make golgi library
cellList = obI.cell.name;
%cellList = [10	129	162	170]
figure
colMap = hsv(256);
Dim = 3;
fsize = max(cat(1,dsObj.subs),[],1);


golgiDir = ['Z:\joshm\LGNs1\Analysis\golgiLibraryDim3\'];
%golgiDir = [MPN 'golgi1\']
if ~exist(golgiDir,'dir'),mkdir(golgiDir),end

for i = 1:length(cellList)
        showCellNames = cellList(i);
        targ = obI.cell.mainObID(i);
        type = sprintf('rtlu_%d%d%d%d',obI.nameProps.rgc(targ),obI.nameProps.tcr(targ),obI.nameProps.lin(targ),obI.nameProps.unk(targ));
        fileName = sprintf('%s_id_%05.0f_dim%d.png',type,showCellNames,Dim);
 if ~exist([golgiDir fileName],'file')
    col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
    col = col(randperm(size(col,1)),:);
    col = [1 1 1];
    I_partCell = showCellSum(obI,dsObj,showCellNames,col,Dim,fsize);
    I_partCell = 256-I_partCell*3;
    image(uint8(I_partCell))
    disp(showCellNames)
    %fileName = sprintf('golgi_%05.0f_dim%d.png',showCellNames,Dim);
    imwrite(uint8(I_partCell),[golgiDir fileName])
    pause(.01)
 end
end