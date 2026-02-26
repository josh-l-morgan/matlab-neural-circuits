function[obI] = makeOBI(TPN,MPN)

%MPN =       'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'

%fileName =  [MPN 'S8_joshmHome_14+05+28_export.txt'];
if ~exist('fileName','var')
textDir = dir([TPN '*.txt']);
    TFN = textDir(1).name;
    
    %[TFN TPN] = uigetfile('.txt');
fileName = [TPN TFN];
end


obI.colStruc = readVastColors(fileName);
obI.nameProps = getNameProps(obI.colStruc.names);


cellIDs = obI.nameProps.cellNum;
obI.cell.name = unique(cellIDs);
for i = 1:length(obI.cell.name)
     obIDs = find((cellIDs == obI.cell.name(i)) & obI.nameProps.ofID);
    obI.cell.obIDs{i} = obIDs;
%     subs = [];
%     for o = 1:length(obIDs)
%         subs = cat(1,subs,
%         
%     end
end

save([MPN 'obI.mat'],'obI');