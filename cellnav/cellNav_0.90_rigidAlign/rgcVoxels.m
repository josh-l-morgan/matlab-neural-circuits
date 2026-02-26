
% find rgc area of syns for target tcr cells

target_cell_num = [3204];
target_cell_num = [3004 3107 3210 3211 3213 3018 3097 3032 3098 3200 3202 3203 3204 3206 3209 3219];

synType = 3;
sms_path = "\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\KxR_P11LGN\CellNav_KxR\Volumes\HighRes2023\Analysis\SMs"
vast_path = "\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\KxR_P11LGN\CellNav_KxR\Volumes\HighRes2023\Merge"
load(vast_path+"\vastSubs.mat");
load(vast_path+"\obI.mat");
load(vast_path+"\dsObj.mat")
cellNum=obI.nameProps.cellNum;
vRes = obI.em.res;
down_sample = 1;
vRes([1 2]) = vRes([1:2]) * 2^obI.em.mipLevel;
data = []
for i=1:length(target_cell_num)
    tcr_voxels = []
    if down_sample
        tcr_voxels=dsObj(cellNum==target_cell_num(i)).subs;
        tcr_voxels = tcr_voxels(:, [2 1 3]);
    else 
        tcr_vastSubs=vastSubs(cellNum==target_cell_num(i));
        for k=1:length(tcr_vastSubs)
            sub = tcr_vastSubs{k};
            tcr_voxels = [tcr_voxels; sub];
        end
        
    end
    tcr_voxels = double(tcr_voxels) .* repmat(vRes/1000,[size(tcr_voxels,1) 1]);
    cell_num = target_cell_num(i);
    load(sms_path+"\sm_cid"+int2str(cell_num)+".mat");
    preIDs=sm.syn.pre;
    synType=sm.syn.synType;
    rgc_cellNum = unique(preIDs(synType==synType));
    rgc_cellNum=rgc_cellNum(rgc_cellNum>0)
    for j=1:length(rgc_cellNum)
        rgc_voxels = []
        
        if down_sample
            rgc_voxels=dsObj(cellNum==rgc_cellNum(i)).subs;
            rgc_voxels = rgc_voxels(:, [2 1 3]);
        else
            rgc_vastSubs=vastSubs(cellNum==rgc_cellNum(j));
            for k=1:length(rgc_vastSubs)
                sub = rgc_vastSubs{k};
                rgc_voxels = [rgc_voxels; sub];
            end
            
        end
        rgc_voxels = double(rgc_voxels) .* repmat(vRes/1000,[size(rgc_voxels,1) 1]);
        data{end+1}.rgcCellNum=rgc_cellNum(j);
        data{end}.rgcVoxels=rgc_voxels;
        data{end}.tcrVoxels=tcr_voxels;
        data{end}.tcrCellNum=cell_num;
        data{end}.synLoc=sm.syn.pos(sm.syn.pre==rgc_cellNum(j),:);
    end
end

save(sms_path+"\rgc_voxels_ds.mat", "data", "-v7.3")