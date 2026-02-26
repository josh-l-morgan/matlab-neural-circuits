
%save vast import data for python manipulation
vastSubsPy = []
vast_path = "\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\KxR_P11LGN\CellNav_KxR\Volumes\HighRes2023\Merge\"
sm_path = "\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\KxR_P11LGN\CellNav_KxR\Volumes\HighRes2023\Analysis\SMs\"
%filename = [MPN 'vastSubsPy.mat']
target_cells = [3204, 3206];
target_cells = [3004 3107 3210 3211 3213 3018 3097 3032 3098 3200 3202 3203 3204 3206 3209 3219];
load(vast_path+"vastSubs.mat","vastSubs")
load(vast_path+"obI.mat","obI")


vRes = obI.em.res;
vRes([1 2]) = vRes([1:2]) * 2^obI.em.mipLevel;

for i=1:length(vastSubs)
    subs = vastSubs{i};
    if ~isempty(subs)
        subs = subs(:,[2 1 3]) .* repmat(vRes/1000,[size(subs,1) 1]);
    end
    vastSubsPy{i}.vox = subs;
    vastSubsPy{i}.cellNum = obI.nameProps.cellNum(i);
end
syns = obI.nameProps.synProp;
cellNum = obI.nameProps.cellNum;
save(vast_path+"vastSubsPy.mat","vastSubsPy",'-v7.3')
save(vast_path+"syns.mat","syns",'-v7.3')

dsRes = obI.em.dsRes;
upsample_factor = dsRes;
res = [obI.em.res(1),obI.em.res(2),1];

for i=1:length(target_cells)
    sm_filename = sm_path+"sm_cid"+int2str(target_cells(i))+".mat";
    cell_data = load(sm_filename);
    taget_cells_data{i}.cellNum = target_cells(i);
    taget_cells_data{i}.pre = cell_data.sm.syn.pre;
    taget_cells_data{i}.post = cell_data.sm.syn.post;
    taget_cells_data{i}.type = cell_data.sm.syn.synType;
    taget_cells_data{i}.pos = cell_data.sm.syn.pos;
end
save(vast_path+"sm_syns.mat","taget_cells_data",'-v7.3')