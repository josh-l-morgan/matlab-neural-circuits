function[] = rleDat2Subs(TPN);


%TPN = GetMyDir;

load([TPN 'rleDat.mat'])

rleDat.parseTime = clock;

%% convert to voxel list

blockDir = [TPN 'blockTemp\'];
if ~exist(blockDir,'dir'),mkdir(blockDir),end

segNum = rleDat.segNum;
blockNum = length(rleDat.rleimg);


poolSize = matlabpool('size');
blockRange = 1:poolSize:blockNum;
blockRange(end) = blockNum+1;

clear fileName
for r = 1:length(blockRange)-1;
        fileName{r} = sprintf('%s%d.mat',blockDir,r);
end


for r = 1:length(blockRange)-1;
    
    disp(sprintf('%d of %d',r,length(blockRange)-1));
    tic
    
    clear stackO
    for i = blockRange(r):blockRange(r+1)-1
        param = rleDat.param(i);
        rleimg = rleDat.rleimg{i};
        obV = rle2subs(rleimg,param,segNum);
        stackO{i} = obV;
    end
    save(fileName{r},'stackO');
    toc
    
end

%% Collect all subs
clear vastSubs
vastSubs{segNum-1} = []; %%!!!!!screwed up by 0 indexing
for r = 1 : length(fileName)
    r
    load(fileName{r});
    for s = 1:length(stackO)
        if ~isempty(stackO{s})
            newSubs = stackO{s};
            
            for o = 1:length(newSubs);
                if ~isempty(newSubs{o})
                    vastSubs{o} = cat(1,vastSubs{o},newSubs{o});
                end
            end
            
        end
    end
end

save([TPN 'vastSubs.mat'],'vastSubs','-v7.3')


%%

% 
% for s = 1:length(vastSubs)
% targ = s;
% subs = vastSubs{targ}(:,1:2);
% if ~isempty(subs)
% subs(:,1) = subs(:,1)-min(subs(:,1))+1;
% subs(:,2) = subs(:,2) - min(subs(:,2))+1;
% maxI = max(subs,[],1);
% showI = zeros(maxI);
% showI(sub2ind(maxI,subs(:,1),subs(:,2))) = 1000;
% image(showI),pause(.01)
% end
% end






