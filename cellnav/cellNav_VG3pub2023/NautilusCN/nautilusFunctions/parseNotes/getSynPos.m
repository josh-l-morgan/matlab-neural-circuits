function[synPos] = getSynPos(selection)


if selection == 1
    
    load('.\data\synapsePositions\track125syn.mat')
    %{
    save('.\data\synapsePositions\track125syn.mat','dat')
    %}
end


groupTags = {'?','pre','post','non','unk','recip','lin'};
isGroup = zeros(size(dat,1),length(groupTags));
pos = zeros(size(dat,1),3);

for i = 1:size(dat,1)
    
    nam = dat{i,2};
    if ~isempty(nam)
        
        for g = 1:length(groupTags)
            isGroup(i,g) = sum(regexp(nam,groupTags{g})>0);
        end
    end
    
    num = dat{i,1};
    if ~isempty(num)
        parsePos = sscanf(num,'(%d, %d, %d)')';
        if length(parsePos) == 3;
            pos(i,:) = parsePos;
        end
    end
    
end
isGroup = isGroup>0;


clear synPos
synPos.groupTags = groupTags;
synPos.prePos = pos(~isGroup(:,1) & isGroup(:,2),:);
synPos.postPos = pos(~isGroup(:,1) & isGroup(:,3),:);
synPos.postRGCPos = pos(~isGroup(:,1) & isGroup(:,3) & ...
    ~(isGroup(:,5) | isGroup(:,4) | isGroup(:,6) | isGroup(:,7) ),:);
synPos.postUnkPos = pos(~isGroup(:,1) & isGroup(:,3) & (isGroup(:,5) | isGroup(:,4)),:);
synPos.preUnkPos = pos(~isGroup(:,1) & isGroup(:,2) & (isGroup(:,5) | isGroup(:,4)),:);
synPos.recipPos = pos(~isGroup(:,1) & isGroup(:,6),:);
synPos.preLinPos = pos(~isGroup(:,1) & isGroup(:,2) & isGroup(:,7),:);
synPos.postLinPos = pos(~isGroup(:,1) & isGroup(:,3) & isGroup(:,7),:);
synPos.allPos = pos(sum(pos,2)>0,:);

%%

conf = dat(:,3);
isConf = zeros(length(conf),1);
hits = [];
parsed = zeros(length(conf),1);
for i = 1:length(conf)
    nam = conf{i};
    if ~isempty(nam)
        if strcmp(class(nam),'double')
            isConf(i) = 1;
            parsed(i) = nam;
            hits = [hits nam];
        else
            nums = sscanf(nam,'%d');
            if sum(nums)
                isConf(i) = 1;
                hits = [hits nums(1)];
                parsed(i) = nums(1);
            end
            %spac = regexp(nam,' ') 
        end  
    end
end


numOutputs = sum(~isGroup(:,1) & isGroup(:,2) & isConf);
numTCOutputs = sum(~isGroup(:,1) & isGroup(:,2) & isConf & ~isGroup(:,5));
sum(hits == 481);
tcrs = parsed(~isGroup(:,1) & isGroup(:,2) & isConf & ~isGroup(:,5));
synPos.tcrs = tcrs; %outputs  tcrs



%% test
%scatter(pos(:,3),pos(:,2),'.','k')
