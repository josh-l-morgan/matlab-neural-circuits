function[] = stackObs(SPN)


%SPN = 'D:\Joshm\S8\export_TestRead\';
TPN = [SPN(1:end-1) '_mat\'];
PPN = [TPN 'obPlanes\'];

OPN = [TPN 'obFiles\'];
if ~exist(OPN,'dir'), mkdir(OPN),end


load([TPN 'tileInfo.mat'])

startStack = tileInfo.range.minZ;
startStack = 1
stopStack = tileInfo.range.maxZ;


%% allocate subs lists within objects
useIds = find(tileInfo.trackIds>0);
idSize = tileInfo.trackIds(useIds);

%% preallocate subs
clear vastSubs
for i = 1:length(useIds)
    o = useIds(i);
    vastSubs{o} = zeros(idSize(i),3,'uint16');
end
    whos vastSubs


%% read in subs
trackZ = cell(max(useIds),1);

subFill = zeros(100000,1);
for z = startStack:stopStack
    z
    newName = sprintf('%05.0f.mat',z);
    if exist([PPN newName],'file')
    load([PPN newName])
    for i = 1:length(vastOb.uniqueIds)
       id = vastOb.uniqueIds(i);
       sub = vastOb.subs{i};
       subSize = size(sub,1);
       vastSubs{id}(subFill(id)+1:subFill(id)+subSize,:) = sub;
       subFill(id) = subFill(id) + subSize;    
       
%     subplot(3,1,1)
% hist(sub(:,1),[min(sub(:,1)):1:max(sub(:,1))])
% subplot(3,1,2)
% hist(sub(:,2),[min(sub(:,2)):1:max(sub(:,2))])
% subplot(3,1,3)
% hist(sub(:,3),[startStack:1:stopStack])
%     pause
       
    end
    end
    

    
end

save([TPN 'vastSubs.mat'],'vastSubs','-v7.3')


%% find zs
% trackZ = cell(max(useIds),1);
% for z = startStack:stopStack
%     newName = sprintf('%05.0f.mat',z);
%     load([PPN newName])
%     for i = 1:length(vastOb.uniqueIds)
%        id = vastOb.uniqueIds(i);
%         trackZ{id} = [trackZ{id} z];
%     end
% end
% 

%% plane to ob files -Too Slow
% 
% for i = 1:length(useIds)
%     i
%     o = useIds(i);
%     subs = zeros(idSize(o),3);    
%     subFilled = 0;
%     zs = trackZ{o};
%     
%     for z = 1:length(zs)
%          newName = sprintf('%05.0f.mat',zs(z));
%         load([PPN newName])
%         targ = find(vastOb.uniqueIds == o);
%         sub = vastOb.subs{targ};
%         subFill = subFilled + size(sub,1);
%         subs(subFilled+1:subFill,:) = sub;
%         subFilled = subFill;
%     end
%     obSubs = subs;
%     obName = sprintf('%05.0f.mat',o) 
%     save([OPN obName],'obSubs')
%     
% end

