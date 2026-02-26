[synPos] = getSynPos(selection)


if selection = 1
    
load('C:\Users\joshm\Documents\MATLAB\jlm_Code\EM\NautilusAnalysis1-5\data\synapsePositions\track125syn.mat')
    
end


groupTags = {'pre','post','non','unk','?'}
isGroup = zeros(size(dat,1),length(groupTags));
for i = 1:size(dat,1)
    
    nam = dat{i,2};
    if isempty(nam)
        isGroup(i,:) = 0;
    else
    for g = 1:length(groupTags)
        isGroup(i,g) = sum(regexp(nam,groupTags{g})>0);
    end
    end
        
    
    
end

[y x] = find(isGroup>1)
