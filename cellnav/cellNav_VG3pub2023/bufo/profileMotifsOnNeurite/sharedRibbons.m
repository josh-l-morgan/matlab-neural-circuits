

global tis


checkRGCLabels = COI.checkRGCLabels; 
checkBPCLabels = COI.checkBPCLabels; 

c = 0;
for i = 1:length(tis.syn.synProp)
    sp = tis.syn.synProp{i};
    if length(sp.postID)>1
        c = c+1
        spr{c} = sp;
    end
end




shareMat = zeros(length(checkBPCLabels),length(checkRGCLabels));

for s = 1:length(spr)

    posts = spr{s}.postID;
    preId = spr{s}.preID;
    if ~isempty(preId) & ~isempty(posts)
        postTypes = cid2type(posts,tis);
        preTypes = cid2type(preId,tis);
        preSubType = preTypes{4};
        postSubTypes = postTypes{4};


        if sum(strcmp(postSubTypes,'vgc'))
            preTypeIdx = find(strcmp(checkBPCLabels,preSubType));
            for p = 1:length(postSubTypes)
                hitTag = find(strcmp(checkRGCLabels,postSubTypes{p}));
                if hitTag
                    shareMat(preTypeIdx,hitTag) =  shareMat(preTypeIdx,hitTag)  + 1;
                end
            end
        end
    end

end

clf
subplot(2,1,1)
cmap = jet(1000);
cmap(1,:) = 0;
colormap(cmap)
image(shareMat*1000/max(shareMat(:)))
maxTriadNumber = max(shareMat(:))

xticks([1:length(checkRGCLabels)])
xticklabels(checkRGCLabels)

yticks([1:length(checkBPCLabels)])
yticklabels(checkBPCLabels)


subplot(2,1,2)
overlapMat = typeVol' .* typeVolBPC;
image(overlapMat * 1000/max(overlapMat(:)))

sum(shareMat(:))
sum(shareMat,1)

49/74





%
%
% clear useR useB
% for i = 1:length(COI.rgcGroupLabel)
%    useR(i) = sum(strcmp(checkRGCLabels,COI.rgcGroupLabel{i}));
% end
% for i = 1:length(COI.bpcGroupLabel)
%    useB(i) = sum(strcmp(checkBPCLabels,COI.bpcGroupLabel{i}));
% end
%
% useR = find(useR);
% useRNames = {COI.rgcGroupLabel{useR}};
% useB = find(useB);
% useBNames = {COI.bpcGroupLabel{useB}};
% rbMat = zeros(length(useB),length(useR));
% clear memRCids memBCids
% for r = 1:length(useR);
%
%
%     rCids = COI.rgcGroupCids{useR(r)};
%     memRCids{r} = rCids;
%     for b = 1:length(useB)
%         bCids = COI.bpcGroupCids{useB(b)};
%         if r == 1
%             memBCids{b} = bCids;
%         end
%
%         for y = 1:length(rCids)
%             for x = 1:length(bCids)
%
%                 hits = sum((syn.pre==bCids(x)) & (syn.post==rCids(y)));
%                 rbMat(b,r) = rbMat(b,r) + hits;
%             end
%         end
%
%
%     end
% end
%
% cmap = jet(1000);
% cmap(1,:) = 0;
% colormap(cmap)
% image(rbMat*1000/max(rbMat(:)))
% maxTriadNumber = max(rbMat(:))
%
%


















