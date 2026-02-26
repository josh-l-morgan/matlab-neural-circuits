function[difRes] = trackNodeCladeDif(results,checkIDs,checkProp)


cg = results.cellGroups;
nodeIDs = results.nodeIDs;

isProp = zeros(1,length(nodeIDs));
nodeProps = zeros(length(nodeIDs),size(checkProp,2));
for i = 1:length(checkIDs);
    targ = find(nodeIDs==checkIDs(i));
    if ~isempty(targ)
        nodeProps(targ,:) = checkProp(i,:);
        isProp(targ) = 1;
    end
end



%%

allDifs = cg*0;
allCounts = cg*0;
for y = 1:size(cg,1)
   groups = cg(y,:);
   gnum = max(groups);
   
   for g = 1:gnum
      isgroup = (groups==g) & (isProp);
      gprops = nodeProps(isgroup,:);
      meanG = mean(gprops,1);
      gdifs = gprops-meanG;
      %gsquares = gdifs.^2;
      allDifs(y,isgroup) = gdifs';
      allCounts(y,isgroup) = length(gdifs);
   end
       
end
%%


relDifs =  repmat(abs(allDifs(1,:)),[size(allDifs,1) 1])-abs(allDifs);
useDifs = (allCounts>1) & (repmat(isProp,[size(allDifs,1) 1]));
meanImprovement = sum(relDifs(useDifs))/sum(useDifs(:));

difRes.allDifs = allDifs;
difRes.allCounts = allCounts;
difRes.relDifs = relDifs;
difRes.useDifs = useDifs;
difRes.meanImprovement = meanImprovement;

% 
% 
% %%
% image(allDifs*100+100);
% image(allCounts*10);
% 
% hasProp = find(isProp);
% for i = 1:length(hasProp);
%     scatter(allCounts(:,hasProp(i)),allDifs(:,hasProp(i)));
%     plot(allDifs(:,hasProp(i)));
%     plot(relDifs(:,hasProp(i)));
%     ylim([-1 1])
%     pause
%     
% end

