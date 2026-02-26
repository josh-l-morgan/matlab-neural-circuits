function[] = ProcessStrats(useStrat)

%%Plot stratifications of bipolar and rgcs relative to depth over apposition

%%load in sorted strat data
load('C:\Documents and Settings\joshm\My Documents\MyWork\Matlab\FixedPairAnalysis\appoStrat\sortedStrat.mat')
% 
% useStrat = [];
% for i = 1:length(sortedStrat)
%     if ~isempty(sortedStrat{i})
%         useStrat(length(useStrat)+1)=i;
%     end
% end



brackMax = 200;
shiftBip = zeros(brackMax*2+1,3);
collectBip = zeros(brackMax* 2+1,3,length(useStrat));
for i = 1:length(useStrat)
   cellSums = sortedStrat{useStrat(i)};
   meanOverlapPos = fix(sum([1:size(cellSums,1)]'.* cellSums(:,2))/sum(cellSums(:,2)));
   startAt = brackMax - meanOverlapPos;
   shiftBip(startAt+1:(startAt + size(cellSums,1)),1:3)  = cellSums; 
   collectBip(:,:,i) = shiftBip;
end

meanBip = mean(collectBip,3);

plot(meanBip(:,1),'r'),hold on
plot(meanBip(:,2),'b')
plot(meanBip(:,3),'g'),hold off
% 
% for i = 1:size(collectBip,3)
%    plot(collectBip(:,1,i),'r'),hold on
%     plot(collectBip(:,2,i),'b')
%     plot(collectBip(:,3,i),'g'),hold off
%     pause(.2)
% end
% hold off












