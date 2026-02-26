%%Make survival curve of presynaptic inputs to a target cell based on con
%%and cellList

targ = find(cellList ==  108)

preHits = con(:,targ);

binHit = 0:max(preHits)
clear survivalC
for i = 1:length(binHit)
    
survivalC(i) = sum(preHits >=binHit(i));
end
plot(binHit,survivalC,'r','LineWidth',5)
hold on
histSyn = histc(preHits,binHit)
bar(binHit,histSyn,'b')
set(gca,'XTick',binHit)
set(gca,'YTick',[0:max(survivalC)])
xlim([0 binHit(end)+1])
hold off


exp(1)^ (-1* dists/distConstant)