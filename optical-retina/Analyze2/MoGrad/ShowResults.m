%% Show results

load('.\Results.mat')
for r = 1 : size(RadRs,2)
   Radius = RadRs(r)
   scatter(Rdend(r,:),Rdot(r,:)),
   ylim([1 4])
   xlim([1 4])
   xlabel('Ratio of inner dendrites over outer')
   ylabel('Ratio of Inner dots over outer')
   pause%(.1) 
end