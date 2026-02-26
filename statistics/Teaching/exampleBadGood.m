

n1 = 20;
n2 = 10;
Tstd1 = 60;
Tstd2 = 40;
Tm1 = 100;
Tm2 = 150;


g1 = randn(n1,1) * Tstd1 + Tm1;
g2 = randn(n2,1) * Tstd2 + Tm2;

%%

g = [ 1 2];

se1 = std(g1)/sqrt(n1);
se2 = std(g2)/sqrt(n2);
m1 = mean(g1);
m2 = mean(g2);


%%
clf
ew = .1;
m = [m1 m2];
errlow = [m1-se1 m2 - se2];
errhigh = [m1+se1 m2+se2];

subplot(1,2,1)
bar(g,m)
hold on
plot([g-ew;g+ew],[errlow;errlow],'linewidth',2,'color','k')
plot([g-ew;g+ew],[errhigh;errhigh],'linewidth',2,'color','k')
plot([g;g],[errlow;errhigh],'linewidth',2,'color','k')

xlim([0 3])
ylim([0 max(g2)])

[tt p] = ttest2(g1,g2);
bString = sprintf('B was %0.1f larger than A (p = %0.4f)',m2-m1,p);
title(bString)


subplot(1,2,2)
gs = [g1*0+1; g2*0+2];
errlow2 = [m1-se1 * 1.96 m2 - se2 * 1.96 ];
errhigh2 = [m1+se1* 1.96  m2+se2* 1.96 ];
mw = .2;
swarmchart(gs, [g1;g2],'markeredgecolor','none','markerfacecolor',[.5 .6 1])
hold on

plot([g-ew;g+ew],[errlow2;errlow2],'linewidth',2,'color','k')
plot([g-ew;g+ew],[errhigh2;errhigh2],'linewidth',2,'color','k')
plot([g;g],[errlow2;errhigh2],'linewidth',2,'color','k')

%plot([g-mw; g+mw],[m; m],'k','linewidth',2)


seC = 1/2 * sqrt((var(g1)/n1 + var(g2)/n2));

md = m2-m1;
CI = [md - 1.96 * seC  md + 1.96 * seC]


xlim([0 3])
ylim([0 max(g2)])

gString = sprintf('g2 was %0.1f to %0.1f greater than g1',CI(1),CI(2));
title(gString)









