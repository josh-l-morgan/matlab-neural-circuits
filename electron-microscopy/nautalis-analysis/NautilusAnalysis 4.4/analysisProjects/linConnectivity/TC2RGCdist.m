
load('MPN.mat')
%nepFileName = sprintf('%snep\\skelNep%d.mat',MPN,skelList);

load([MPN 'nep\skelNep125.mat'])
load([MPN 'nep\useSynPos125.mat'])
pos = nep.nodePos;
edges = nep.edges;



points1 = useSynPos.synPos125toTCR;
points2 = useSynPos.fromRGC_dat;


points1 = points1(sum(points1,2)>1,:);
points2 = points2(sum(points2,2)>1,:);





%% measure topo distances

[edges2 pos2] = doubleEdge(edges,pos,2);%% Double nodes
pp = shortestPoint2Point(edges2,pos2,points1,points2);

maxMin = pp.minMax;


%%
colVal = ceil(maxMin);
colVal(colVal<1) = 1;
colVal(colVal>100) = 100;
cmap =  jet(100);


scatter3(pos2(:,1),pos2(:,2),pos2(:,3),2,'k','filled')
%axis('color','k')

hold on
%scatter3(points1(:,1),points1(:,2),points1(:,3),10,'r','filled')
%scatter3(points2(:,1),points2(:,2),points2(:,3),10,'g','filled')

scatter3(points1(:,1),points1(:,2),points1(:,3),50,cmap(colVal,:),'filled')
hold off
pause(1)

view([0 0])
axis('off')

%% plot

minRange = [0:5:100];
histMin = histc(maxMin,minRange);
ratMin = histMin/sum(histMin);
bar(minRange,ratMin,'k');

checkLength = 5;
sum(maxMin<=checkLength)/length(maxMin)





