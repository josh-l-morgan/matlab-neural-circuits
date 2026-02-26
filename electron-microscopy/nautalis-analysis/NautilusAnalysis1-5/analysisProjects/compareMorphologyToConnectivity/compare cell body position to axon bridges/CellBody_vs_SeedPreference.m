
load([SPN 'traceCB.mat'])

seedCells = [108 201]

centroids = traceCB.CBcentroids;
cbID = traceCB.cbID;

seedTarg = find(cbID == seedCells(1));
seedPos = centroids(seedTarg,:);

cbDists = sqrt((centroids(:,1)-seedPos(1)).^2 + ...
    (centroids(:,2)-seedPos(2)).^2 +...
(centroids(:,3)-seedPos(3)).^2 );
    


cbDists; %% distances from cell bodies to seed cells
cbID; %% cell bodies mapped to network


conPref = seedPreferences(seedCells, useList);

cellPref = conPref.cellPrefNoExclusion(1,:)./sum(conPref.cellPrefNoExclusion,1);
prefNan = isnan(cellPref)


%%
outOfNetDists = cbDists(cbID==0);
inNetDists = cbDists(cbID>0);
inNetIDs = cbID(cbID>0);



scatter(outOfNetDists,ones(size(outOfNetDists)),'r')
hold on
scatter(inNetDists,ones(size(inNetDists))+1,'b')
hold off

xlim([0 100])
ylim([0 3])




%%

subplot(1,1,1)

for i = 1:400
   sumin(i) = sum(inNetDists<=i);
   sumout(i) = sum(outOfNetDists<=i);
    
end

infrac = sumin./(sumin+sumout);

plot(sumin,'b')
hold on
plot(sumout,'r')
hold off


plot(infrac)

%% nearness list

dsDif = (429-374) + (456 - 449);
fullDif = dsDif * 32;

cbCenters = centroids;
[sortDists idx] = sort(cbDists,'ascend');
sortCenters = cbCenters(idx,:);
sortCentersFull(:,1) = sortCenters(:,1) / .016 - fullDif;
sortCentersFull(:,2) = sortCenters(:,2)  / .016;
sortCentersFull(:,3) = sortCenters(:,3)  / .03;
%sortCenters = sortCenters(sortDists<100,:);


sortCentersDS(:,1) = sortCentersFull(:,1) / 32 ;
sortCentersDS(:,2) = sortCentersFull(:,2)  / 32;
sortCentersDS(:,3) = sortCentersFull(:,3)  / 16;


sortCentersFull = round(sortCentersFull);
sortCentersDS = round(sortCentersDS);

sortIDs = cbID(idx)';


%% distance to 201

targOther = find(sortIDs == 201);

dist2other = sqrt((sortCenters(:,1) - sortCenters(targOther,1)).^2 + ...
    (sortCenters(:,2) - sortCenters(targOther,2)).^2 + ...
    (sortCenters(:,3) - sortCenters(targOther,3)).^2);


%% Navigate CB




dsRes = [832 896 610]
fullRes = [26624 28672 9802]

scaleFac = round(fullRes./dsRes)

pos108 = [12944, 18328, 3948]

changePos = [369, 597, 331];

dsPos = round(changePos./scaleFac)
fullPos = round(changePos.*scaleFac)








