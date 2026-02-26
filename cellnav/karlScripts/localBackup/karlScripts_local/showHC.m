function map=showHC()
HCcolMap=load('HCcolmap.mat');
HCcolMap=HCcolMap.HCcolmap;
pCols=round(HCcolMap./255,3);
pColsSrt=pCols([2:2:end 1:2:end],:);
figure();
scatter(1:15,repmat(1,[15 1]),repmat(100,[15 1]),pColsSrt,'filled');
hold on
scatter(1:15,repmat(2,[15 1]),repmat(100,[15 1]),pCols,'filled');
ylim([0 3])
map=pCols;
end

