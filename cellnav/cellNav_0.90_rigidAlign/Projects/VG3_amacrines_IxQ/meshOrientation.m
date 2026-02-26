global glob tis

useCid = 4;
fileName = sprintf('%s%d.mat',glob.useFvDir,useCid);
fv = loadFV(fileName);


p = patch(fv);
view(30,-15);
axis vis3d;
colormap copper
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud


pos = fv.vertices;
fv.faces

edges = cat(1,[fv.faces(:,1) fv.faces(:,2)], [fv.faces(:,2) fv.faces(:,3)],...
    [fv.faces(:,1) fv.faces(:,3)]);
maxE = max(edges,[],1);
indE = sub2ind(maxE, edges(:,1), edges(:,2));
[y x] = ind2sub(maxE, indE);
edges = [y x];

dif1 = fv.vertices(edges(:,1),1) - fv.vertices(edges(:,2),1);
dif2 = fv.vertices(edges(:,1),2) - fv.vertices(edges(:,2),2);
dif3 = fv.vertices(edges(:,1),3) - fv.vertices(edges(:,2),3);

ePos  = (fv.vertices(edges(:,1),:) + fv.vertices(edges(:,2),:))/2;

hist(ePos(:,1))

%% range
rng = [floor(min(pos(:,1))) :2: ceil(max(pos(:,1)))];
dif = zeros(length(rng)-1,3);
for z = 1 : length(rng)-1
    isBin = ((pos(:,1)>rng(z) ) & (pos(:,1)<=rng(z+1)));
    dif(z,1) = sum(abs(dif1(isBin)));
    dif(z,2) = sum(abs(dif2(isBin)));
    dif(z,3) = sum(abs(dif3(isBin)));
    varDif1(z) = var(pos(isBin,1));
    varDif2(z) = var(pos(isBin,2));
    varDif3(z) = var(pos(isBin,3));

end

ang = atan(dif(:,1)./mean(dif(:,(2:3)),2));
rat = dif(:,1)./(dif(:,1) + mean(dif(:,2:3),2));

plot(varDif2)















