function[] = showSkelChunks3D(subs,vox2node)


offset = min(subs,1);
for i = 1:3
    subs(:,i) = subs(:,i) - offset(i) + 1;
end



%% Render 3D
hold off
clf

axis off
pause(.01)
hold on
nCol = hsv(1000);
nodes = unique(vox2node)
for n = 1:length(nodes)
    nSubs = subs(vox2node==nodes(n),:);
    nColi = nCol(ceil(rand*size(nCol,1)),:);
    renderCon(nSubs,[],nColi,.2)
    pause(.01)

end
view([0 0])

ax = gca;
ax.Clipping = 'off';
hold off


