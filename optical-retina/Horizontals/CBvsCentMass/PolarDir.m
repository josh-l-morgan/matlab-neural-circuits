

scatter(OffSet,Direct)

hist(Direct./OffSet)

scatter(OffSet,Direct./OffSet)


%{
Differences between first and last				in MICRONS
'ER'	'OffSet'	'DistsCB'	'DistsTer'	'Direct'
%}

O=ones(size(CS,1),1);
scatter([O ;O*2],[CS(:,2) ;CS(:,5)]) % OffSet to Direct

scatter([O ;O*2],[CS(:,3) ;CS(:,4)]) % OffSet to Direct
hold on
plot([CS(:,3) CS(:,4)]')
hold off

scatter(CS(:,3), CS(:,4))