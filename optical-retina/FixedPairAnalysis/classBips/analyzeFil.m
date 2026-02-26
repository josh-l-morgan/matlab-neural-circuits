%%Plot with Radii

[TFN TPN] = GetMyFile;
load([TPN TFN]);
%%
aXYZ = round(FilStats.aXYZ * 10);
aRad = double(FilStats.aRad )* 1;
aEdges = FilStats.aEdges+1;

aXYZ(:,1) = aXYZ(:,1)-min(aXYZ(:,1))+1;
aXYZ(:,2) = aXYZ(:,2)-min(aXYZ(:,2))+1;
aXYZ(:,3) = aXYZ(:,3)-min(aXYZ(:,3))+1;

showFil = zeros(fix(max(aXYZ(:,1)))+1,fix(max(aXYZ(:,2)))+1);
[ys xs] = size(showFil);
%% Draw arbor
for i = 1:size(aRad,2)

    [dY dX] = find(fspecial('disk',round(aRad(i) * 2)));
    dY = dY + mean(aXYZ(i,1),aXYZ(i,2))-round(aRad(i));
    dX = dX + mean(aXYZ(i,2),aXYZ(i,2))-round(aRad(i));
    dY(dY<1) = 1; dX(dX<1) = 1;
    dY(dY>ys) = ys; dX(dX>xs) = xs;
    rDisk = sub2ind([ys xs], dY,dX);

    showFil(rDisk) = showFil(rDisk)+1;

end

image(showFil*100)

%%
% for i = 1:size(aEdges,1)
%     plot([aXYZ(aEdges(i,1),1) aXYZ(aEdges(i,2),1)],[aXYZ(aEdges(i,1),2) aXYZ(aEdges(i,2),2)],...
%         'LineWidth',aRad(aEdges(i,1))*20);
%    %scatter(aXYZ(i,1),aXYZ(i,2),'.','r','MarkerSize',10)
%    hold on
% end
% hold off