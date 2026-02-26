   
clear pts
y  = rand(100,1)*100;
x = rand(100,1)*100;
pt = [y x];

for z = 1:10

    pt(:,1) = pt(:,1)+randn*10;
    pt(:,2) = pt(:,2)+randn*10;
    theta = randn * 2;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    pt = pt * R;
    pts{z} = pt;
    
end

tforms(1).tform = fitgeotform2d([0 0;0 1;1 0], [0 0;0 1;1 0],'similarity');
for z = 2:10
    tforms(z).tform = fitgeotform2d(pts{z}, pts{z-1},'similarity');
end

newPts = pts;

for z = 1:length(tforms)
    tform = tforms(z).tform;
    for n = z:length(tforms)
        if ~isempty(newPts{n})
            [nY nX] = transformPointsForward(tform,newPts{n}(:,1),newPts{n}(:,2));
            newPts{n}(:,1) = nY;
            newPts{n}(:,2) = nX;
        end
    end
end

oldP = [];
newP = [];
for z = 1:length(tforms)
    if ~isempty(newPts{z})
        currentPts = newPts{z}(:,1:2);
        newP = cat(1,newP,[currentPts ones(size(currentPts,1),1) * z]);
        currentPts = pts{z}(:,1:2);
        oldP = cat(1,oldP,[currentPts ones(size(currentPts,1),1) * z]);
    end
end



clf
ax = cla;
ax.NextPlot = 'add';

scatter3(newP(:,1),newP(:,2),newP(:,3),'r','.');
scatter3(oldP(:,1),oldP(:,2),oldP(:,3),'g','o');
