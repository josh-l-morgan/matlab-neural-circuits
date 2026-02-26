%%Find translations for improving alignment 
%%between each plane of segmentations stored as dsObj

fig = figure;

global tis glob
if 1
load('MPN.mat')
load([MPN 'Merge\dsObj.mat'])
end

%%Make list of which objects are in each plane = zHas

numZ = 480;
for i = 1:numZ
    zCents{i} = [];
end

numObs = length(dsObj);
mins = [inf inf inf];
maxes = [0 0 0];
clear zHas
for i = 1:numObs
    subs = double(dsObj(i).subs);
    if ~isempty(subs)
        mins = min(mins,min(subs,[],1));
        maxes = max(maxes,max(subs,[],1));
        zs = unique(subs(:,3));
        for z = 1:length(zs)
            zSub = subs(subs(:,3)==zs(z),[1 2]);
            mid = median(zSub,1);
            zCents{zs(z)} = cat(1,zCents{zs(z)},[i mid]);
        end
    end
end


%% Align
shifts = zeros(numZ,2);
clear shifted
shifted{1} = zCents{1}(2:3);
for z = 2:length(zCents);
    
    p1 = zCents{z-1};
    p2 = zCents{z};
    inBoth = intersect(p1(:,1),p2(:,1));
    if ~isempty(inBoth)
        y1 = zeros(length(inBoth),1);
        y2 = y1; x1 = y1; x2 = y1;
        for p = 1:length(inBoth)
            t1 = find(p1(:,1)==inBoth(p));
            t2 = find(p2(:,1)==inBoth(p));
            y1(p,1) = p1(t1,2);
            y2(p,1) = p2(t2,2);
            x1(p,1) = p1(t1,3);
            x2(p,1) = p2(t2,3);
        end
    
        shiftY = median(y1-y2);
        shiftX = median(x1-x2);
        shifts(z,:) = shifts(z-1,:) + [shiftY shiftX];
        shifted{z} = p2(:,2:3) + shifts(z,:);
    
    end
end

oldP = [];
newP = [];
for z = 1:length(zCents)
    oldP = cat(1,oldP,[zCents{z}(:,2:3) ones(size(zCents{z},1),1) * z]);
    newP = cat(1,newP,[shifted{z} ones(size(shifted{z},1),1) * z]);
end

clf
hold on
scatter3(oldP(:,1),oldP(:,2),oldP(:,3),'r','.');
scatter3(newP(:,1),newP(:,2),newP(:,3),'k','.');

shiftZ.type = 'translation';
shiftZ.shifts = shifts;
save([MPN 'shiftZ.mat'],'shiftZ');







