function[] = alignSegmentation02()
%%Find translations for improving alignment
%%between each plane of segmentations stored as dsObj

fig = figure; %Make temporary figure

clear zHas zCents shifted shifts

global tis glob
if 1
    load('MPN.mat')  %Get volume directory
    load([MPN '\dsObj.mat']) %load o x n x 3 list of object subs
end

subsMax = double(max(cat(1,dsObj.subs),[],1)); %Get max of subs
subsMin = double(min(cat(1,dsObj.subs),[],1)); %Get min of subs

%%Make list of which objects are in each plane = zHas
%
% numZ = 480;
% for i = 1:numZ
%     zCents{i} = [];
% end
zCents{1} = []; % Make array to store object centers in

numObs = length(dsObj); %Find the total number of objets
mins = [inf inf inf]; % Make matrix for tracking miniums
maxes = [0 0 0]; % make matrix for tracking maxes
for i = 1:numObs %Start running objects
    subs = double(dsObj(i).subs); %Get subs for each object
    if ~isempty(subs) %if there are subs
        % mins = min(mins,min(subs,[],1)); % find mins of object
        % maxes = max(maxes,max(subs,[],1)); % find maxes of object
        zs = unique(subs(:,3)); %Make list of all planes that object is in
        for z = 1:length(zs) %Step through list of planes
            zSub = subs(subs(:,3)==zs(z),[1 2]); %Get subs for current plane
            mid = median(zSub,1); %Find median position of subs in current plane
            if zs(z)>length(zCents) % if current plane is larger then the size of the zCents array, add new position to the array
                zCents{zs(z)} = [];
            end
            zCents{zs(z)} = cat(1,zCents{zs(z)},[i mid]); % Add new center point to the plane that records [object id, y, x]
        end
    end
end

numZ = length(zCents); % count the number of planes spanned by objects

%% Find shifts
shifts = zeros(numZ,2); % Create matrix for recording final shifts
clear shifted
if ~isempty(zCents{1}) %if the first plane has a shift record
    shifted{1} = zCents{1}(2:3); % record shift in shifted matrix
else
    shifted{1} = [];
end

clear tforms
tforms(1).tform = fitgeotform2d([0 0;0 1;1 0], [0 0;0 1;1 0],'similarity');
for z = 2:length(zCents);    %Run through all sections using z to and compairng to z-1
    tforms(z).tform = tforms(1).tform;
    p1 = zCents{z-1}; %Get previous list of centers (z-1)
    p2 = zCents{z}; %Get list of centers for current plane
    if ~isempty(p1) & ~isempty(~p2) % If points exist in both planes
        [inBoth idxa idxb] = intersect(p1(:,1),p2(:,1)); % Find objects common to both planes
        if ~isempty(inBoth) %if the same objects are in both planes
            fixedPts = p1(idxa,1:2);
            movePts = p2(idxb,1:2);
            tforms(z).tform = fitgeotform2d(movePts, fixedPts,'similarity');
        end
    end
end

%% Propogate shifts
if 0 %Brute force test
    % newCents = zCents;
    %
    % for z = 1:length(tforms)
    %     tform = tforms(z).tform;
    %     for n = z:length(tforms)
    %         if ~isempty(newCents{n})
    %             [nY nX] = transformPointsForward(tform,newCents{n}(:,2),newCents{n}(:,3));
    %             newCents{n}(:,2) = nY;
    %             newCents{n}(:,3) = nX;
    %         end
    %     end
    % end
    %
    % oldP = [];
    % newP = [];
    % for z = 1:length(tforms)
    %     if ~isempty(newCents{z})
    %         currentPts = newCents{z}(:,2:3);
    %         newP = cat(1,newP,[currentPts ones(size(currentPts,1),1) * z]);
    %         currentPts = zCents{z}(:,2:3);
    %         oldP = cat(1,oldP,[currentPts ones(size(currentPts,1),1) * z]);
    %     end
    % end

else
    for z = 2:length(tforms)
        previous = tforms(z-1).tform;
        current = tforms(z).tform;
        new = current;
        new.A = previous.A * current.A;
        new.Scale = new.A(1,1);
        new.Translation = [new.A(1,3) new.A(2,3)];
        tforms(z).tform = new;
    end

    %%Remove negative
    maxZ = subsMax(3);
    minY = subsMin(1);
    minX = subsMin(2);
    cornerY = [subsMin(1); subsMin(1);subsMax(1);subsMax(1)];
    cornerX =[subsMin(2); subsMax(2);subsMax(2);subsMin(2)];
    newCorners = [];
    for z = 1:maxZ
        tform = tforms(z).tform;
        [nY nX] = transformPointsForward(tforms(z).tform,cornerY,cornerX);
        newCorners = cat(1,newCorners,[nY nX nY*0+z]);
        minX = min(minX,min(nX));
        minY = min(minY,min(nY));
    end

    shiftTform  = tforms(1).tform;
    shiftTform.A(1,3) = -minY+10;
    shiftTform.A(2,3) = -minX+10;
    for z = 1:length(tforms)
        tforms(z).tform.A = tforms(z).tform.A * shiftTform.A;    
    end


    %%Translate points for visualization
    oldP = [];
    newP = [];
    for z = 1:length(tforms)
        if ~isempty(zCents{z})
            currentPts = zCents{z}(:,2:3);
            [nY nX] = transformPointsForward(tforms(z).tform,currentPts(:,1),currentPts(:,2));
            oldP = cat(1,oldP,[currentPts ones(length(nY),1) * z]);
            newP = cat(1,newP,[nY nX ones(length(nY),1) * z]);
        end
    end

end

%% Make minimum shift 0
clf
hold on
scatter3(oldP(:,1),oldP(:,2),oldP(:,3),'k','.');
scatter3(newP(:,1),newP(:,2),newP(:,3),'r','.');
scatter3(newCorners(:,1),newCorners(:,2),newCorners(:,3),'g','.')

shiftZ.type = 'similarity';
shiftZ.shifts = shifts;
shiftZ.tforms = tforms;
save([MPN 'shiftZ.mat'],'shiftZ');







