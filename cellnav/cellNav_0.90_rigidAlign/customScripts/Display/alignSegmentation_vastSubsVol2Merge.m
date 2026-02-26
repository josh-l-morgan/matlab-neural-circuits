function[shiftZ] = alignSegmentation_vastSubsVol2Merge(vastSubs)
%%Find translations for improving alignment
%%between each plane of segmentations stored as dsObj
if 1
    fig = figure; %Make temporary figure

    clear zHas zCents shifted shifts

    global tis glob
    % if 1
    %     load('MPN.mat')  %Get volume directory
    %     %load([MPN '\dsObj.mat']) %load o x n x 3 list of object subs        
    %     load([MPN '\vastSubs.mat']) %load o x n x 3 list of object subs
    % 
    % end

    allSubs = vastSubs{:};
    subsMax = double(max(allSubs,[],1)); %Get max of subs
    subsMin = double(min(allSubs,[],1)); %Get min of subs
    maxZ = subsMax(3);

    %%Make list of which objects are in each plane = zHas
    %
    % numZ = 480;
    % for i = 1:numZ
    %     zCents{i} = [];
    % end
    zCents{1} = []; % Make array to store object centers in

    numObs = length(vastSubs); %Find the total number of objets
    mins = [inf inf inf]; % Make matrix for tracking miniums
    maxes = [0 0 0]; % make matrix for tracking maxes
    for i = 1:numObs %Start running objects
        disp(sprintf('Running object %d of %d',i,numObs))
        subs = double(vastSubs{i}); %Get subs for each object
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

    %numZ = length(zCents); % count the number of planes spanned by objects
end

%% Find shifts
clear tformsRaw
showRes = 0;
blankT = rigidPtFit();
tformsRaw(1).tform = blankT;
%tformsRaw(1).tform = fitgeotform2d([0 0;0 1;1 0], [0 0;0 1;1 0],'similarity');
%tformsRaw(1).tform  = pcregistericp([0 0;0 1;1 0], [0 0;0 1;1 0]);
for z = 2:length(zCents);    %Run through all sections using z to and compairng to z-1
    disp(sprintf('Finding alignment for section %d of %d',z,length(zCents)));
    tformsRaw(z).tform = tformsRaw(1).tform;
    p1 = zCents{z-1}; %Get previous list of centers (z-1)
    p2 = zCents{z}; %Get list of centers for current plane
    if ~isempty(p1) & ~isempty(~p2) % If points exist in both planes
        [inBoth idxa idxb] = intersect(p1(:,1),p2(:,1)); % Find objects common to both planes
        if ~isempty(inBoth) %if the same objects are in both planes
            fixedPts = p1(idxa,2:3);
            movePts = p2(idxb,2:3);
            %tformsRaw(z).tform = pcregistericp(movePts,fixedPts);
            %tformsRaw(z).tform = fitgeotform2d(movePts, fixedPts,'similarity');
            tformsRaw(z).tform = rigidPtFit(movePts,fixedPts, showRes);
        end
    end
end


%% Select shifts
tformsFilt = tformsRaw;
performance = zeros(length(tformsFilt),1);
for z = 1:length(tformsFilt)
    performance(z) = tformsFilt(z).tform.betterFrac;
    if performance(z) < 0
        tformsFilt(z).tform = blankT;
    else
        %tformsFilt(z).tform.A(1:2,1:2) = [1 0 ; 0 1];
    end
end


%% Propogate shifts
%Brute force test
newCents = zCents;

for z = length(tformsFilt):-1:2
    tform = tformsFilt(z).tform;
    for n = length(tformsFilt):-1:z
        if ~isempty(zCents{n})
            newPts = newCents{n}(:,2:3);
            newPts =  rigidPts(newPts,tform.A);
            % newPts = newPts + tform.T;
            % newPts =  newPts * tform.R;
            newCents{n}(:,2:3) = newPts;
        end
    end
end

oldP = [];
bruteP = [];
for z = 1:length(tformsFilt)
    if ~isempty(newCents{z})
        currentPtsN = newCents{z};
        bruteP = cat(1,bruteP,[currentPtsN ones(size(currentPtsN,1),1) * z]);
        currentPtsO = zCents{z};
        oldP = cat(1,oldP,[currentPtsO ones(size(currentPtsO,1),1) * z]);
    end
end

%%Create blank
tforms = tformsFilt;
for z = 1:length(tforms)
    tforms(z).tform.n = 0;
    tforms(z).tform.T = [0 0];
    tforms(z).tform.R = [1 0;0 1];
end
%%Propogate transforms
for z = 2:length(tforms)
    previous = tforms(z-1).tform;
    current = tformsFilt(z).tform;
    new = current;
    new.A = previous.A * current.A;
    %new.A = current.A * previous.A;
    new.R = new.A(1:2,1:2);
    new.T = new.A(1:2,3)';
    tforms(z).tform = new;
end




%%Remove negative
if 0
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
end

%%Translate points for visualization
newP = [];
for z = 1:length(tforms)
    if ~isempty(zCents{z})
        newPts = zCents{z}(:,2:3);
        newPts =  rigidPts(newPts,tforms(z).tform.A);
        % newPts = newPts + tforms(z).tform.T;
        % newPts = newPts * tforms(z).tform.R;
        newP = cat(1,newP,[zCents{z}(:,1) newPts ones(size(newPts,1),1) * z]);
    end
end



%% Show results
clf
hold on
ax = gca;
ax.Clipping = 'off';
for r = 1:3;

    %useO = find(rand(numObs,1)<.1);
    useO = unique(ceil(rand(520,1)*numObs));
    [hit idxHit] = intersect(oldP(:,1),useO);

    idxHit = [];
    for o = 1:length(useO)
        idxHit = cat(1,idxHit,find(oldP(:,1)==useO(o)));
    end

    %idxHit = find(oldP(:,1) == ceil(rand*numObs));

    cla(ax);
    view(rand*360,00)
    ax.NextPlot = 'add';
    scatter3(oldP(idxHit,2),oldP(idxHit,3),oldP(idxHit,4),4,'v','MarkerFaceColor','none','markerFaceAlpha',.2,'markeredgecolor',[0 0 0],'markeredgealpha',.2)
    scatter3(newP(idxHit,2),newP(idxHit,3),newP(idxHit,4),4,'*','MarkerFaceColor','none','markerFaceAlpha',.2,'markeredgecolor',[1 0 0],'markeredgealpha',.2)
    scatter3(bruteP(idxHit,2),bruteP(idxHit,3),bruteP(idxHit,4),4,'o','MarkerFaceColor','none','markerFaceAlpha',.2,'markeredgecolor',[0 0 1],'markeredgealpha',.2)
    pause(1)
    %scatter3(newCorners(:,1),newCorners(:,2),newCorners(:,3),'g','.')
end


%% Save
shiftZ.type = 'rigid';
shiftZ.tforms = tforms;


As = zeros(length(shiftZ.tforms(z)),3,3);
for z = 1:length(shiftZ.tforms)
    As(z,:,:) = shiftZ.tforms(z).tform.A;
end

shiftZ.As = As;










