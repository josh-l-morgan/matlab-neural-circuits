

clear all
load('MPN.mat')
load([MPN 'obI.mat']);
mot = getMotifs(obI);
synStruct = getSynMat(obI);
%anaMot = analyzeMotifs(obI);

%% find *t experiments
names = obI.nameProps.names;
tag = '125 part ';
clear nameTag tagged
c = 0;
for i = 1:length(names)
    nam = names{i};
    hit = regexp(nam, tag);
    if sum(hit) & (nam(end) == 't')
        c = c+1;
        nameTag{c} = nam;
        tagged(c) = i;
    end
end

%% 
allEdges = obI.nameProps.edges;
tri = mot.tri;
numTri = size(tri.triCell,1);

is125 = (allEdges(:,1)==125) | (allEdges(:,2) == 125);
tags125 = allEdges(is125,3);
anchors = obI.colStruc.anchors;
dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
anchors(:,1) = anchors(:,1)*dSamp(1);
anchors(:,2) = anchors(:,2)*dSamp(2);
anchors(:,3) = anchors(:,3)*dSamp(3);
anchors = round(anchors);
anchors(anchors<1) = 1;
% scatter(anchors(tags125,2),anchors(tags125,1),'.','k')
% hold on
%scatter(anchors(tagged,2),anchors(tagged,1),'s','r')


for i = 1:length(tagged)
    nameTag{i};
    targ = find(allEdges(:,3) == tagged(i));
    edge = allEdges(targ,1:2);
    anch = anchors(tagged(i),:);
    tNum = 0;  % number of triads found
    triDists = [];
    for t = 1:numTri % run the triad
        
        hit = (tri.triCell(t,1) == edge(2)) & ...
            (tri.triCell(t,2) == edge(1));
        if hit %% if it is the right RGC to LIN identity
            tNum = tNum+1;
            
            %% find nearest RGC in
            prePos = tri.synPos{t,1};
            
%             scatter(prePos(:,2),prePos(:,1),'o','r')
%             scatter(anch(2),anch(1),'x','g')
%             pause(.01)
            %hold off
            
            dist = sqrt((prePos(:,1)-anch(1)).^2 + ...
                (prePos(:,2)-anch(2)).^2 + ...
                (prePos(:,3)-anch(3)).^2);
            minDist = min(dist);
            sameIn = find(dist==minDist,1);
            pos12 = prePos(sameIn,:);
        
        
            %% find nearest LIN to TC
            
            prePos = tri.synPos{t,3};
            dist = sqrt((prePos(:,1)-pos12(1)).^2 + ...
                (prePos(:,2)-pos12(2)).^2 + ...
                (prePos(:,3)-pos12(3)).^2);
            minDist = min(dist);
            sameIn = find(dist==minDist,1);
            pos23 = prePos(sameIn,:);
            
            
            
            %% find nearest RGC to TC
        
            prePos = tri.synPos{t,2};
            dist = sqrt((prePos(:,1)-pos12(1)).^2 + ...
                (prePos(:,2)-pos12(2)).^2 + ...
                (prePos(:,3)-pos12(3)).^2);
            minDist = min(dist);
            sameIn = find(dist==minDist,1);
            pos13 = prePos(sameIn,:);
            
        
            %% Get triad distances
            eshift = tri.edgeOrder;    
            triOrd = '12 - 13, 12 - 23, 13 - 23';
            
            dist12_13 = sqrt((pos12(1)-pos13(1))^2 + ...
                (pos12(2)-pos13(2))^2 + ...
                (pos12(3)-pos13(3))^2);
            
            dist12_23 = sqrt((pos12(1)-pos23(1))^2 + ...
                (pos12(2)-pos23(2))^2 + ...
                (pos12(3)-pos23(3))^2);
            
            dist13_23 = sqrt((pos13(1)-pos23(1))^2 + ...
                (pos13(2)-pos23(2))^2 + ...
                (pos13(3)-pos23(3))^2);
            
            triDists(tNum,:) = [dist12_13 dist12_23 dist13_23] * .2;
            
                 
        end
        
    end
    
    if tNum
        triTagged(i).is = 1;
        triTagged(i).triDists = triDists;
    else
        triTagged(i).is = 0;
    end
    
    
end

hold off


%% Pick

allMax = 5; %Maximum distance of link for triad

taggedNum = length(triTagged);
anyTri = sum([triTagged.is]);

clear allTri closeTri closest
for i = 1:length(triTagged)
    
    
   if triTagged(i).is
   dists = triTagged(i).triDists
   passThresh = sum(sum(dists<=allMax,2)==3,1);
   passThresh = sum((dists(:,1)<=allMax),1);
   
   allTri(i) = size(dists,1);
   closeTri(i) = sum(passThresh);
   closest(i) = min(max(dists,[],2),[],1);
   closestLinRun(i) = min(dists(:,2));
   
   else
       allTri(i) = 0;
       closeTri(i) = 0;
       closest(i) = inf;
       closestLinRun(i) = inf;
   end
    
end


%% show closest
histRange = 0:6;
sampNum = length(allTri);
histAll = hist(allTri,histRange)/sampNum;
histClose = hist(closeTri,histRange)/sampNum;

bar(histRange,[histAll;histClose]')

bar(histRange,histClose,'k')


ce = closest(closest<10000);
closeRange = [0:2:20 300];
histClosest = histc(ce,closeRange)/length(ce);
bar(closeRange,histClosest,'k')
xlim([-5 25])


%% Show closest LIN run
histRange = 0:6;
sampNum = length(allTri);
histAll = hist(allTri,histRange)/sampNum;
histClose = hist(closeTri,histRange)/sampNum;

bar(histRange,[histAll;histClose]')

bar(histRange,histClose,'k')


ce = closestLinRun(closestLinRun<10000);
closeRange = [0:1:50  1000];
histClosest = histc(ce,closeRange)/length(ce);
bar(closeRange,histClosest,'k')
xlim([-5 60])






%% plot for all triads
isTri = cat(1,triTagged.is);
 allTriDist = cat(1,triTagged(isTri>0).triDists);
 
 useTri = (allTriDist(:,1)< 5) & (allTriDist(:,2)<1000)
 useDists = allTriDist(useTri,2);
 

ce = useDists;
closeRange = [0:1:600];
histClosest = histc(ce,closeRange)/length(ce);
bar(closeRange,histClosest,'k')
xlim([-5 300])



















