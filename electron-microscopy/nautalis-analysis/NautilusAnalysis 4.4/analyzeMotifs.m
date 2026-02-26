function[anaMot] = analyzeMotifs(obI);


mot = getMotifs(obI);
syn = getSynMat(obI);

%% Questions
%{

1) Filter motifes by euclidean distance

How many RGC inputs that are not triads

Possible motifs = 
    classic triad = RGC/TC + RGC/LIN + LIN/TC
    simple rgc input = RGC/LIN
    simple tc output = LIN/TC
    simple lin to lin = LIN/LIN
    reciprocal = LIN1/LIN2 + LIN2/LIN1
    rgc driven reciprocal = RGC/LIN1 + RGC/LIN2 + LIN1/LIN2 + LIN2/LIN1
    reciprocal driving TC = LIN1/TC + LIN2/TC + LIN1/LIN2 + LIN2/LIN1

What is the distance distribution from inputs to outputs?

What is the influence of each input over each output?

What is the probability that the LIN will innervate the same TCs as an RGC
input

What is the overlap of the output field (minus local triad?)


Globally promiscuous, ultrastructurally specifi, microciruit specificity?

Get distance between synapses converging on same target. 

Get motifs for 

%}


%% For each RGC input, does it form a triad and, if so, how far away
syn;
toRGC = syn.post == 125;
preIsRGC  = syn.preClass == 1;
synUse = find(toRGC & preIsRGC);


%% Check for triad
clear tri
maxDist = 5;
%triCell = mot.tri.triCell;
post125 = syn.post((syn.pre==125)& (syn.post>0));
targ2 = 125;
triad.maxDist = 5;
for i = 1:length(synUse) %for each rgc input
    targ = synUse(i);
    pos12 = syn.synPos(targ,:);
    ID = syn.obID(targ);
    targ1 = syn.pre(targ);
    
    post1 = unique(syn.post((syn.pre==targ1)& (syn.post>0)));
    targs3 = intersect(post1,post125);
    targs3 = setdiff(targs3,targ2);
    
    
    tri(i).targ1 = targ1;
    tri(i).targ2 = targ2;

    numTri = length(targs3);
    tri(i).numTri = numTri;
    targ2Dist = [];
    for t = 1:numTri;
        targ3 = targs3(t);
        pos13 = syn.synPos((syn.pre == targ1) & ((syn.post) == targ3),:);
        pos23 = syn.synPos((syn.pre == targ2) & ((syn.post) == targ3),:);
        
        dists = sqrt((pos13(:,1)-pos12(1)).^2 + (pos13(:,2)-pos12(2)).^2 ...
            + (pos13(:,3)-pos12(3)).^2);
        d = find(dists == min(dists));
        pos13 = pos13(d,:);
        dist13 = min(dists);
        
        dists = sqrt((pos23(:,1)-pos12(1)).^2 + (pos23(:,2)-pos12(2)).^2 ...
            + (pos23(:,3)-pos12(3)).^2);
        d = find(dists == min(dists));
        pos23 = pos23(d,:);
        dist23 = min(dists);
                
        tri(i).closest(t).targ3 = targ3;
        tri(i).closest(t).pos13 = pos13;
        tri(i).closest(t).pos23 = pos23;
        tri(i).closest(t).dist13 = dist13;
        tri(i).closest(t).dist23 = dist23;
        tri(i).targ2Dist(t) = dist23;
    end
    tri(i).numClose = sum(tri(i).targ2Dist<=maxDist);
    if numTri == 0
        tri(i).minDist = inf;
    else
        tri(i).minDist = min(tri(i).targ2Dist);
    end
end

triad.tri = tri;

%%
closeList = cat(1,tri.numClose);
triList = cat(1,tri.numTri);
minDists = cat(1,tri.minDist);
rgcMat = [triList closeList minDists]



tri(triList==0).targ1































