



%% Get data
mot = getMotifs(obI);
synMat = getSynMat(obI);
synStruct = getSynMat(obI);
synPos = getSynPos(1);


%% use synmat to get post tcr
is125pre = synMat.pre == 125;
post125class = synMat.postClass(is125pre)

synMat.typeNames

[sum(post125class==1) ...
sum(post125class==2) ...
sum(post125class==3) ...
sum(post125class==4)]

is125post = synMat.post == 125;
pre125class = synMat.preClass(is125post);

[sum(pre125class==1) ...
sum(pre125class==2) ...
sum(pre125class==3) ...
sum(pre125class==4)]


%% use syn pos to get post tcr

fromRGC_dat = dsAnchors(synPos.postRGCPos,obI,[2 1 3]);
fromLIN_dat = dsAnchors(synPos.postLinPos,obI,[2 1 3]);


%% motifs

%anaMot = analyzeMotifs(obI);

rgcs = mot.cel.types.rgcs;
tcrs = mot.cel.types.tcrs;
lins = mot.cel.types.lins;

%% syns
syns = synMat;
from125 = synMat.pre == 125;
to125 = synMat.post == 125;
toRGC  = from125 & (synMat.postClass == 1);
toTCR  = from125 &(synMat.postClass== 2);
toLIN = from125 & (synMat.postClass == 3);
toUNK = from125 & (synMat.postClass == 4);
fromRGC = to125 & (synMat.preClass == 1);
fromLIN = to125 & (synMat.preClass == 3);
fromUNK = to125 & (synMat.preClass == 4);

sum(toRGC)
sum(toTCR)
sum(toLIN)
sum(toUNK)

sum(fromRGC)
sum(fromLIN)
sum(fromUNK)

synTarg = synMat.pre(to125);
synTargRGC = intersect(synTarg,rgcs);
synTargTCR = intersect(synTarg,tcrs);
synTargLIN = intersect(synTarg,lins);

synPos125toTCR = cat(1,synMat.synPos(toTCR,:));
synPos125toLIN = cat(1,synMat.synPos(toLIN,:));

synPos125fromRGC = cat(1,mot.syn.synPos{fromRGC,:});
synPos125fromLIN = cat(1,mot.syn.synPos{fromLIN,:});
synPos125fromUNK = cat(1,mot.syn.synPos{fromUNK,:});

%% triads
tri = mot.tri;
triCell = tri.triCell;
triCell = triCell(triCell(:,1) == 125,:);
triTarg = triCell(:,3);
triTargRGC = intersect(triTarg,rgcs);
triTargTCR = intersect(triTarg,tcrs);
triTargLIN = intersect(triTarg,lins);


prim125 = find(tri.triCell(:,1)==125);
useTri = prim125;
%for u = 1:length(prim125);%1:size(tri.triCell,1);

triPoints = drawTriads(tri,useTri);

%synCol = cat(1,triPoints.lineCol,triPoints.ballCol);
% synRad = cat(1,triPoints.lineRad,triPoints.ballRad);
% synType = cat(1,triPoints.lineRad*0+1,triPoints.ballRad*0+2);
% useSyn = cat(2,triPoints.lineGroup, triPoints.ballGroup);
primCell = unique(triPoints.cellGroup(:,1));
secCell = unique(triPoints.cellGroup(:,2));
tertCell = unique(triPoints.cellGroup(:,3));



%% diads

diads = mot.di.diCell;
diads = diads(diads(:,1) == 125,:);
diTarg = diads(:,3);
diTarg = setdiff(diTarg,triTarg);

diTargRGC = intersect(diTarg,mot.cel.types.rgcs);
diTargTCR = intersect(diTarg,mot.cel.types.tcrs);
diTargLIN = intersect(diTarg,lins);


%% color relationship to 124

postTCR = intersect(tcrs,isPost);
triTarg = intersect(tertCell,tcrs);
diTarg = setdiff(diTarg,triTarg);



