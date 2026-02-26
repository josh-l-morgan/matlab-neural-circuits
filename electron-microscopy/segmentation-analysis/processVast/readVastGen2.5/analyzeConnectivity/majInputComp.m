
%% load data
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
TPN = [MPN 'skel\'];
load([TPN 'cellArbors'])
load([MPN 'obI.mat'])

%% Get syn data and select synapses
targCell = 108;
synapses = obI.nameProps.edges;
preWithTarg = unique(synapses(synapses(:,1)==targCell,2));
preWithTarg = preWithTarg(preWithTarg>0);
synWithPre = [];
for i = 1:length(preWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==preWithTarg(i)) & ...
        (synapses(:,1)>0)));
end
synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);
synPos = ceil(obI.colStruc.anchors(synObj,:)/8);

conArbor.targCell  = targCell;


%% find closest point and total length within 2x dist of closest point

preList = unique(synPre);
preList = preList(preList>0);
postList = unique(synPost);
postList = postList(postList>0);
lookDist = 30;

goodCheck = synObj*0 + 1;
showPlots = 0;
synMat = zeros(length(preList),length(postList));
lengthMat = zeros(length(preList),length(postList),lookDist);

for pre = 1 : length(preList)
    for post = 1: length(postList)
        synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
    end
end



%% Find majors

synCounts = synMat(synMat>0);
hist(synCounts,[1:1:max(synCounts)])
majThresh = 2 * median(synCounts);
percentMaj = 100 * sum(synCounts>majThresh)/length(synCounts);
majMat = synMat>=majThresh;

%% show input distributions

for i = 1:size(synMat,2)
    subplot(size(synMat,2),1,i)
    subplot(1,1,1);
    countSyn = synMat(:,i);
    countSyn = sort(countSyn(countSyn>0),'descend');
    bar(countSyn);
    pause
    
end



    count108 = synMat(:,3);
    count108 = sort(count108(count108>0),'descend');
    hist108 = hist(count108,[1:4:20]);
    hist108 = hist108/max(hist108);
    bar(count108);
    sharedMat108 = getMinShared(inMat);
    
    
    blankMat = synMat;
    blankMat(:,3) = 0;
    countBlanked = blankMat(:);
    countBlanked = sort(countBlanked(countBlanked>0),'descend');
    histBlanked = hist(countBlanked,[1:4:20]);
    histBlanked = histBlanked/max(histBlanked);
    bar(countBlanked)

    bar([histBlanked' hist108'])











