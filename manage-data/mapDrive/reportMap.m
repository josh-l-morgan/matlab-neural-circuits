function[report] = reportMap(dat);


nextPath = dat.nextPath;
paths = dat.paths(1:nextPath,1);
folders = dat.folders(1:nextPath,1);
totBytes = dat.totBytes(1:nextPath,1);
parent = dat.parent(1:nextPath,1);

pathCount = size(paths,1);
parentCount = histc(parent,[1:pathCount]);

firstPathLength = length(paths{1});

%%Sum to branch
sumBytes = totBytes;
sumCaps = dat.containsCappedPaths;
checkTip = ones(pathCount,1);
tips = find((parentCount==0) & checkTip);
childrenLeft = parentCount;
while sum(checkTip)

    for f = 1:size(tips,1)
        tp = tips(f);
        if parent(tp)
        sumBytes(parent(tp)) = sumBytes(parent(tp)) + sumBytes(tp);
        sumCaps(parent(tp)) = sumCaps(parent(tp)) + dat.containsCappedPaths(tp);
        childrenLeft(parent(tp)) = childrenLeft(parent(tp)) - 1;
        end
        checkTip(tp) = 0;        
    end
    tips = find((childrenLeft==0) & checkTip);
    tipLength = length(tips);
end

%%Sort
listReport = cell(size(paths,1),3);
sumBytesCappedBig = sumBytes + sumCaps * 10^15;
[dummySort sortByteIdx] = sort(sumBytesCappedBig,'descend');
sortBytes = sumBytes(sortByteIdx);
sortPaths = paths(sortByteIdx);
sortCapped = sumCaps(sortByteIdx);
for i = 1:size(sortPaths,1)
    listReport{i,1} = sprintf('%0.9f tb',sortBytes(i)/10^12);
    listReport{i,2} = sortCapped(i);
    listReport{i,4} = sortPaths{i};
    listReport{i,3} = sortPaths{i}(firstPathLength:end);
end

report.listReport = listReport;
report.sumBytes = sumBytes;
report.sumCaps = sumCaps;
report.sortByteIdx = sortByteIdx;




