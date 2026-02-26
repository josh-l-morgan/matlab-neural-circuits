function[report] = sortMap(dat);

nextPath = dat.nextPath;
paths = dat.paths(1:nextPath,1);
totBytes = dat.totBytes(1:nextPath,1);

firstPathLength = length(dat.SPN);

sumBytes = totBytes;

%%Sort
listReport = cell(size(paths,1),3);
[sortBytes sortByteIdx] = sort(sumBytes,'descend');
sortPaths = paths(sortByteIdx);
for i = 1:size(sortPaths,1)
    listReport{i,1} = sprintf('%0.9f TB',sortBytes(i)/10^12);
    listReport{i,2} = sortPaths{i};
    listReport{i,3} = sortPaths{i}(firstPathLength:end);
end

report.listReport = listReport;



