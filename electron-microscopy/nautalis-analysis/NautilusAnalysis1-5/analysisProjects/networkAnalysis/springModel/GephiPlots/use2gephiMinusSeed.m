%function[csvText] = use2gephi(useList)


%% make node list
% 
% lineNum = 1
% colHead = 

useList = obI2cellList_minusSeed(obI,[108 201])

mat2gephDir = 'D:\LGNs1\Analysis\connectivitySpreadsheets\'
nodeFileName = [mat2gephDir 'nodesNoSeed.csv']
edgesFileName = [mat2gephDir 'edgesNoSeed.csv']

%%

lineNum = 1;
clear N
N{1,1} = 'ID';
N{1,2} = 'Class';
N{1,3} = 'Color';
N{1,4} = 'Label';

for i = 1:length(useList.postList)
    lineNum = lineNum+1;
    N{lineNum,1} = num2str(useList.postList(i));
    N{lineNum,2} = 'tcr';
    N{lineNum,3} = '[1 0 0]';
        N{lineNum,4} = num2str(useList.postList(i));

end

for i = 1:length(useList.preList)
    lineNum = lineNum+1;
    N{lineNum,1} = num2str(useList.preList(i));
    N{lineNum,2} = 'rgc';
    N{lineNum,3} = '[0 0 1]';
    N{lineNum,4} = num2str(useList.preList(i));

end



delete(nodeFileName)
[nrows ncols] = size(N);
fid = fopen(nodeFileName, 'w');
for row=1:nrows
    lineStr = char();
    for col = 1:ncols
        lineStr = [lineStr N{row,col} ','];
    end
    lineStr = lineStr(1:end-1);
    fprintf(fid, '%s\n',lineStr);
end
fclose(fid);



%%

lineNum = 1;
clear E
E{1,1} = 'source';
E{1,2} = 'target';
E{1,3} = 'type';
E{1,4} = 'weight';
E{1,5} = 'ID';
E{1,6} = 'Label';

for y = 1:size(useList.con,1)
    for x = 1:size(useList.con,2)
        if useList.con(y,x)>0
        
           lineNum = lineNum+1;
    E{lineNum,1} = num2str(useList.preList(y));
    E{lineNum,2} = num2str(useList.postList(x));
    E{lineNum,3} = 'Directed';
    E{lineNum,4} = num2str(useList.con(y,x));
        end
        
    end
end


delete(edgesFileName)
[nrows ncols] = size(E);
fid = fopen(edgesFileName, 'w');
for row=1:nrows
    lineStr = char();
    for col = 1:ncols
        lineStr = [lineStr E{row,col} ','];
    end
    lineStr = lineStr(1:end-1);
    fprintf(fid, '%s\n',lineStr);
end
fclose(fid);







