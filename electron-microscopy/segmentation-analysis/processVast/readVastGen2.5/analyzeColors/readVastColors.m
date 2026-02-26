function[objInfo] = readVastColor(fileName);


if ~exist('fileName')
[TFN TPN] = uigetfile('.txt');
fileName = [TPN TFN];
end

f = fileread(fileName);

returns = cat(2,[0],regexp(f,'\r'));


for i = 1:length(returns)-1
    lines{i} = f(returns(i)+1:returns(i+1));
end

header = f(1:returns(1));


%%
matCount = 0;

for i = 1:length(lines)
    lin = lines{i};
    [a,count,errmsg,nextindex]=sscanf(lin, '%d   %d    %d %d %d %d   %d %d %d %d   %d %d %d   %d %d %d %d   %d   %d %d %d %d %d %d   ');
    
    if ~isempty(a)
        
    matCount = a(1);

    if (matCount>0) 

      data(matCount,:)=a';
      
      anchors(matCount,1) = a(12)';
            anchors(matCount,2) = a(11)';
      anchors(matCount,3) = a(13)';

      col1(matCount,:) = a(3:5);
      col2(matCount,:) = a(7:9);
      boundBox(matCount,:) = a(19:24);
    quotes = regexp(lin,'"');
    comments = regexp(lin,'%');
    
    if isempty(comments) & ~isempty(quotes)
        spaces = regexp(lin,' ');
        ids{matCount,1} = str2num(lin(1:spaces(1)));
        objName = lin(quotes(end-1)+1:quotes(end)-1);
        objName = objName(objName ~= '*');
        ids{matCount,2} = objName;
    end
    end
    end
end

%%

for i = 1:size(ids,1)
    idNum = ids{i,1};
    if idNum
        objNames{idNum} = ids{i,2};
    end
end

objInfo.lines = lines;
objInfo.names = objNames;
objInfo.ids = ids;
objInfo.anchors = anchors;
objInfo.boundBox = boundBox;
objInfo.col1 = col1;
objInfo.col2 = col2;


