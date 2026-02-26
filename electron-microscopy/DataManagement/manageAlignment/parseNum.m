function[numList] = parseNum(nam,tags)


%%

tagNum = 1;
getOn = 0;
gotNum = char;
startNum = 0;
numList = [];
for i = 1:length(nam)
    if isempty(tagNum)
        getOn = 1;
    elseif strcmp(nam(i:i+length(tags{tagNum})-1),tags{tagNum})
        getOn = 1;
    end
        
    if getOn
    d = str2num(nam(i));
    if ~isempty(d) 
        startNum = 1;
        gotNum(length(gotNum)+1) = nam(i);
    end
    
    if startNum & isempty(d)
        startNum = 0;
        numList(tagNum) = str2num(gotNum);
        tagNum = tagNum+1;
        gotNum = char;
    end
    end
    
    if tagNum>length(tags)
        break
    end
end


