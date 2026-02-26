function[ids] = readVastColor(fileName);


if ~exist('fileName')
[TFN TPN] = uigetfile('.txt');
fileName = [TPN TFN];
end

f = fileread(fileName);

returns = cat(2,[0],regexp(f,'\r'));


for i = 1:length(returns)-1
    lines{i} = f(returns(i)+1:returns(i+1));
end

matCount = 0;
for i = 1:length(lines)
    lin = lines{i};
    quotes = regexp(lin,'"');
    comments = regexp(lin,'%');

    
    if isempty(comments) & ~isempty(quotes)
        matCount = matCount + 1;
        spaces = regexp(lin,' ');
        ids{matCount,1} = str2num(lin(1:spaces(1)));
        ids{matCount,2} = lin(quotes(end-1)+1:quotes(end)-1);
    end
    
end


