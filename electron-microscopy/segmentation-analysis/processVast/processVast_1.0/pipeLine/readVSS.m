

fileName = 'C:\Users\joshm\Documents\myWork\LGNs1\jlmHomeSeg\S8_joshmHome_14+01+17C.vss'


f = fileread(fileName);
fclose(fileName)

returns = cat(2,[0],regexp(f,'\r'));


for i = 1:length(returns)-1
    lines{i} = f(returns(i)+1:returns(i+1));
end
