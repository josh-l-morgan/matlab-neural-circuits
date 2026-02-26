function[y x z] = find3(volT,num)

if ~exist('num','var'),num = inf;end
ind = find(volT);
num = min(num,length(ind));
[y x z] = ind2sub(size(volT),ind(1:num));
