function [zerocross] = zerocrossing(input)
%estimates zeroscrosssing before peak

% if max(input) > abs(min(input)) %ON cells
if input(end-10) > 0 %ON cells
    ind = find(input == max(input));
    input(ind:end)=[];
    prior = find(input<0,1,'last');
% elseif max(input) < abs(min(input)) %OFF cells
elseif input(end-10) < 0 %OFF cells
    ind = find(input == min(input));
    input(ind:end)=[];
    prior = find(input> 0,1,'last');
    
else
    disp ('NO UNIQUE EXRTEMUM')
end

y1 = input(prior);
y2 = input(prior+1);
m = (y2-y1);
zerocross = (m*prior - y1)/m;

