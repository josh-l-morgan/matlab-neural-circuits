function [ picked ] = pick( listIn,pickNum )
%pick 

if ~exist('pickNum')
    pickNum = 1;
end

if isa(listIn,'cell') 
    picked = listIn{fix(rand(pickNum,1) * length(listIn))+1};
else
    picked = listIn(fix(rand(pickNum,1) * length(listIn))+1);

end

 