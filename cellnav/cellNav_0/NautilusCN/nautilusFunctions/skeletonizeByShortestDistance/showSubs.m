function[maxObj objVol] = showSubs(subs,prop);

if ~exist('prop')
    prop = ones(1,size(subs,1))*300;
end


minSubs = min(subs,[],1);
for i = 1:3
   moveObj(:,i) = subs(:,i)-minSubs(i)+1; 
end
maxSubs = max(moveObj,[],1);

objVol = zeros(maxSubs,'uint8');
objVol(sub2ind(maxSubs,moveObj(:,1),moveObj(:,2),moveObj(:,3))) = prop;
maxObj = max(objVol,[],3);
image(maxObj), pause(.01)