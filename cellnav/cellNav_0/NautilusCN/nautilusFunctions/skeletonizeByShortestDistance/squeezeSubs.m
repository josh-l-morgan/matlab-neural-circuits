function[moveObj maxSubs] = squeezeSubs(subs);

minSubs = min(subs,[],1);
for i = 1:3
   moveObj(:,i) = subs(:,i)-minSubs(i)+1; 
end
maxSubs = max(moveObj,[],1);

