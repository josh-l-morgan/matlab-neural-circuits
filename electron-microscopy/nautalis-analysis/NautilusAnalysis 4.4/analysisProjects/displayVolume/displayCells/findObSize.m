function[fsize] = findObSize(subStruc)


fsize = [0 0 0 ];

for i = 1:length(subStruc)
   subs = subStruc(i).subs;
   maxSubs = max(subs,[],1);
   fsize = max([fsize; maxSubs],[],1);    
end