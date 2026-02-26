function[prev dists] = passMax(con1,vProp,conDat,seed)


vNum = size(con2,1);

con2 = con1;
con2 = cat(1,con1(1,:)*0+1,con1+1); %shift to lookup.

dist1 = conDat.dists;

dists = ones(vNum+1,1) * vNum*2;
dists(seed+1) = 0;
prev = dists * 0;


%%

for r = 1:100000;
    tic
   for i = 1:26
       newVox = con2(:,i);
       newDists = dists(newVox)+dist1(i);
       closer = newDists < dists;
       dists(closer) = newDists(closer);
       prev(closer) = newVox(closer);
   end
   
   if ~sum(~prev)
       break
   end
   toc
end


