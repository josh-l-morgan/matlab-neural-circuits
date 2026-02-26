function[preIn postIn] = parseBinDat(dat);

if ~exist('dat','var')
    load('dat.mat')
end

preIn = [];
postIn = [];

c = 0;

for i = 1: size(dat,1)
   cons = dat(i,1:3);
   numAx = dat(i,4);
   for n = 1:numAx
      c = c + 1;
       for d = 1:3
         L = length(preIn);
         numCon = cons(d);
         preIn(L+1:L + numCon) = c;
         postIn(L+1:L + numCon) = d;
       end
   end
end