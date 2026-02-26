clear all
SPN = 'D:\Joshm\S8\export_14+04+12\';


TPN = [SPN(1:end-1) '_mat\'];
if ~exist(TPN),mkdir(TPN),end

PPN = [TPN 'obPlanes\'];
if ~exist(PPN),mkdir(PPN),end

dSamp = 25;
Isize = [15000 15000];

for z = 100:100:10000


zName = sprintf('%05.0f.mat',z)

pass = 1;
try 
    load([PPN zName])

catch err
    pass = 0;
end

if pass
subs = [];
for s = 1:length(vastOb.subs)
    subs = cat(1,subs,vastOb.subs{s});
end
subs = ceil(subs/dSamp);
%maxSubs = max(subs,[],1);
maxSubs = Isize/dSamp;
showI = zeros(maxSubs(1:2));

inds = sub2ind(maxSubs(1:2),subs(:,1),subs(:,2));
uinds = unique(inds);
hinds = hist(inds,uinds);
showI(uinds) = hinds;

image(showI),pause(.01)
end
end

