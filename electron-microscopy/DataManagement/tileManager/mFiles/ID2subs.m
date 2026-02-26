function[subs] = ID2subs(IDs);

holder = 10^12;
w = fix(IDs/10^9)-1000;
s = mod(fix(IDs/10^6),1000);
r = mod(fix(IDs/10^3),1000);
c = mod(IDs,1000);

subs = [w s r c];