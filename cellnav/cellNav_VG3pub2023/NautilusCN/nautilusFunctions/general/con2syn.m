function[syn] = con2syn(con);

%%Convert connectivity matrix to list of synapses

%%

[y x] = find(con>0);
vals = con(:);
vals = vals(vals>0);

syn = zeros(sum(vals),2);
c = 0;
for i = 1:length(y)
    syn(c+1:c+vals(i),1) = y(i);
    syn(c+1:c+vals(i),2) = x(i);
    c = c+ vals(i);
end