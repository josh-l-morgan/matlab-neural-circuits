function[maxCon] = sumMaxCon(con);

%% gives percentage of synapses participating in maximum connectivity of each node (by row)

maxCon = sum(max(con,[],2));
maxCon = maxCon/sum(con(:))*100;