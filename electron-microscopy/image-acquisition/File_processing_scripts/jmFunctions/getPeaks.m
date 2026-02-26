function[pos vals] = getPeaks(vec)

vec = vec(:);
dif = vec(2:end) - vec(1:end-1);
peaks = ~xor(dif(2:end)>0, dif(1:end-1)<0); 
peaks = [0; peaks];
pos = find(peaks);
vals = vec(pos);