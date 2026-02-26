





s = 10; %Standard deviation
m = 100; %Mean
z = 1.96;
f = .5; %Fraction of a standard deviation


Nr = (z/f)^2;
N = ceil(Nr)



%%%
reps = 100000;
samps = randn(N,reps)*s+m;
means = mean(samps,1);
dif = abs(means-m);
goodFrac = mean(dif<=(s*f));









