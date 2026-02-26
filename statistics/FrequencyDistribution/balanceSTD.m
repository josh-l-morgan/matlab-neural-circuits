

%3.476

std1 = 10
n = 10000;

dat = rand(n,1) * std1 * 2 - std1;

stdr = std(dat)

cRange = [3:.001:4];
stdr = cRange * 0;
for i = 1:length(cRange)
    
    c = cRange(i);

    dat = (rand(n,1)-.5) * std1 * c;
    stdr(i) = std(dat);


end

stde = std1-stdr;
plot(stde)

cRange(abs(stde) == min(abs(stde)))


%% bimodal

% for s = 1, c = 0.867, .... for s = 2, c = 0; .... for s = 1.7, c = 0.527
std1 = 10
s = 1.5
bn = 1000000;
m1 = 10;

clf
cRange = [0.6:.001:.71];
stdr = cRange * 0;
hRange = [m1-std1*5:std1/10:m1+std1*5];
for i = 1:length(cRange)
    
    c = cRange(i);

 dat = randn(bn,1)*std1*c + (rand(bn,1)>.5) * s *std1 - s *std1/2 + m1;
    stdr(i) = std(dat);
    stdr(i) 
    h = hist(dat,hRange);
    plot(h)
    drawnow



end

stde = std1-stdr;
plot(stde)

cRange(abs(stde) == min(abs(stde)))




