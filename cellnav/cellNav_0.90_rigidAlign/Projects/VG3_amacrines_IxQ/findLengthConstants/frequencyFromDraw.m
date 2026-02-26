
%If you draw from a population of frequencies with an equal distribution
%between 0 and 1, 99% of the time that you pull a 0, the source frequency
%will be <= X

%N = 319;
N = 196;

freqs = [0:.001:1];
reps = 10000
hRange = [0:N];
h = zeros(length(freqs),length(hRange));
for f = 1:length(freqs)
    disp(sprintf('testing frequency %0.2f (%d of %d)',freqs(f),f,length(freqs)))
    r = rand(reps,N);
    hits = sum(r<=freqs(f),2);
    h(f,:) = hist(hits,hRange)/reps;
    
end

colormap jet(100)
image(h*1000)
plot(freqs,h(:,1))

hNorm1 = h ./ repmat(sum(h,1),[size(h,1) 1]);
cumH1 = cumsum(hNorm1,1);
cumH2 = cumsum(hNorm1,1,'reverse');
thresh = 0.975;
threshH1 = cumH1>thresh;
threshH2 = cumH2>thresh;
lookupFreq = repmat(freqs',[1 length(hRange)]);
lookupFreq2 = lookupFreq;
lookupFreq2(threshH2==0)= 0;
lookupFreq1 = lookupFreq;
lookupFreq1(threshH1==0) = inf;
highBound = min(lookupFreq1,[],1);
lowBound = max(lookupFreq2,[],1);

%%median
threshHm = cumH1> 0.5;
lookupFreqM = lookupFreq;
lookupFreqM(threshHm) = 0;
medBound = max(lookupFreqM,[],1);


colormap(jet(100))
image(hNorm1 *400/max(hNorm1(:)))
hold on
plot(hRange,lowBound*1000,'w')
hold on
plot(hRange,highBound*1000,'w')
hold off

drawFreq.h = h;
drawFreq.N = N;
drawFreq.freqs = freqs;
drawFreq.hRange = hRange;


save('.\drawFreq.mat','drawFreq')

%%  
%obs = [0 1 2 4 5 6 19 24 28 41 83 130]
obs = [13 69 48 25 17 18 6];

obsH = h(:,obs+1);

for i = 1:length(obs)
   plot(freqs,obsH(:,i))
   hold on
    
end
hold off

res = cat(2,obs',obs'/N*100,medBound(obs+1)'*100,lowBound(obs+1)'*100,highBound(obs+1)'*100)



