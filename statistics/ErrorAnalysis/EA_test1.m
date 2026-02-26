
m = 100; %average biological value of control
bv = 100; %normal, standard deviation of biological value
bvExp = 1; %exponent of biological variation
es = [10 100 500 1000]; % Experimental effect size
ev  = bv; %normal, standard deviation of experimental effect size
dp = .01; % average detection probability
dt = 0; % detection threshold. 
dv = 10; % normal, standard deviation of detection

db = 1; % detection background level
n = 100;


clear bio val

bio{1} = zeros(n,1) + m + (randn(n,1) * bv).^bvExp;
for i = 1:length(es)
    baseValues = zeros(n,1) + m 
    bio{1+i} = poissrnd(likelyDetected) zeros(n,1) + m + (randn(n,1) * bv).^bvExp + es(i) + randn(n,1) * ev;
end

clf

sp = subplot(2,1,1);
sp.NextPlot = 'add';
for i = 1:length(bio)
    swarmchart(sp,bio{i}*0+i,bio{i})
end
 


for i = 1:length(bio)
    likelyDetected = bio{i};
    likelyDetected(likelyDetected<dt) = 0; %applyl threshold
    likelyDetected = max(0,likelyDetected); %remove negatives
    likelyDetected = likelyDetected * dp; %apply average detection probability
    detected = poissrnd(likelyDetected);
    reading = detected/dp;
    reading(isnan(reading)) = 0;
    reading(reading<0) = 0;
    readNoise = randn(n,1) * dv;
    reading = db +  reading .* readNoise;
    val{i} = reading;

end

sp2 = subplot(2,1,2);
sp2.NextPlot = 'add';
for i = 1:length(val)
    swarmchart(sp2,val{i}*0+i,val{i})
end






