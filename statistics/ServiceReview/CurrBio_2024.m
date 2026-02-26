%% Boot strap intervals



%%Input data
dat.n1 = 52; % Incidence of binocular responces in interneruons
dat.v1 = 33;
dat.n2 = 49; % Incidence of binocular response in relay cells
dat.v2 = 20;

%%Model parameters
reps = 100000;

%%Pull data for model
n1 = dat.n1;
n2 = dat.n2;
v1 = dat.v1;
v2 = dat.v2;
nT = n1+n2;
vT = v1+v2;
vFrac = vT/nT;


%%Run model
v1r = zeros(reps,1);
v2r = zeros(reps,1);
for r = 1:reps
    v1r(r,1) = sum(rand(n1,1)<=vFrac);
    v2r(r,1) = sum(rand(n2,1)<=vFrac);
end

%%Convert model results percent difference
v1rp = v1r/n1;
v2rp = v2r/n2;
vDr = v1rp - v2rp;
vDr = sort(vDr);


%%Convert observed
v1p = v1/n1;
v2p = v2/n2;
vD = v1p-v2p;

%%Calculate stats
alpha = 0.05;
lowerB = vDr(round(reps*alpha/2));
uppperB = vDr(round(reps - (reps*alpha/2)));













