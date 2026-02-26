colormap gray(100)


%% Generate Raw Data
P=0.1;
Dur=1000;
Reps=10;
Rp=rand(Reps,Dur);



%% Run Data
BurstProb=.01;
BurstDur=50;
Burst=rand<=BurstProb*BurstDur;
for r = 1:Reps
    for d= 1:Dur
        if rand<=BurstProb
           Burst=BurstDur+randn*BurstDur;
           Pb=P+(1-P)*rand;
        end
        Burst=Burst-1;
        if Burst>0
            Pu=Pb;
        else
            Pu=P;
        end
        Raw(r,d)=Rp(r,d)<=Pu;
        
    end
end

[AutoAve,Difs]=PlotPos(Raw);