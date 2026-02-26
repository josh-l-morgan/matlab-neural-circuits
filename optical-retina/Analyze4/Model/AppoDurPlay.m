clear all


TotTime=1000;
MeanStart=TotTime/4;
MeanDur=300;
DurStd=100;


Num=10000;

% Form=randn(Num,1)*TotTime/10+MeanStart;
Form=rand(Num,1)*TotTime;
hist(Form)
Dur=randn(Num,1)*DurStd+MeanDur;
% Dur=ones(Num,1)*500;
Dur(Dur<0)=0;
hist(Dur,0:max(Dur)/20:max(Dur))
Lost=Form+Dur;
hist(Lost)

NumOn=zeros(TotTime,1);
for i = 1:TotTime
   On=Form<=i & Lost>i;
   NumOn(i)=sum(On);    
end
bar(NumOn)


%% Probe
StartProbe=TotTime/2;
NumProbe=15;
Int=10;
Times=StartProbe+(0:NumProbe-1)*Int;

for i = 1:NumProbe
    Time=Times(i);
    ProbeOn(:,i)=Form<=Time & Lost>Time;
end

%% Stats
Starters=ProbeOn(ProbeOn(:,1),:);
NumStart=mean(ProbeOn(:,1));
mStarters=mean(Starters,1);
plot(mStarters) 
ylim([0 1])


%% Translate probability of elimination / time into mean duration
%%Prob at time = proportion of time

TimeProbed=(1:NumProbe-1)*Int;
GuessDur=TimeProbed./(1-mStarters(2:NumProbe))




%% conclusions
%{
Short intervals are subject to noise.  However, the more durations there are 
that are shorter then the interval, the more more the prediction of average
duration will be biased by the shape of the duration distribution.

%}

