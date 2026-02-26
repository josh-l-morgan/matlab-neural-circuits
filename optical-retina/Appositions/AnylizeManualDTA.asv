
%%Compare observation and events for manually counted Pre/Post synaptic
%%events

%{
Data Structure (dat)
1= first ob, 2 = last ob,	
3 = first dot, 4 = last dot,	
5 = first pre,	6 = last pre,	
7 = NumTimes,	file	cell	id


%}


Gains=dat(:,1)-dat(:,3); % time after start that puncta appeared
AfterGain = dat(:,2) - dat(:,3); % number of observations of puncta present
Losses = dat(:,2) - dat(:,4); % last time puncta observed prior to end
BeforeLoss = dat(:,1) - dat(:,4); % number of observations of puncta present

pGains=dat(:,5)-dat(:,3);  % time after first puncta that pre appeared
%pAfterGain = dat(:,6) - dat(:,3) + 1; % number of observations of pre present
pLosses = dat(:,6) - dat(:,4) ; % last time pre observed prior to loss of puncta
%pBeforeLoss = dat(:,4) - dat(:,5); % number of observations of pre present

%% Count observations
%%Gb = gain before
chkB=[-6:0]; % Check Before
chkA = [0:6];
clear SumGa SumGb SumLa SumLb pSumGa pSumGb pSumLa pSumLb
for i = 1:length(chkB)
   SumGb(i)=sum((Gains<=chkB(i)).*dat(:,7));
   SumGa(i)=sum((AfterGain(Gains<0)>=chkA(i)).*dat(Gains<0,7));
   SumLa(i)=sum((Losses(Losses>0)>=chkA(i)).*dat(Losses>0,7));
   SumLb(i)=sum((BeforeLoss(Losses>0)<=chkB(i)).*dat(Losses>0,7)); 
   
   pSumGb(i)=sum((pGains(Gains<0)<=chkB(i)).*dat(Gains<0,7));
   pSumGa(i)=sum((pGains(Gains<0)<=chkA(i)).*dat(Gains<0,7));
   pSumLa(i)=sum((pLosses(Losses>0)>=chkA(i)).*dat(Losses>0,7));
   pSumLb(i)=sum((pLosses(Losses>0)>=chkB(i)).*dat(Losses>0,7)); 
      
end


%% Fix numbers
TotGain=sum(dat(Gains<0,7)); 
TotLoss=sum(dat(Losses>0,7));
pSumLb=SumLb-(TotLoss-pSumLb);
pSumGa=SumGa-(TotGain-pSumGa);

SumGo=[SumGb(1:6) SumGa]
SumGp=[pSumGb(1:6) pSumGa]
SumLo=[SumLb(1:6) SumLa]
SumLp=[pSumLb(1:6) pSumLa]
Times=[chkB(1:6) chkA]

AllDat=[SumGp ; SumGo ; SumLp ; SumLo; Times]

%% Plot
plot(Times,SumGp./SumGo,'g')
hold on
plot(Times,SumLp./SumLo,'r')
ylim([0 1.1])
hold off


%% find Puncta Lifetimes
%%What does the survival curve of new puncta look like.


DotLifeOb = zeros(1,10); 
DotLifeSurv= zeros(1,10);
for i = 1: size(dat,1)
   sDat = dat(i,:);
   if sDat(1)~=sDat(3)  % if there is gain
       DotLifeOb(1:sDat(2)-sDat(3)+1)= DotLifeOb(1:sDat(2)-sDat(3)+1) + sDat(7);
       DotLifeSurv(1:sDat(4)-sDat(3)+1) = DotLifeSurv(1:sDat(4)-sDat(3)+1) +sDat(7);
   end
end

PrecentSurv = DotLifeSurv./DotLifeOb * 100;

plot(DotLifeOb/max(DotLifeOb),'r')
hold on
plot(DotLifeSurv/max(DotLifeSurv),'g')
plot(DotLifeSurv./DotLifeOb,'b')
hold off
    



res = [PrecentSurv ; DotLifeSurv; DotLifeOb]



    
