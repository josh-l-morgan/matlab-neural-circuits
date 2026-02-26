
%% Run a range of Radii and SDs
RadRs=.5:.01:1;
StanDevRs=FindDevs(.1);
NumTrials=size(RadRs,2)*size(StanDevRs,2);

clear Rdend Rdot
c=0;
clear RecRad RecSD
for r = 1:size(RadRs,2);
    for s = 1: size(StanDevRs,2);
        c=c+1;
        ['Running R=' num2str(RadRs(r)) ', S=' num2str(StanDevRs(s))];
        [Rdend(r,s) Rdot(r,s)]=MoGradF(RadRs(r),StanDevRs(s));
        RecRad(r,s)=RadRs(r);
        RecSD(r,s)=StanDevRs(s);
        PercentDone=c/NumTrials * 100
    end
end
% 
% Use=Rdend<10 & Rdot<10;
% Rdend=Rdend(Use);
% Rdot=Rdot(Use);
% RecRad=RecRad(Use);
% RecSD=RecSD(Use);

%% display Results

%%Relationship between Gaussian and dendritic ratio
scatter(RecSD(:),Rdend(:))

for r = 1 : size(RadRs,2)
    RadRs(r)
   scatter(Rdend(r,:),Rdot(r,:)),
   ylim([1 4])
   xlim([1 4])
   pause%(.1) 
end
% 
% for s = 1 : size(StanDevRs,2)
%    Rdend(s)
%    scatter(RecRad(:,s),Rdot(:,s)),
%    ylim([1 4])
%    pause%(.1) 
% end


MoGradRes.RecRad=RecRad;
MoGradRes.RecSD=RecSD;
MoGradRes.Rdend=Rdend;
MoGradRes.Rdot=Rdot;

scatter(Rdend, Rdot)
scatter(RecRad,Rdot)
scatter(RecSD,Rdot)