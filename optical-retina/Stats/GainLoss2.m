
%{

input: dat
1=TetPos, 2=Age, 3=hours, 4 = min, 
5=length, 6= dots, 7 = gain, 8 = loss

%}

Pos=dat(:,1)>0;
Time=dat(:,3)+dat(:,4)/60;
Age = dat(:,2);
Length = dat(:,5)/1000;

RateGain=dat(:,7)./dat(:,5)./Time * 1000; % gain/mm/h
RateLoss=dat(:,8)./dat(:,5)./Time * 1000; % loss/mm/h
ProbLoss=dat(:,8)./dat(:,6)./Time * 1000;

DD=dat(:,6)./dat(:,5); %densities

Turn=(dat(:,7)+dat(:,8))./(dat(:,6)*2);
Turn24=(Turn./Time)*24;
PlotSE(Pos,Turn24)


scatter(Pos,RateGain)
xlim([-1 3])
ranksum(RateGain(Pos),RateGain(~Pos))

scatter(Pos,RateLoss)

scatter(Pos,ProbLoss)
ranksum(ProbLoss(Pos),ProbLoss(~Pos))


scatter(Pos,DD)
ranksum(DD(Pos),DD(~Pos))

%% Split Experiments

subplot(1,3,1)
scatter(Pos(Age==8),DD(Age==8),'r'), hold on
scatter(Pos(Age==9),DD(Age==9),'g'), 
scatter(Pos(Age==10),DD(Age==10),'b'), hold off
xlim([-1 2])

subplot(1,3,2)
scatter(Pos(Age==8),RateGain(Age==8),'r'), hold on
scatter(Pos(Age==9),RateGain(Age==9),'g'), 
scatter(Pos(Age==10),RateGain(Age==10),'b'), hold off
xlim([-1 2])
% 
% scatter(Pos(Age==8),RateLoss(Age==8),'r')
% scatter(Pos(Age==10),RateLoss(Age==10),'b')

subplot(1,3,3)
scatter(Pos(Age==8),ProbLoss(Age==8),'r'), hold on
scatter(Pos(Age==9),ProbLoss(Age==9),'g'), 
scatter(Pos(Age==10),ProbLoss(Age==10),'b'), hold off
xlim([-1 2])


subplot(1,3,1)
PlotSE(Pos,DD)

subplot(1,3,2)
PlotSE(Pos,RateGain)

subplot(1,3,3)
% PlotSE(Pos,RateLoss)
% PlotSE(Pos,RateGain-RateLoss)
PlotSE(Pos,ProbLoss)



ranksum(ProbLoss(Pos),ProbLoss(~Pos))

ranksum(RateGain(Pos & Age==8),RateGain(~Pos& Age==8))
ranksum(RateGain(Pos & Age==10),RateGain(~Pos& Age==10))

ranksum(DD(Pos),DD(~Pos))

%% Anova
% x = [Pos Age];
% [d p stats] = manova1(x,RateGain);
% 

%% Plot




%% Predict P21

AddDays = 12;

for p = 0:1
   meanRateGain(p+1) = mean(RateGain(Pos==p));
   meanProbLoss(p+1) = mean(dat(Pos==p,8)./dat(Pos==p,6));
   meanDD(p+1)=mean(DD(Pos==p)); 
   meanDnum(p+1) = mean(dat(Pos==p,6));
   meanGain(p+1) = mean(dat(Pos==p,7));
   meanLength(p+1) = mean(dat(Pos==p,5));
end

newAge=9;
newDnum=meanDnum;
newRdnum=dat(:,6);
for i = 1: AddDays * 12
    newAge(i+1)=newAge(i)+1/12;
    newDnum(i+1,:)=newDnum(i,:);
    %newDnum(i+1,:) = newDnum(i+1,:) - newDnum(i,:) .* meanProbLoss;
    newDnum(i+1,:)=newDnum(i+1,:) + meanGain;
    newRdnum(:,i+1) = newRdnum(:,i)+ RateGain .* Length * 2;
    
end
newAge


for i = 1:size(newRdnum,1)
    newRdd(i,:) = newRdnum(i,:)/Length(i);
    if Pos(i), col = 'r'; else, col = 'k'; end
    col
    plot(newAge,newRdd,col); 
    hold on
end
hold off

plot(newAge,mean(newRdd(~Pos,:),1),'k');
hold on
plot(newAge,mean(newRdd(Pos,:),1),'r')
hold off




PlotSE(Pos,newRdd(:,size(newRdd,2)))








