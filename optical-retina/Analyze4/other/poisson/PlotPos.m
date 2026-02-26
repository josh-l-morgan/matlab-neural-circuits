function[AutoAve,Difs]=PlotPos(Raw)

colormap gray(100)
[Reps,Dur]=size(Raw);

%% Generate Raw Data

subplot(3,1,1)
image(~Raw*100)

%% Plot Autocorrelation
[R,T]=find(Raw);
Add=zeros(Reps,Dur*2-1);
Num=Add;
for i = 1:length(T)
   Add(R(i),Dur-T(i)+1:2*Dur-T(i))= Add(R(i),Dur-T(i)+1:2*Dur-T(i))+Raw(R(i),:); 
   Num(R(i),Dur-T(i)+1:2*Dur-T(i))=Num(R(i),Dur-T(i)+1:2*Dur-T(i))+1;
end
Auto=Add./Num;
AutoAve=mean(Auto,1);
subplot(3,1,2)
plot(AutoAve)


%% Plot interval distribution
Difs=[];
for i = 1:Reps
   T=find(Raw(i,:)); 
   Dif=T(2:length(T))-T(1:length(T)-1);
   Difs(length(Difs)+1:length(Difs)+length(Dif))=Dif;
    
end

subplot(3,1,3)
hist(Difs,1:1:max(Difs))

