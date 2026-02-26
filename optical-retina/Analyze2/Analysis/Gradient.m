

TPN=GetMyDir
load([TPN 'Dots.mat'])


%% Identify cell center with user's help
viewDots=zeros(Dots.ImSize,'uint8');
for i = 1:Dots.Num
    viewDots(Dots.Vox(i).Ind)=1;
   
end

colormap gray(255)
MaxVD=sum(viewDots,3);
image(MaxVD*300),pause(.1)

'Select Center of Cell Body'
[x y] = ginput(1)

OrthoVD=sum(viewDots,1);
OrthoVD=squeeze(OrthoVD);
image(OrthoVD'*300),pause(.1)

'Select Center of Cell Body'
[a z] = ginput(1)

Dots.CellCenter.Pos=[y x z];


clear viewDots

%% Find distance of each dot to cell center

for i = 1: Dots.Num
    DistAll=Dots.CellCenter.Pos-Dots.Pos(i,:);
    DistAll(1:2)=DistAll(1:2)*.103;
    DistAll(3)=DistAll(3)*.3;
    Dots.CellCenter.DistAll(i,:)=DistAll;
    Dots.CellCenter.XYZDist(i)=sqrt(DistAll(1)^2 + DistAll(2)^2 + DistAll(3)^2);
    Dots.CellCenter.XYDist(i)=sqrt(DistAll(1)^2 + DistAll(2)^2); 
end
    
%% Find stats with sliding window    
XYDist=Dots.CellCenter.XYDist;
%%Grab all dots that Passed Crit
Pass=Dots.Pass;
Win=10; %define window size
d=0;
clear NumDots Volumes ItSums MeanBrights DFOfs DFOfTopHalfs DFs Dists   AllBrights
%AllBrights=zeros(1,1000);
for i = 0 :5: max(XYDist)
    d=d+1;
    %Select puncta 
    Include= (XYDist < (i+Win/2)) & (XYDist > (i-Win/2)) & Pass;
    if sum(Include(:)>0)
        IDs=find(Include);
        Dists(d)=i;
        Area= pi * (i+Win/2)^2 -  (pi*(max(0,(i-Win/2)))^2);
        NumDots(d)=sum(Include(:))/Area;
        Volumes(d)=mean(Dots.Vol(Include));
        ItSums(d)=mean(Dots.ItSum(Include));
        MeanBrights(d) =mean(Dots.MeanBright(Include));
        DFOfs(d)=mean(Dots.DFOf(Include));
        DFOfTopHalfs(d) = mean(Dots.DFOfTopHalf(Include));
        DFs(d)=mean(Dots.DF(Include));
    end
end


subplot(2,1,1)

viewDots=zeros(Dots.ImSize,'uint8');
Show=find(Pass);
for i = 1:size(Show,2)
    viewDots(Dots.Vox(Show(i)).Ind)=Dots.MeanBright(Show(i));
   
end

colormap gray(255)
MaxVD=sum(viewDots,3);
image(MaxVD*300),pause(.1)


subplot(2,1,2)
plot(Dists,NumDots)
plot(Dists,Volumes)
plot(Dists,ItSums)
plot(Dists,MeanBrights)
plot(Dists,DFOfs)
plot(Dists,DFOfTopHalfs)
plot(Dists,DFs)


PlotStats=Dots.MeanBright;
scatter(Dots.CellCenter.XYDist,PlotStats,'.')
hold on
scatter(Dots.CellCenter.XYDist(Pass),PlotStats(Pass),'.','r')
hold off
    

    
    
    
    
    
    
    
    