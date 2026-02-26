%function[] = anaCl(TPN DPN)
%% Find clustering
%Concider both absolute and dendritic distances. 
'start'

clear all
%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions
ClusterDist=5;

TPN=GetMyDir;

%% load segments
'loading data'

%{
%% Get Mask

load([TPN '\data\D.mat']) %Load dendrite mask
%sumD=sum(D,3);
%maskD=sumD>0;
%ID=maskD;

[Dy Dx Dz]=find3(D); % vectorize by find 3
Dy=Dy*xyum; Dx=Dx*xyum; Dz=Dz*zum; %scale to microns
Dv=[Dy Dx Dz]; clear Dy Dx Dz % Concatinate D find results into vector

clear D
%}

load([TPN '\data\AllSegCut.mat'])
AllSeg=AllSegCut; clear AllSegCut

Mids=mean(AllSeg,3);  %find midpoints
for i=1:size(Mids,1)  %get segment lengths
   Length(i)=dist(AllSeg(i,:,1),AllSeg(i,:,2)); 
end
Nodes=[AllSeg(:,:,1) ; AllSeg(:,:,2) ; Mids]; %list all nodes (tips and midpoints)

%% Get all OK dot positions in um
load([TPN '\Dots.mat'])
DP=Dots.Pos;
DP(:,1:2)=DP(:,1:2)*.1;
DP(:,3)=DP(:,3)*.3;
DP=DP(Dots.OK,:);


%% create Nearest Node list (nearest node for each dot and distance)
[ASy ASx ASz] = size(Mids);
Doted=zeros(ASy,1);
NN=zeros(size(DP,1),3);
DotToNN=zeros(size(DP,1),1);
for i = 1:size(DP,1)
    Ndist=dist(Mids,DP(i,:)); %find dist from dot to all midpoints
    Near=min(Ndist); %find shortest distance
    [NearSeg NearEnd] =find(Ndist==Near,1); %get node at that distance
    NN(i,:)=Mids(NearSeg,:); %assign that node to NearestNode list for dots
    Doted(NearSeg) = Doted(NearSeg)+1;  %Recored via index that seg has dot
    DotToNN(i)=Near; %record that distance for posterity    
end



%% Start Clutering
%%Important variables Length, Mids, Doted, AllSeg
mpDots=zeros(size(Mids,1),1);
mpLength=zeros(size(Mids,1),1);
for i = 1 : size(Mids,1)
   clear Branch
   Check=ones(size(Mids,1),1);  %create index of checked
   Branch(1).tip=AllSeg(i,:,1);  %Create structure array for growing skel
   Branch(1).length=Length(i)/2;
   Branch(1).done=0;
   Branch(1).dots=Doted(i)/2;
   Branch(2).tip=AllSeg(i,:,2);
   Branch(2).length=Length(i)/2;
   Branch(2).done=0;
   Branch(2).dots=Doted(i)/2;
   Check(i)=0; %Clear starting point
   while Branch(size(Branch,2)).done==0; %run while the last branch remains undone. 
   for p= 1:size(Branch,2)
       while (Branch(p).done==0) %run until done
          clear NewPoints
          %Generate list of next points 
          np=0;
          next=find(AllSeg(:,1,1)==Branch(p).tip(1) & AllSeg(:,2,1)==Branch(p).tip(2) & AllSeg(:,3,1)==Branch(p).tip(3) & Check);
          for n=1:size(next,1), np=np+1; NewPoints(np,1:3)=(AllSeg(next(n),:,2)); npIndex(np)=next(n);
            Check(next(n))=0; 
          end %clear segment from check list
           SN=size(next);
          next=find(AllSeg(:,1,2)==Branch(p).tip(1) & AllSeg(:,2,2)==Branch(p).tip(2) & AllSeg(:,3,2)==Branch(p).tip(3) & Check);
          SN=SN+size(next);
          for n=1:size(next,1), np=np+1; NewPoints(np,1:3)=(AllSeg(next(n),:,1)); npIndex(np)=next(n); 
            Check(next(n))=0; 
          end %clear segment from check list
           SN;np;

            %%Enter first New Point
          if ~exist('NewPoints') %if there were no new points 
              Branch(p).done=1; %finish
          else %if there is something new
              Branch(p).tip=NewPoints(1,:); %add new tip
              Branch(p).length=Branch(p).length+Length(npIndex(1));  %increace branch length
              Branch(p).dots=Doted(npIndex(1)); %add dots
              if Branch(p).length > ClusterDist, Branch(p).done=1; end %finish if branch is long enough
              
              %%Form new branch if necessary
              if size(NewPoints,1)>1 %if Branch Branched
                for nb=2:size(NewPoints,1)
                    NewNum=size(Branch,2)+1; %position for new branch
                    Branch(NewNum).tip=NewPoints(nb,:);
                    Branch(NewNum).length=Length(npIndex(nb));
                    Branch(NewNum).done=0;
                    Branch(NewNum).dots=Doted(nb);
                end %ran all new branches
              end %if there was a branching
              
          end %if there were new points
          
        end   %while branch isnt done
    
    end  %run Branches
    end %while the last branch is not done
    
    %%Gather midpoint stats
    mpBranches=size(Branch,2);
    for b=1:size(Branch,2)
        mpDots(i)=mpDots(i)+Branch(b).dots;
        mpLength(i)=mpLength(i) + Branch(b).length;
    end
      
    
    if mod(i,100)==0
        PercentDoneClimbing=double(i)/size(Mids,1)*100
    end
    
end %run all midpoints

%%find Density
mpDD=mpDots./mpLength;

%%Draw density map
MapSize=fix(max(Mids))+1;
DDmap=zeros(MapSize(1),MapSize(2));
Lmap=DDmap; Dmap=DDmap;
for i=1:size(Mids,1)
   DDmap(fix(Mids(i,1))+1,fix(Mids(i,2))+1)=mpDD(i);
   Dmap(fix(Mids(i,1))+1,fix(Mids(i,2))+1)=mpDots(i);
   Lmap(fix(Mids(i,1))+1,fix(Mids(i,2))+1)=mpLength(i);
end

%%%%(((((((((((((((((((((((((
%%look for branches
for i = 1:size(AllSeg,1)
    Found1=findv(AllSeg,AllSeg(i,:,1));
    SizeFound(i)=size(Found1,1);
end
hist(SizeFound)


for i = 1 :size(Mids,1)

next=find(AllSeg(:,1,1)==AllSeg(i,1,1) & AllSeg(:,2,1)==AllSeg(i,2,1) & AllSeg(:,3,1)==AllSeg(i,3,1) );
sizeNext(i)=size(next,1);
end
hist(sizeNext)


%%Current Progress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



pause

%% run all mid points
'Running Segments'
MidStep=1;
Mid2=Mids(1:MidStep:size(Mids,1),:);

for i = 1: size(Mid2,1)  
    DotDist=dist(Dots,Mid2(i,:));
    MidDist=dist(Mids,Mid2(i,:));
    NNDist=dist(NN,Mid2(i,:));
    %VolDist=dist(Dv,Mids(i,:));
  
        
    %% Run all Bins
    for b = 1:size(DBins,2)
        LocalDots(i,b)=sum(DotDist<=DBins(b));
        LocalLength(i,b)=sum(Length(MidDist<=DBins(b)));
        LocalN(i,b)=sum(NNDist<=DBins(b));
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    if mod(i,100)==0, PercentDone=i/size(Mid2,1)*100,end
end
LocalDD=LocalDots./LocalLength;
LocalND=LocalN./LocalLength;
%% run from prespective of dots

'Running Segments'
for i = 1: size(Dots,1)  
    DotToDotDist=dist(Dots,Dots(i,:));
    DotToMidDist=dist(Mids,Dots(i,:));
    %VolDist=dist(Dv,Mids(i,:));
  
        
    %% Run all Bins
    for b = 1:size(DBins,2)
        DotToDots(i,b)=sum(DotToDotDist<=DBins(b));
        DotToLength(i,b)=sum(Length(DotToMidDist<=DBins(b)));
        DotToLength(i,b)=max(DotToLength(i,b),DBins(b)); %set minimum length
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    if mod(i,100)==0, PercentDone=i/size(Dots,1)*100,end
end



DotToDD=DotToDots./DotToLength;

%}



%% Calculate Volumes for each Bin

'Checking Volume'
VolumeOn=0;
if VolumeOn == 1

    for i = 1: size(Mid2,1)  
        VolDist=dist(Dv,Mid2(i,:));        
        %% Run all Bins
        for b = 1:size(DBins,2)
            LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
        end
        if mod(i,10)==0, DendPercentDone=i/size(Mid2,1)*100,end
    end



%save('.\temp\LocalVol.mat','LocalVol')

%%Get SurfaceArea
%%assuming cylinders
%a=sqrt(Vol/(L*pi))*L*2*pi;
LocalA=sqrt(LocalVol./(LocalLength*pi)).* (LocalLength * pi * 2);
minArea=LocalLength * pi *2 * .05;  % conciders minimum area for length as 0.1um diam process
LocalA(LocalA<minArea)=minArea(LocalA<minArea); %area is at least minimum for tiny process
LocalDA=LocalDots./LocalA;


%%Run Volume for Dots

    for i = 1: size(Dots,1)  
        VolDist=dist(Dv,Dots(i,:));        
        %% Run all Bins
        for b = 1:size(DBins,2)
            DotToVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
        end
        if mod(i,10)==0, DotPercentDone=i/size(Dots,1)*100,end
    end



%save('.\temp\LocalVol.mat','LocalVol')

%%Get SurfaceArea
%%assuming cylinders
%a=sqrt(Vol/(L*pi))*L*2*pi;
DotToA=sqrt(DotToVol./(DotToLength*pi)).* (DotToLength * pi * 2);
minArea=DotToLength * pi *2 * .05;  % conciders minimum area for length as 0.1um diam process
DotToA(DotToA<minArea)=minArea(DotToA<minArea); %area is at least minimum for tiny process
DotToDA=DotToDots./DotToA;




Local.Mid.LocalVol=LocalVol;
Local.Mid.LocalA=LocalA;
Local.Mid.LocalDA=LocalDA;

Local.Dot.DotToVol=DotToVol;
Local.Dot.DotToA=DotToA;
Local.Dot.DotToDA=DotToA;

end %if Volume on


%% Get some useful numbers
hist(LocalDD(:,3),0:.02:3)
pause(1)
hist(DotToDD(:,3),0:.05:3)
hist(DotToDots(:,3),0:.02:10)

for i=1:size(DBins,2)-1
   diffAll(:,i)=abs(LocalDD(:,i)-LocalDD(:,i+1)); 
end
MeanDiffs=mean(diffAll,1)

Local.DBins=DBins;
Local.Mid.MidStep = MidStep;

Local.Dot.DotToDots=DotToDots;
Local.Dot.DotToLength=DotToLength;
Local.Dot.DotToDD=DotToDD;


Local.Mid.LocalDots=LocalDots;
Local.Mid.LocalLength=LocalLength;
Local.Mid.LocalDD=LocalDD;

%%Nearest Neighbor percent
clear Neighbors

for b= 1: 10
    for n=1:10
        Neighbors(b,n)=sum(DotToDots(:,b)>n)/size(DotToDots,1)*100;
    end
    plot(Neighbors(b,:)),hold on
end
hold off 
Local.Dot.Neighbors=Neighbors;


%save([TPN '\data\Local.mat'],'Local')

%% Display

LocalPlot=LocalDD; %decide who to plot

for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldLD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)=LocalDots(i);
    end        
end
fieldC=fieldLD*0; %make field count
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)= fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)+1;
    end        
end
fieldDD=fieldLD*0;
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)=fieldDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)+LocalPlot(i,b)/fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b);
    end        
end


%%Make Dot map
fieldDP=zeros(fix(max(Dots(:,1)))+2,fix(max(Dots(:,2))+2));
for i = 1: size(Dots,1)
   fieldDP(round(Dots(i,1))+1,round(Dots(i,2))+1)=fieldDP(round(Dots(i,1))+1,round(Dots(i,2))+1)+1;
end
mask=max(fieldDD,[],3)>0;
xs=max(size(mask,2),size(fieldDP,2));
ys=max(size(mask,1),size(fieldDP,2));
CellPic=zeros(ys,xs,3);
CellPic(1:size(mask,1),1:size(mask,2),1)=mask;
CellPic(1:size(fieldDP,1),1:size(fieldDP,2))=fieldDP;
CellPic(:,:,1)=mask*20;
CellPic(:,:,3)=CellPic(:,:,3)*75;
CellPic=uint8(CellPic);
subplot(6,6,[1:5,7:11,13:17,19:23,25:29])
image(CellPic*10)



%% Play movie
r=1:255;b=255:-1:1;
r=r/255; b=b/255;
rb=r'; rb(:,3)=b';
rb(2:size(rb,1)+1,:)=rb;
rb(1,:)=0;


mask=max(fieldDD,[],3)>0;
image(mask*100)
load('.\temp/cmap.mat')
set(gcf,'Colormap',cmap)
pause(.01)
clear field

meanBright=mean(LocalPlot(:,20));
Times=3;

while 1==1
for b = 1:size(DBins,2)
    Nine=sort(LocalPlot(:,b));
    Ninety=Nine(round(.95*size(Nine,1)));
    scaleDD=255/(meanBright*Times); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    field=fieldDD(:,:,b);
    field=field*scaleDD;
    field(mask)=field(mask)+2;
    subplot(6,6,[1:5,7:11,13:17,19:23,25:29])
    
    image(field)
    subplot(6,6,31:35)
    Diam=([1:size(field,2)]<=DBins(b)*2)*1000;
    Diam(:,:,2)=Diam*255;     Diam(:,:,3)=Diam(:,:,1)*255;
    image(uint8(Diam))
    subplot(6,12,[12,24,36,48])
    Inten=1:256;
    image(Inten')
    htext = uicontrol('Style','text','String',Times);
    DBins(b)
    pause
end
end




'finished'
%%Notes
%{
Could segment each bin
individual variation in density as bins increase
montecarlo for nearest node
run for surface area (and save)
how good is area assumption for multiple diameters

%}

