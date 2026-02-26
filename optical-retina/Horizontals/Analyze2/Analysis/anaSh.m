function[]=anaSh(TPN)
% finds distribution of dots in xy relative to central point

'Starting analysis of Dots and Dendrites'
if exist([TPN 'data\Results.mat'])
        load([TPN 'data\Results.mat'])
end
load([TPN 'data\DotStats'])

%% Get directory name
DPNd=[TPN 'data\'];
Cname=TPN;
%choose source data based on availability



%% Vectorize Dot Data
'Vectorizing Dot Data'
xyum=.103;
zum=.3;

load([DPNd 'DotStats.mat'])
BC=DotStats(:,:,3); %BC is the dot position in microns
BC=BC(DotStats(:,3,1)>.5,:); %dont count dots with low DFOf

save([DPNd 'BC.mat'],'BC')
TotalDots=size(BC,1);
save([DPNd 'TotalDots.mat'],'TotalDots')

%% Assign Dend properties
'Finding Dendrite Segment Properties'

%%Load Segments
if exist([TPN 'dataFix\AllSeg.mat'])
    load([TPN 'dataFix\AllSeg.mat'])
else
    load([TPN 'data\AllSeg.mat'])
end

%%Find Segment Lengths
SegLength=sqrt((AllSeg(:,1,1)-AllSeg(:,1,2)).^2 +(AllSeg(:,2,1)-AllSeg(:,2,2)).^2 +(AllSeg(:,3,1)-AllSeg(:,3,2)).^2);
save([DPNd 'SegLength.mat'],'SegLength')
TotalLength=sum(SegLength);
save([DPNd 'TotalLength.mat'],'TotalLength')
AverageDensity=TotalDots/TotalLength;
save([DPNd 'AverageDensity'], 'AverageDensity')

%%Find Segment Midpoints
SegMP=(AllSeg(:,:,1)+AllSeg(:,:,2))/2;

%% Find cell body center
if exist([TPN 'temp\Imax.mat'])
    load([TPN 'temp\Imax.mat'])
else
    DI=imread([TPN 'images\RawDendMax.tif']);
    Dsum=sum(sum(DI));
    Dt=find(Dsum==max(Dsum),3);
    Imax=DI(:,:,Dt);
end
    
load([DPNd 'Threshold.mat'])
Ithresh=Imax>Thresh;
image(Ithresh*300),pause(.01)
[ys,xs,zs]=size(Imax);


SE=strel('disk',3);
for i=1:100
    Ithresh=imerode(Ithresh,SE);
    [a,b]=bwlabel(Ithresh); 
    if b==1, break, end
end %run i
image(Ithresh*300)
[y x]=find(Ithresh);
ym=mean(y); xm=mean(x);

CB=[ym xm];  %Cell body yx position
Results.XY.CellBody(1,:)=CB;
CB=CB*xyum;  %Scale Cell body center to microns
Results.XY.CellBody(2,:)=CB;


%% Find arbor center
for i=1:3; Acent(i)=sum(SegMP(:,i) .* SegLength)/TotalLength; end



%% Bin with dist from center
'Binning Data with xy'  %concider depth restrictions

xyBin=5;
cent=CB; %Define central reference point

%%make list of pixel distances
c=0;
pd=zeros(xs*ys,1);
for x=1:xs, for y=1:ys
    c=c+1; pd(c)=sqrt((y-cent(1)/xyum)^2 + (x-cent(2)/xyum)^2);
end,  end
pd=pd*xyum; %change pixels to microns

%Find distance of Seg Midpoint to Cell Center
SegDist=sqrt((SegMP(:,1)-cent(1)).^2+(SegMP(:,2)-cent(2)).^2);

BCDist=sqrt((BC(:,1)-cent(1)).^2+(BC(:,2)-cent(2)).^2);
xyDend=0; xyDot=0;
for i =  0:1:max(SegDist) %run all bins
   sarea(i+1)=((sum(pd>i & pd<(i+xyBin)))*xyum^2);
   xyDend(i+1)=sum(SegLength(SegDist>i & SegDist<(i+xyBin)));
   xyDot(i+1)=sum(BCDist>i & BCDist<(i+xyBin));
end  
xyDotDend=xyDot./xyDend; %make density with distance 
xyDend=xyDend./sarea;
xyDot=xyDot./sarea;
Results.XY.Dend=xyDend;
Results.XY.Dot=xyDot;
Results.XY.DotDend=xyDotDend;
Results.XY.xyBin=xyBin;

%% run the same bins for each arbor
for a = 1: size(Results.Arbor)
    Top=Results.Arbor(a).Top;
    Bottom=Results.Arbor(a).Bottom;
    BCDistA=BCDist((BC(:,3)>Top) & (BC(:,3)<Bottom));  %extract arbors dots
    SegA=(SegMP(:,3)>Top) & (SegMP(:,3)<Bottom); %ID appropriate segs
    SegLengthA=SegLength(SegA); %extract appropriate lengths
    SegDistA=SegDist(SegA);
    
    xyDend=0; xyDot=0;
    for i =  0:1:max(SegDist) %run all bins
       sarea(i+1)=((sum(pd>i & pd<(i+xyBin)))*xyum^2);
       xyDend(i+1)=sum(SegLengthA(SegDistA>i & SegDistA<(i+xyBin)));
       xyDot(i+1)=sum(BCDistA>i & BCDistA<(i+xyBin));
    end  
    xyDotDend=xyDot./xyDend; %make density with distance 
    xyDend=xyDend./sarea;
    xyDot=xyDot./sarea;
    Results.XY.Arbor(a).Dend=xyDend;
    Results.XY.Arbor(a).Dot=xyDot;
    Results.XY.Arbor(a).DotDend=xyDotDend;
    
end
    
 

plot(xyDotDend)
pause(.01)
if isdir([TPN 'images'])==0, mkdir([TPN 'images']); end %create directory to store steps
save([TPN 'data\Results.mat'],'Results');


%% Enter Data to Structured Array
%{
load('./Dat.mat')
for i= 1:size(Dat,2)
    if strcmp(Cname,Dat(i).name); targ=i; break
    else targ=size(Dat,2)+1; end
end

%Get information
Dat(targ).type=input('What is the Cell type?  ', 's')
Dat(targ).age=input('What is the Cell age? ')
Dat(targ).notes=input('Enter notes here ==> ', 's')


Dat(targ).name=Cname;
Dat(targ).ImageInfo.xyum=xyum;
Dat(targ).ImageInfo.zum=zum;
Dat(targ).CellStats.TotalLength=TotalLength;
Dat(targ).CellStats.TotalDots=TotalDots;
Dat(targ).CellStats.AverageDensity=AverageDensity;
Dat(targ).CellStats.DotsBinDepth=DotBD;
Dat(targ).CellStats.DendBinDepth=DendBD;
Dat(targ).CellStats.DotPerDendDepth=DotPerDendDepth;
Dat(targ).CellStats.AverageingBinWidthInMicrons=bin2*zum;
for i=1:size(P,1)
    Dat(targ).Arbor(i).PeakFind=PeakFind;
    Dat(targ).Arbor(i).BoarderFind=BoarderFind;
    Dat(targ).Arbor(i).Length=ALength(i);
    Dat(targ).Arbor(i).Dots=ADots(i);
    Dat(targ).Arbor(i).DotDend=ADotDend(i);
    Dat(targ).Arbor(i).Top=left(i)*zum;
    Dat(targ).Arbor(i).Peak=P(i)*zum;
    Dat(targ).Arbor(i).Bottom=right(i)*zum;
end

save('./Dat.mat','Dat')


%% Draw Dot and Dend

'Drawing Dots and Dendrites'

Sc=(1/xyum)/2;
DD=uint8(zeros(round(max(max(AllSeg(:,1,:)))*Sc),round(max(max(AllSeg(:,2,:)))),round(max(max(AllSeg(:,3,:)))*Sc)));


%%Draw Segments
SkelRes=.1;
for i=1:size(AllSeg,1)
    Dist=sqrt((AllSeg(i,1,1)-AllSeg(i,1,2))^2 + (AllSeg(i,2,1)-AllSeg(i,2,2))^2 + (AllSeg(i,3,1)-AllSeg(i,3,2))^2); %find distance
    Length(i)=Dist;
      devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
    for d=1:devs+1
        sy=AllSeg(i,1,1)+((AllSeg(i,1,2)-AllSeg(i,1,1))/devs)*(d-1);
        sx=AllSeg(i,2,1)+((AllSeg(i,2,2)-AllSeg(i,2,1))/devs)*(d-1);
        sz=AllSeg(i,3,1)+((AllSeg(i,3,2)-AllSeg(i,3,1))/devs)*(d-1);
        DD(round(sy*Sc)+1,round(sx*Sc)+1,round(sz*Sc)+1)=1; %draw Skel
    end
end
clear Dist

%%DrawNodes
for i=1:size(AllSeg,1)
    DD(round(AllSeg(i,1,1)*Sc)+1,round(AllSeg(i,2,1)*Sc)+1,round(AllSeg(i,3,1)*Sc)+1)=2;
    DD(round(AllSeg(i,1,2)*Sc)+1,round(AllSeg(i,2,2)*Sc)+1,round(AllSeg(i,3,2)*Sc)+1)=2;
end

%%Draw Dots
for i=1:size(BC,1)
    DD(round(BC(i,1)*Sc)+1,round(BC(i,2)*Sc)+1,round(BC(i,3)*Sc)+1)=3;
end

imwriteNp(TPN,DD,'DD')


%}
'Done'