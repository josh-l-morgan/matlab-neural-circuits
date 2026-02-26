function[] = anaDD(TPN)
% Dot Dend Does something useful with puncta and dendritic data from FD
%%and Skel please

'Starting analysis of Dots and Dendrites'


%% Get Path to Big Centroid
%Get directory name
DPNd=[TPN 'data\'];
Cname=TPN;
%choose source data based on availability
if exist([DPNd 'DotStats.mat'])
    BCFN='DotStats.mat'
else
    BCFN='BigCentroid.mat'
end


save([DPNd 'SourceData-' BCFN(1:size(BCFN,2)-4)]) %record source of Data


%% Vectorize Dot Data
'Vectorizing Dot Data'
xyum=.103;
zum=.3;


if strcmp(BCFN,'BigCentroid.mat')

    clear BCy BCx BCz BC
    if exist('BigCentroid','var')==0, load([DPNd BCFN]); end %load Big Centroid
    [BCy BCx BCz]=find3(BigCentroid>1);
    clear BigCentroid

    BCy=BCy*xyum; BCx=BCx*xyum; BCz=BCz*zum;
    BC(:,1)=BCy; BC(:,2)=BCx; BC(:,3)=BCz;
    
elseif strcmp(BCFN,'DotStats.mat') 
    
    load([DPNd BCFN])
    BC=DotStats(:,:,3);
    BC=BC(DotStats(:,3,1)>.5,:) %eliminate Dots with low DF0f
   
end %which to vectorize
   
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

%% Bin with Depth
'Binning Data with Depth'

%%Depth
bin=zum;%zum;
step=zum;
bin2=7; %size to average
filt=ones(bin2,1);

maxDepth=fix(max(max(SegMP(:,3)), max(BC(:,3))))+1;
minDepth=max(fix(min(min(SegMP(:,3)), min(BC(:,3)))),1);
minDepth = minDepth-mod(minDepth,zum); %Set min depth to factor of zum
DotBinDepth=zeros(round((maxDepth-minDepth)/step+1),1);
DendBinDepth=zeros(round((maxDepth-minDepth)/step+1),1);

b=0;
for d=minDepth:step:maxDepth
    b=b+1;
    DotBinDepth(b)=sum(BC(:,3)>=(d-bin/2) & BC(:,3) <= (d+bin/2));
    DendBinDepth(b)=sum(SegLength(SegMP(:,3)>=(d-bin/2) & SegMP(:,3) <= (d+bin/2)));
end %end d, run all depth bins

%%filter to smooth
DotBD=filter(filt,bin2,DotBinDepth)
DendBD=filter(filt,bin2,DendBinDepth);


DotPerDendDepth=DotBD./DendBD;

subplot(3,5,1)
plot(DendBD,'r')
subplot(3,5,6)
plot(DotBD,'g')
subplot(3,5,11)
plot(DotPerDendDepth)
pause(.01)
%}

pause(2)
%% Identify Strata

%%manual Identify Strat
PeakFind='Manual'
'Select peaks and click return'
[x y]=ginput
P=round(x);
peaks=DendBinDepth*0;
if size(x,1)>1
    'Identify Valley '
    [x y]=ginput;
    Valley=round(x);
else, Valley=0; %stray value for Valley
end



BoarderFind='full width half maxima';
DendBDmin=imregionalmin(DendBD);
%%find peak edges using full width half maxima and reginal minimus
for i =1:size(P,1)
    top=DendBD(P(i));
    left(i)=1; right(i)=size(DendBD,1);
    for l=1:P(i)-1
        if (DendBD(P(i)-l)< top/2) | sum((P(i)-l)==Valley)
            left(i)=P(i)-l; break;   end %search for left side
    end %end searching for left side
    for r=1:size(DendBD,1)-P(i)
        if (DendBD(P(i)+r)< top/2) | sum((P(i)+r)==Valley)
            right(i)=P(i)+r; break;   end %search for left side
    end %end searching for left side
end %find peak width


%}

%%find arbor stats
APLength=DendBD*0;
APDots=DendBD*0;
APDotDend=DendBD*0;
for i=1:size(P,1)
    ALength(i)=sum(DendBinDepth(left(i):right(i)));
    APLength(left(i):right(i))=ALength(i);
    ADots(i)=sum(DotBinDepth(left(i):right(i)));
    APDots(left(i):right(i))=ADots(i);
    ADotDend(i)=ADots(i)/ALength(i);
    APDotDend(left(i):right(i))=ADotDend(i);
end %run all arbors


%%Look at Valley

if size(P,1)>1 %if bi
    vdm=DendBD*0;
    vr=2; %valley width = vr *2 +1
    for i=min(right):max(left) %run valley
        vdm(i)=mean(DotPerDendDepth(i-vr:i+vr));
    end
    vcent=find(vdm==min(vdm(min(right):max(left))));
    
    VPLength=DendBD*0;
    VPDots=DendBD*0;
    VPDotDend=DendBD*0;

    %%Find edges of valley
    vleft=vcent-vr;
    vright=vcent+vr;
    %%Find stats
    VLength=sum(DendBinDepth(vleft:vright));
    VPLength(vleft:vright)=VLength;
    VDots=sum(DotBinDepth(vleft:vright));
    VPDots(vleft:vright)=VDots;
    VDotDend=VDots/VLength;
    VPDotDend(vleft:vright)=VDotDend;
end %if bi
    


%%Draw stats
subplot(3,5,2)
plot(APLength,'r')
subplot(3,5,7)
plot(APDots,'g')
subplot(3,5,12)
plot(APDotDend,'b')
title(Cname)
if exist('VPDotDend','var')
    hold on
    plot(VPDotDend,'b')
    hold off
end
%}

if isdir([TPN 'images'])==0, mkdir([TPN 'images']); end %create directory to store steps
saveas(gcf,[TPN 'images/DD'],'ai') %save figure as illustrator file in images
save([DPNd 'DepthFigBin.mat'], 'bin2') %save size of figure bin 

pause(2)


%% Enter Data to Structured Array
    if exist([TPN 'data\Results.mat'])
        load([TPN 'data\Results.mat'])
    end
    %Get information
    %Results.type=input('What is the Cell type?  ', 's')
    %Results.age=input('What is the Cell age? ')
    %Results.notes=input('Enter notes here ==> ', 's')
    

    Results.location=Cname;
    Back=find(Cname=='\');
    Results.name=Cname(Back(size(Back,2)-1)+1:Back(size(Back,2))-1)
    Results.ImageInfo.xyum=xyum;
    Results.ImageInfo.zum=zum;
    Results.CellStats.TotalLength=TotalLength;
    Results.CellStats.TotalDots=TotalDots;
    Results.CellStats.AverageDensity=AverageDensity;
    Results.Depth.DotsBinDepth=DotBD;
    Results.Depth.DendBinDepth=DendBD;
    Results.Depth.DotPerDendDepth=DotPerDendDepth;
    Results.Depth.AverageingBinWidthInMicrons=bin2*zum;
    for i=1:size(P,1)
        Results.Arbor(i).PeakFind=PeakFind;
        Results.Arbor(i).BoarderFind=BoarderFind;
        Results.Arbor(i).Length=ALength(i);
        Results.Arbor(i).Dots=ADots(i);
        Results.Arbor(i).DotDend=ADotDend(i);
        Results.Arbor(i).Top=left(i)*zum;
        Results.Arbor(i).Peak=P(i)*zum;
        Results.Arbor(i).Bottom=right(i)*zum;
    end
    
    %save Results
    save([TPN 'data\Results.mat'],'Results')


%{
if exist('./Dat.mat')
    load('./Dat.mat')
    for i= 1:size(Dat,2)
        if strcmp(Cname,Dat(i).name); targ=i; break
        else targ=size(Dat,2)+1; end
    end
    Dat(targ)=Results;

    save('./Dat.mat','Dat')

end %if Dat exists
%}
    
    
    
%% Draw Dot and Dend
%{
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

save([TPN 'pics\DD.mat'])

%}

%% Finish

[TPN(size(TPN,2)-6:size(TPN,2)-1)]
DotDendAt=uint16(clock)
save([TPN 'data/DotDendAt.mat'],'DotDendAt')

'Done DotDend'
