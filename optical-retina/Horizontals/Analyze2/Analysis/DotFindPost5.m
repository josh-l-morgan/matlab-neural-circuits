%function[]=anaDF()
%global DPN DFN TPN

%% Dot Finder 
%F3 has been modified to handle 8bit single tiff stacks
%large files will be broken into smaller blocks and then recombined
'start', StartTime=clock, 
tic


%% Get file names
DPN=GetMyDir

%Get directory name
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)
if isdir('./history')==0, mkdir('./history'); end %create directory to store steps
save(['./history' TPN(f3:f2-1)],'TPN') %record path in history folder

colormap gray(255) %standard grey colormap
if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps
 


%% Enter Variables

%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions
channels=3;
BlockSize=200; %aproximate size of individual processing blocks
BlockBuffer=30; %amount of block to be cut off as edge buffer

%%Dot criteria
step=2; %Sensitivity=grey value step of iterative threshold (2 was standard)
MaxDot=7^3;  %Maximum Dot Volume (in pixels)= maximum dot size for iterative threshold (6^3 was standard)
MinDot=3;  %Minimum Dot Volume (in pixels)= minimum dot size for iterative threshold   (10 was standard)
PunctaThreshold=1; %Minimum Number of Steps Passed = centroids less than puncta threshold are zeroed. 
EdgeOfPeak=.5; %Determine Dot Edge = ratio of edge brightness to peak brightness (.5 was standard)
minFilledVolume=3; %minimum number of pixels in final object contour (20 was standard)
RoundThreshold=50; %minimum roundness threshold (60 was standard)


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'

%%Read IgmTsum if it exists

%{
if exist([TPN 'data\IgmTsum.mat'])
    load([TPN 'data\IgmTsum.mat'])
elseif exist([TPN 'pics\IgmTsum'])
    TsumDir=dir([TPN 'pics\IgmTsum']);
    TsumDir=TsumDir(3:size(TsumDir));
   
    for i=1:size(TsumDir,1)
        IgmTsum(:,:,i)=imread([TPN 'pics\IgmTsum\' TsumDir(i).name]);
    end
    save([TPN 'data\IgmTsum.mat'],'IgmTsum')
end
%}

%%Find Images
d=dir(DPN); d=d(3:size(d,1)); %get number of files in directory
planes=size(d,1); %find number of planes
I(:,:,:)=imread([DPN d(1).name]); [ys,xs,zs]=size(I);    %Read first pic
Csums=sum(sum(I)); %find sums of channels in each plane
if Csums(3)>0, channels=2; else, channels=1; end
IgRaw=zeros(ys,xs,planes,'uint8');


%%Figure out channels
if size(I,3)==1, channels=1; %if only one channel
elseif size(I,3)==2, channels=1; %if only two channels
elseif sum(sum(I(:,:,3)))==0, channels=1; %if third channel is blank
else channels=2; %if third channel is not blank
end

for i=1:planes
    I(:,:,:)=imread([DPN d(i).name]);
    IgRaw(:,:,i)=uint8(I(:,:,channels)); %ColorSeperate
end
%%PreSample
%IgRaw=IgRaw(1330:1430,950:1050,:);

clear I
'Save Raw files'
save([TPN 'temp\IgRaw.mat'],'IgRaw')

%% Find Image Stats
[Rys,Rxs,Rzs]=size(IgRaw); %get sizes
if exist('IgRaw') ==0, load([TPN 'temp/IgRaw.mat']); end %load green
lowpoint=min(IgRaw(:));highpoint=max(IgRaw(:)); %find top and bottom
GstatsH=hist(IgRaw(1:min(10000000,Rys*Rxs*Rzs)),1:1:255); %find mode (mean background)
for g=1:255  %look for 95% point
if sum(GstatsH(1:g))/sum(GstatsH)>=.95, Gmode=g, break, end
end %end g = looking for 95% point


clear IbRaw IrRaw

%% Create final data output
BigCentroid=zeros(Rys,Rxs,Rzs,'uint8');
save([TPN 'temp/BigCentroid.mat'],'BigCentroid');
clear BigCentroid

BigFilled=zeros(Rys,Rxs,Rzs,'uint8');
save([TPN 'temp/BigFilled.mat'], 'BigFilled');
clear BigFilled

BigIT=zeros(Rys,Rxs,Rzs,'uint8');
save([TPN 'temp/BigIT.mat'], 'BigIT');
clear BigIT


%% SubSample
%%Find Real Block Size
'subsampling'

if Rxs>BlockSize
    NumBx=round(Rxs/BlockSize);
    Bxc=fix(Rxs/NumBx); %Block x dim closest to BlockSize for the image
else Bxc=Rxs; NumBx=1; end 
if Rys>BlockSize
    NumBy=round(Rys/BlockSize);
    Byc=fix(Rys/NumBy); %Block x dim closest to BlockSize for the image
else Byc=Rys; NumBy=1; end
if Rzs>BlockSize
    NumBz=round(Rzs/BlockSize);
    Bzc=fix(Rzs/NumBz); %Block x dim closest to BlockSize for the image
else Bzc=Rzs; NumBz=1; end

%% Run Blocks
for Bz=1:NumBz, for By=1:NumBy, for Bx=1:NumBx
PercentBlocksDone=((Bz-1)*NumBy*NumBx+Bz   +  (By-1) * NumBx + By  + Bx)/(NumBz*NumBx*NumBy)

%Find real territory
Tystart=(By-1)*Byc+1;
Txstart=(Bx-1)*Bxc+1;
Tzstart=(Bz-1)*Bzc+1;
if By<Byc, Tyend=By*Byc, else Tyend=Rys, end
if Bx<Bxc, Txend=Bx*Bxc, else Txend=Rxs, end
if Bz<Bzc, Tzend=Bz*Bzc, else Tzend=Rzs, end

%Find buffered Borders (extend to image boarders for last blocks in row and column)
ystart=Tystart-BlockBuffer; ystart=max(1,ystart);
yend=Tyend+BlockBuffer; yend=min(Rys,yend);
xstart=Txstart-BlockBuffer; xstart=max(1,xstart);
xend=Txend+BlockBuffer; xend=min(xend,Rxs);
zstart=Tzstart-BlockBuffer; zstart=max(1,zstart);
zend=Tzend+BlockBuffer; zend=min(zend,Rzs);

ViewTrans=logical(zeros(Rys,Rxs));
ViewTrans(ystart:yend,xstart:xend)=1;
image(max(ViewTrans,[],3)*1000),pause(.01)


%%Subsample Green Channel
if exist('IgRaw') ==0, load([TPN 'temp/IgRaw.mat']); end %load green
Ig=single(IgRaw(ystart:yend,xstart:xend,zstart:zend)); clear IgRaw %subsample
%image(max(Ig,[],3)*255/max(Ig(:))) %Image green

%Find New sizes
[ys,xs,zs]=size(Ig);

%% Median filter 
'median filtering'
Igm=Ig*0;
for i=1:zs        
    Igm(:,:,i)=medfilt2(Ig(:,:,i),[3,3]); 
end
    
%%Image median filtered data 
%image(max(Igm,[],3)*255/max(Igm(:))) %Image green
pause(.01)

'median filter done'

%% FIND DOTS Green Channel%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'iterative threshold',
centroid=Igm*0; %set up matrix to map centroids
IgmTsum=Igm*0; %set up matrix to sum passed thresholds
highpoint=max(Igm(:));
steps=fix((highpoint-Gmode)/step)+2;

 for i = fix(highpoint)+2:-step:Gmode %run thresholds through all relevant intensities

       %IterativeThresholdPercentDone=uint8((i-Gmode)/(highpoint+2-Gmode)*100)

       %%run threshold
       clear IT
       IT(:,:,:)=Igm>i;
       [Igl,labels]=bwlabeln(IT(:,:,:),6); % label each area to check
       if labels<65536, Igl=uint16(Igl); end  %shrink bitdepth if possible
       IglH=hist(Igl(Igl>0),1:1:labels+2);  %count pix for each labeled object

       %identify pixels for each puncta  
       for p=1:labels  %run all lables
           if IglH(p)< MaxDot & IglH(p)>MinDot  %% Morphology Filter ! Puncta size criteria !
            %% Find Peaks to make centroid
                if sum(sum(sum(logical(centroid(Igl==p))))) ==0 
                    [tsy tsx tsz]=find3(Igl==p);
                    centroid(round(mean(tsy)),round(mean(tsx)),round(mean(tsz)))=1;
                end % if region doesnt alrea
           
           else %end if the right size
               Igl(Igl==p)=0; %clear wrong size object
         end% end check size        
           
       end  %end do each label

       %%Add all passing labeled objects to IgTsum
       IgmTsum(Igl>0)=IgmTsum(Igl>0)+1;  

end  %end iterative threshold

%image(sum(centroid,3)*100),pause(.01)

IgmTsum=uint8(IgmTsum);
clear IT Igl


%% Find dot countour %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
'finding dot contour',pause(.01)

%%find positions of dots, give dots an ID 
centroid(IgmTsum==0)=0; %safegard to eliminate centroids outside of threshold area (shouldnt be necessary)
[Py,Px,Pz]=find3(centroid>0);
numdots=size(Py,1);
filled=uint16(zeros(ys,xs,zs));
Dotmap=uint8(zeros(ys,xs,zs));  %make 3D map where positions are labled with dot ID
for i=1:numdots;   Dotmap(Py(i),Px(i),Pz(i))=i; end

for n=1:numdots%run dots
    %PercentDoneFindingContour=100*n/numdots
    
    %% Characterize puncta core (find IgmTsum val at centroid and peak of nearby pixels
    peak=Igm(Py(n),Px(n),Pz(n)); %value at centroid
    countTS=0; %make counter for summed thresholds
    %grab all values with 5pix(cardinal) from centroid
    for yc=Py(n)-2:Py(n)+2, for xc=Px(n)-2:Px(n)+2, for zc=Pz(n)-2:Pz(n)+2;
               if yc>0 & yc<=ys,if xc>0 & xc<=xs, if zc>0 & zc<=zs
               countTS=countTS+1;
               TS(countTS)=IgmTsum(yc,xc,zc); %record values
               end,end,end %end walls
    end,end,end %end search yxz
    TS=sort(TS,'descend'); %decending list of values near puncta
    peakTS=median(TS(1:4));  %take median of top four values withing 2pix of centriod
    Cthresh=min((peakTS*EdgeOfPeak),IgmTsum(Py(n),Px(n),Pz(n))); %insure centroid is included

    %% Identify region as countour based on peak
    ContourL=single(IgmTsum>=max(Cthresh,1)); %Threshold (minimum threshold is 1 pass)
    [ContourL]=single(bwlabeln(ContourL,6)); %label countour based on connectivty
    Ctarget=ContourL(Py(n),Px(n),Pz(n)); %Pic object overlapping center
    
    if sum(sum(sum(Dotmap(ContourL==Ctarget)>0)))>1 % if contour includes more than 1 centroid
        %% sort pixes in area according to who is closest
        [ay,ax,az]=find3(ContourL==Ctarget); numa=size(ay,1); %find each pix in countour
        ContourL=single(Dotmap).*single(ContourL==Ctarget); %restrict map to area
        [pdy,pdx,pdz]=find3(ContourL>0); numpd=size(pdy,1); % find centroids in area
        AtoPD=1;
        for a=1:numa %run all pixels
            for pd=1:numpd %end check all centroids
                AtoPD(pd)=sqrt((ay(a)-pdy(pd))^2+(ax(a)-pdx(pd))^2+(az(a)-pdz(pd))^2); %find distances
            end %end chack all centroids
            closest=find(AtoPD==min(AtoPD));
            if size(closest,1)>1  %if equal distance
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!INSERT VALUE CHECKER!!!!!!!!!!!1
            end %end if equal distance
            if Dotmap(pdy(closest),pdx(closest),pdz(closest))==n %if its the dot currently running
                filled(ay(a),ax(a),az(a))=Dotmap(pdy(closest),pdx(closest),pdz(closest)); %enter ID
            end %if current dot wins
        end %end run all pixels (a)
    else
        filled(ContourL==Ctarget)=n;
          
    end % if there is more than one centroid within countour
        
end %run contour for each dot

%image(max(filled,[],3)*10000),pause(.01) %image filled contours


%% Size and roundness filtering of final contour
%{
'running contour size and roundness filter'
%buffer filled
filled2=single(zeros(ys+10,xs+10,zs+10));
filled2(6:ys+5,6:xs+5,6:zs+5)=filled;
filled=filled2; clear filled2
sym=filled; %matrix for recording roundness
%make spheres (distance map)
look=5;cent=look+1;
sphere=single(zeros(cent+look,cent+look,cent+look));
for yr=1:look*2+1,for xr=1:look*2+1,for zr=1:look*2+1
    sphere(yr,xr,zr)=sqrt((yr-cent)^2+(xr-cent)^2+(zr-cent)^2);
end,end,end %end make sphere

for n=1:max(filled(:))
    %PercentDoneCheckingRoundness=n/max(filled(:))*100

    %% extra final size filter
    [filledpix]=find3(filled==n);
    size(filledpix);
    if size(filledpix,1)<minFilledVolume,filled(filled==n)=0;
    else %check roundness
            
    %% check roundness
    [fy,fx,fz]=find3(filled==n); %find area
    fym=round(mean(fy));fxm=round(mean(fx));fzm=round(mean(fz)); %find center
    checkf=single(filled(fym-look:fym+look,fxm-look:fxm+look,fzm-look:fzm+look)>0); %sample filled
    %image(max(checkf,[],3)*300),pause(.01)
    for r=1:5
        heck=(checkf.*(sphere>r-1 & sphere<=r));
        %image(sum(heck,3)*30),pause(.1)
        hell(r)=sum(heck(:));
        crap(r)=sum(sum(sum(sphere>r-1 & sphere<=r)));
        shit(r)=abs(.5-hell(r)/crap(r));
    end
    roundness=mean(shit)*200; %imageable roundness rateing (1:100)
    sym(filled==n)=roundness; %assign roundness value
    
    end %if big enough check roundness
end %size and symetry countour filter
filled=filled(6:ys+5,6:xs+5,6:zs+5); %unbuffer filled
sym=sym(6:ys+5,6:xs+5,6:zs+5);

%apply sym filter to filled
filled=filled.*(sym>RoundThreshold);

%image(sum((filled(:,:,:)>0),3)*300),pause(.01)
clear ContourT ContourL ContourD Igm ay ax az
clear  round sym Dotmap Igm 
%}

%% Recombine Blocks
%%Unbuffer
centro=uint8(centroid(Tystart-ystart+1:ys-(yend-Tyend), Txstart-xstart+1:xs-(xend-Txend), Tzstart-zstart+1:zs-(zend-Tzend)));
fille=uint8(filled(Tystart-ystart+1:ys-(yend-Tyend), Txstart-xstart+1:xs-(xend-Txend), Tzstart-zstart+1:zs-(zend-Tzend)));
IgmTsu=uint8(IgmTsum(Tystart-ystart+1:ys-(yend-Tyend), Txstart-xstart+1:xs-(xend-Txend), Tzstart-zstart+1:zs-(zend-Tzend)));
clear centroid filled 

%%Enter into Big figure
load([TPN 'temp/BigCentroid']); 
BigCentroid(Tystart:Tyend,Txstart:Txend,Tzstart:Tzend)=centro;
save([TPN 'temp/BigCentroid'],'BigCentroid')
clear BigCentroid centro

load([TPN 'temp/BigFilled']);
BigFilled(Tystart:Tyend,Txstart:Txend,Tzstart:Tzend)=fille;
save([TPN 'temp/BigFilled'],'BigFilled')
clear BigFilled fille

load([TPN 'temp/BigIT']);
BigIT(Tystart:Tyend,Txstart:Txend,Tzstart:Tzend)=IgmTsu;
save([TPN 'temp/BigIT'],'BigIT')
clear BigIT IgmTsu

%% Go to next Block
clear filled centroid IgmTsum
end,end,end %%End Bz, By, Bz Block translation


%% Recombine Blocks
load([TPN 'temp/BigCentroid']);
save([TPN 'data/BigCentroid'],'BigCentroid'); 
imwriteNp(TPN,BigCentroid,'BigCentroid')
clear BigCentroid


load([TPN 'temp/BigFilled']);
save([TPN 'data/BigFilled'], 'BigFilled');
imwriteNp(TPN,BigFilled,'BigFilled')
subplot(1,1,1)
image(max(BigFilled,[],3)*1000)

load([TPN 'temp/BigIT']);
save([TPN 'data/BigIT'], 'BigIT');
imwriteNp(TPN,BigIT,'BigIT')
subplot(1,1,1)
image(max(BigIT,[],3)*1000)



%% Vectorize BigFilled to Make Dots
load([TPN 'data\BigFilled.mat'])
Puncta=zeros(1,4);
PidCount=0; %start counter for puncta ID
for i=1:max(BigFilled(:)) % run all puncta
    clear y x z p Dis
    [y x z]=find3(BigFilled==i);
    p=[y x z]; 
    p(:,4)=0;
    Checked=ones(size(p,1),1);
    %% ID puncta
    while sum(Checked) %while there are more pixels to check
        Targ=find(Checked,1);  %find first unchecked pixel
        PidCount=PidCount+1;  %increase puncta ID counter
        p(Targ,4)=PidCount;
        while sum(p(:,4)==PidCount & Checked) %find all pixels for that ID
                Targ=find(p(:,4)==PidCount & Checked,1);
                for d=1:3;    Dis(:,d)=abs(p(:,d)-p(Targ,d))>1;   end  %find distances
                Con=~sum(Dis,2);  %finds all connected 
                p(Con,4)=PidCount;  %convert all new connected to current ID
                Checked(Targ)=0;  % Clear checked pixel       
        end %while more from group not checked
    end %while Checked (when there are more to check)
    Puncta=cat(1,Puncta,p);
    PercentDoneVectorizing=(single(i)/single(max(BigFilled(:))))*100
end %end i , run all puncta
Puncta=Puncta(2:size(Puncta,1),:);

%%Create Dots
Dots.Num=max(Puncta(:,4));
[Fys Fxs Fzs] = size(BigFilled);
clear BigFilled
siz=[Fys Fxs Fzs];
Dots.ImSize=siz;
for i = 1: Dots.Num
   Dots.Vox(i).Pos=Puncta(Puncta(:,4)==i,1:3); 
   Dots.Vox(i).Ind=sub2ind(siz,Dots.Vox(i).Pos(:,1),Dots.Vox(i).Pos(:,2),Dots.Vox(i).Pos(:,3));
   Dots.Vol(i)=size(Dots.Vox(i).Pos,1);   
   Dots.Pos(i,:)=mean(Dots.Vox(i).Pos,1);   
end

load([TPN 'data\BigIT.mat'])
for i=1:Dots.Num
    Dots.ITMax(i)=max(BigIT(Dots.Vox(i).Ind));
    Dots.ItSum(i)=sum(BigIT(Dots.Vox(i).Ind));
end
clear BigIT

load([TPN 'temp\IgRaw.mat'])
for i = 1: Dots.Num
    Dots.Vox(i).RawBright=IgRaw(Dots.Vox(i).Ind)
    Dots.MeanBright(i)=sum(Dots.Vox(i).RawBright)/Dots.Vol(i)
end
clear IgRaw

save([TPN 'Dots.mat'],'Dots')

%% Finish

TotalHours=toc/60/60
[TPN(size(TPN,2)-6:size(TPN,2)-1)]
DotFindAt=uint16(clock)
save([TPN 'data/DotFindAt.mat'],'DotFindAt')

%clear all
'Done DotFind'
profile off



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Notes:

DO:
Do percent over brightness predicted from dendrite GFP TdTomato ratio.
Concider finding length by stringing maxes together. 
factor in aspect ratio
jitter on dendrite instead of just IPL
should somehow account for channel bleed through


DID:
Identified dots by Iterative threshold
Identified dot volume

%}

