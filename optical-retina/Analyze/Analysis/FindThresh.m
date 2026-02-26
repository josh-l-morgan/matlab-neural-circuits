function[]=anaFT(DPN)

%% Dend Finder 
%F3 has been modified to handle 8bit single tiff stacks
%large files will be broken into smaller blocks and then recombined
%Strategy: 1)erode into wire frame.  
%          2)measure length segments (or convolve with lenght look up table)
%          3)for each segment find radius based on smallest xy diameter

'start'
clear all
tic
colormap gray(255) %standard grey colormap

%Get file names


[DFN,DPN]=uigetfile('*.tif','DialogTitle','Choose first Image of Data Stack');


%Get directory name
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)


if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps
if isdir('./historyT')==0, mkdir('./historyT'); end %create directory to store steps

save(['./historyT' TPN(f3:f2-1)],'TPN') %record path in history folder

%%Enter Variables

%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions
live=1; %if live tissue then 1, 2 prevents threshold climbing in fixed tissue
save([TPN 'data\ImageName.mat'],'DPN')


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'reading image'

if exist([TPN 'temp/Imax.mat']),
    load([TPN 'temp/Imax.mat'])
else
    ns=size(DFN,2); %find size of name

    %%Get number of planes and spacer zeros
    d=dir(DPN); %get number of files in directory
    planes=size(d,1)-2; %find number of planes
    pdiddy=fix(log10(planes))+1;

    %%Figure out channels
    nameD=[DPN DFN(1:ns-4-pdiddy) num3str(1,pdiddy) DFN(ns-3:ns)];
    Ic(:,:,:)=imread(nameD); %read
    if size(Ic,3)==1, channels=1; %if only one channel
    elseif size(Ic,3)==2, channels=2; %if only two channels
    elseif sum(sum(Ic(:,:,3)))==0, channels=2; %if third channel is blank
    else channels=3; %if third channel is not blank
    end
    Imax=Ic(:,:,channels);

    clear I Ic
    for i=1:planes
        nameD=[DPN DFN(1:ns-4-pdiddy) num3str(i,pdiddy) DFN(ns-3:ns)];
        Ic(:,:,:)=imread(nameD); %read
        I=Ic(:,:,channels);%ColorSeperate
        I=medfilt2(I,[3,3]); %median filter
        Imax=max(Imax,I); %Find max
        PercentRead=i/planes*100
    end

    %%PreSample
    %Irm=Irm(800:900,800:900,:);

    clear I
    'Save Imax'
    save([TPN 'temp/Imax.mat'],'Imax')
    imwriteNp(TPN,Imax,'Imax')
end %end read if max not available


[ys,xs,zs]=size(Imax); %get sizes

%%Image median filtered data 
subplot(1,1,1), image(Imax*255/max(Imax(:))) %Image green
pause(.01)

'Done reading and median filtering'


%% FIND DENDRITES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'finding threshold'

%% Find statistics
stats=sort(Imax(1:min(10000000,ys*xs*zs))); %make one dimensional list of image values
dev=std(stats);
statsH=hist(stats,1:1:255);
low=find(statsH==max(statsH),1); %take the median value of the dimest half of image

for g=1:255  %look for 95% point
if sum(statsH(1:g))/sum(statsH)>=.90, Thresh=g, break, end
end %end g = looking for 95% point

myThresh=mean(stats)+dev*3  %Arbitrary but not unreasonable. 
save([TPN 'data/myThresh.mat'],'myThresh')
%clear stats
profile off

%% User Threshold

load([TPN 'temp/Imax.mat'])
image(Imax), pause(.01)


thresh=0
while thresh ~= 255
thresh=input('try threshold (255 to exit) ==>')
pause(.1),image((Imax>thresh)*300),pause(.1)
end

Threshold=input('what is the final threshold?  ');
save([TPN 'data/Threshold.mat'], 'Threshold');
pause(.1),image((Imax>thresh)*300),pause(.1)

clear all
