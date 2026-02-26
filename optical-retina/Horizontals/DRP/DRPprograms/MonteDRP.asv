%% Monte DRP

%% Give  both 3D and 2D versions of DRP
%Helps to filter the crap out of image prior to thresholding ( so long as
%objects dont touch

%%Resolutions
% 60x z2 1024 = .103 FV1000
% 25x z1 1024 = .4944 ? not tested
clear all
colormap gray(256)
reps = 1000;

TPN = GetMyDir;
load([TPN 'Dat.mat'])
ImageInfo=Dat.Image.ImageInfo;
Rebin=Dat.Rebin;

minObSize=ImageInfo{1};
Dim = ImageInfo{2};
Objective = ImageInfo{3};
Zom = ImageInfo{4};
Resolution = ImageInfo{5};
zum = ImageInfo{6};
Bin = ImageInfo{7};
Look = ImageInfo{8};
StandardRes = ImageInfo{9};

xyum=(StandardRes/Zom)*(60/Objective)*(2048/Resolution);


%% Get image
'getting image'
TPNi=[TPN]
Idir=dir(TPNi);
In={};
for i = 1: length(Idir)
    name=Idir(i).name;
    LN=length(name);
    if LN>=3
        if length(name)>=4
            if sum(name(LN-3:LN)=='.tif')==4;
                In{length(In)+1}=name;
            end
        end
    end
end

I = imread([TPNi '\' In{1}]);
I = max(I,[],3);
siz = [size(I,1) size(I,2) length(In)];

if siz(3)>1
    I = zeros(siz,'uint8');
    for i = 1:siz(3)
        I(:,:,i)=imread([TPNi '\' In{i}]);
    end
end
Imax=max(I,[],3);
image(Imax'*1000);

I = I > median(I(:)); %remove any background
I = uint16(I);


%% Get Pos
[I numOb] =bwlabeln(I,26);
Pos=[];
for i = 1: numOb
    ObInd=find(I==i);
    if length(ObInd)>=minObSize  %if big enough
        [y x z] = ind2sub(siz,ObInd);
        Pos(size(Pos,1)+1,:)=mean([y x z],1);
    end %% if big enough
    ['imaging ' num2str(i) ' of ' num2str(numOb)]
end

realPos=Pos;
Ranges=[max(Pos,[],1);min(Pos,[],1)];


%% Create distance matrix
maxDist=Look;%max(Dists(:));
DUsed=(0:Look)';

LookPix=(1./[xyum xyum]) * Look;
Mid=fix(LookPix)+1;
DmatSize=Mid*2+1;
Dmat=zeros(DmatSize);
[dmy dmx]=ind2sub(DmatSize(1:2),find(Dmat(:,:)==0));
Dmat(:)=sqrt((dmy-Mid(1)).^2 + (dmx-Mid(2)).^2 );
Dmat=Dmat*xyum;



%% REpeat
for r = 1:reps
    ['running ' num2str(r) ' of ' num2str(reps)],pause(.01)

%% Randomize
rPos=rand(size(Pos));
for i = 1:size(Ranges,2)
    rPos(:,i)=rPos(:,i)*(Ranges(1,i)-Ranges(2,i))+Ranges(2,i);
end
Pos=round(rPos);

%% Scale Positions
Pos(:,1:2)=Pos(:,1:2)*xyum;
Pos(:,3)=Pos(:,3)*zum;

%% Find Dists to centers
clear Centers I
Dnum=size(Pos,1);
Dists=zeros(Dnum);
for d = 1:Dnum
    Dists(d,:)=sort(dist(Pos(:,1:2), Pos(d,1:2)));  %use rounded positions
end


%% Nearest Neighbor
mDists=mean(Dists,1);
Near=Dists(:,2);
hNear=hist(Near,2.5:1:max(Near(:)));

%% Density Recovery Profile
' Running DRP'

Cnum2=zeros(Dnum,size(DUsed,1));  %hist of cell bodies relative to each cellbody
Anum2=Cnum2; Anum3=Anum2; Cnum3=Cnum2;%hist of area pixels relative to each cell body
DistsA=Dnum;

for c = 1:Dnum
    %['running ' num2str(c) ' of ' num2str(Dnum)],pause(.01)
    cPos=Pos(c,:); %% find center
    cPos=round(cPos./[xyum xyum zum]);
    dPos=siz-cPos; %% find difference from size (direction 2)
    pPos=[cPos;dPos];
    Vmat=Dmat*0;

    startP = max(1,Mid-cPos(1:2)+1);
    stopP = min([DmatSize ;Mid+dPos(1:2)],[],1);
    Vmat(startP(1):stopP(1),startP(2):stopP(2))=1;

    %%find 3D cells and Volumes

    for d = 2:length(DUsed)
        
        Anum2(c,d)=sum(sum(sum(Vmat(Dmat<DUsed(d) & Dmat>DUsed(d-1)))));
        Cnum2(c,d)=sum((Dists(c,:)<=DUsed(d)) & (Dists(c,:)>DUsed(d-1)));
%         Smat=Vmat*100;
%         Smat(Dmat<DUsed(d) & Dmat>DUsed(d-1))= Smat(Dmat<DUsed(d) & Dmat>DUsed(d-1))+100;
%         image(Smat),DUsed(d),pause
    end

end

%% Convert to area/volume
pix=(xyum*xyum)/(1000)^2;
vox=pix * zum/1000;
Anum2=Anum2*pix;

Anum=Anum2;
Cnum=Cnum2;


Csum=mean(Cnum,1)';
Asum=mean(Anum,1)';

CsumA(r,:)=Csum;
AsumA(r,:)=Asum;
CoA=Csum./Asum;
CoAA(r,:) =CoA;

%% Reshape
Rebin=2;
Bins=(1:Rebin:Look-Rebin);
clear Cb Ab 
numBin=length(Bins);
for d = 1 : numBin
    %Bins(d):Bins(d)+Rebin-1
    Cb(:,d)=sum(Cnum(:,Bins(d):Bins(d)+Rebin-1),2);
    Ab(:,d)=sum(Anum(:,Bins(d):Bins(d)+Rebin-1),2);
end

CA=sum(Cb,1)./sum(Ab,1);
%bar(CA),pause(.01)

CAs(r,:)=CA;

end  %% Run all



%save([TPN 'CAs_b2n1000.mat'],'CAs')
%save([TPN 'CAs_b2n1000_AllData.mat'])
%load([TPN 'CAs_b2n1000.mat'])
% 


%% Plot with errorbars
N=size(CAs,1); 
clear E
for b = 1:numBin
    E(b)=std(CAs(:,b)/sqrt(N));    
end
bar(Bins+Rebin/2,mean(CAs,1))
hold on
errorbar(Bins+Rebin/2,mean(CAs,1),E,'+')
hold off

%% Show Dist
% Showb=zeros(fix(max(CAs(:)))+1,size(CAs,2));
% for i = 1:size(CAs,1)
%    for x = 1: size(CAs,2)
%       Showb(1:fix(CAs(i,x)),x)=Showb(1:fix(CAs(i,x)),x)+1; 
%    end    
% end
% image(Showb*255/size(CAs,1))



%% Reshape
Rebin=2;
Bins=(1:Rebin:Look-Rebin);
numBin=length(Bins);
clear CAsR
for r = 1:reps
clear Cb Ab 
Csum=CsumA(r,:);
Asum=AsumA(r,:);
for d = 1 : numBin
    %Bins(d):Bins(d)+Rebin-1
    Cb(:,d)=sum(Csum(Bins(d):Bins(d)+Rebin-1),2);
    Ab(:,d)=sum(Asum(Bins(d):Bins(d)+Rebin-1),2);
end

CA=sum(Cb,1)./sum(Ab,1);
%bar(CA),pause(.01)

CAsR(r,:)=CA;

end %run all
bar(Bins,mean(CAsR),'r')
hold on

clear Cb Ab 
Csum=mean(Dat.Raw.Cnum2)
Asum=mean(Dat.Raw.Anum2);
for d = 1 : numBin
    %Bins(d):Bins(d)+Rebin-1
    Cb(:,d)=sum(Csum(Bins(d):Bins(d)+Rebin-1),2);
    Ab(:,d)=sum(Asum(Bins(d):Bins(d)+Rebin-1),2);
end

CA=sum(Cb,1)./sum(Ab,1);
bar(Bins,CA),pause(.01)
hold off

%% Test significance of CA

Num=size(CAs,1);
for i = 1: size(CA,2)
    P(i)=sum(CAs(:,i)<CA(i))/Num;
end
%bar(Bins,P)
[Bins ; P]




