clear all
TPN = GetMyDir;
TPNm = [TPN 'I\'];
dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));

if ~exist([TPN 'data']), mkdir([TPN 'data']),end
if ~exist([TPN 'temp']), mkdir([TPN 'temp']),end


%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions

%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'reading image'


if size(dTPNm,1)==1
    Ir=tiffread2([TPNm dTPNm(1).name]);
    D=zeros(size(Ir(1).data,1),size(Ir(1).data,2),size(Ir,2));
    for i = 1: size(Ir,2)
        D(:,:,i)=Ir(i).data;
    end
else
    
    clear I Ic
    D(:,:)=imread([TPNm dTPNm(1).name]); %read
    D(1,1,size(dTPNm,1))=0;
    for i=1:size(dTPNm,1)
        D(:,:,i)=imread([TPNm dTPNm(i).name]);
        PercentRead=i/size(d,1)*100
    end
end


%Find New sizes
[ys,xs,zs]=size(D);
sumD=sum(D,3);
colormap gray(256)
image(sumD*50/(mean(sumD(:))))


%% Initial Threshold

Gmode=mode(D(:));
Gstd=std(D(:));

T=D>(Gmode+Gstd);
image(sum(T,3)*500)

%% Dilate Erode
C = GetCon(T,26);
T(C==0)=0;

C = GetCon(T,26);
T=T|(C>0);
    image(sum(T,3)*50),pause(.1)
C=GetCon(~T,26);
T=T&~(C>0);
    image(sum(T,3)*50),pause(.1)


%% Trim singles


image(sum(T,3)*60),pause(.01)

P=bwperim(T,6);
image(sum(P,3)*30),pause(.01)

for i = 1 : 10
    P=bwperim(T,6);
    TmP=T & ~ P;
    Touch=GetCon(TmP);
    T=T & ~(P & (Touch>0));
    image(sum(T,3)*60),pause
end

for i = 1: size(T,3)
    image(T(:,:,i)*1000),pause
end




















