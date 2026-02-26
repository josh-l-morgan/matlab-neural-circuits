function[]=axConvex(TPN,xyum)

%  TPN = GetMyDir;
% DPN=[TPN 'I\'];
TPNm = [TPN 'mask\'];
dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));

if ~exist([TPN 'data']), mkdir([TPN 'data']),end
if ~exist([TPN 'temp']), mkdir([TPN 'temp']),end

colormap gray(256)
Dil=5; %Dilation diameter
Pt=1; % time to pause at each image


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'reading image'


if size(dTPNm,1)==1
    Ir=tiffread2([TPNm dTPNm(1).name]);
    D=zeros(size(Ir(1).data,1),size(Ir(1).data,2),size(Ir,2));
    for i = 1: size(Ir,2)
        D(:,:,i)=Ir(i).data;
    end
    clear Ir
else

    clear I Ic
    D(:,:)=imread([TPNm dTPNm(1).name]); %read
    D(1,1,size(dTPNm,1))=0;
    c=0;
    for i=1:size(dTPNm,1)
        nam=dTPNm(i).name;
        naml=length(nam);
        if nam(naml-3:naml)=='.tif'
            c=c+1;
            D(:,:,c)=imread([TPNm nam]);
        end
        PercentRead=i/size(dTPNm,1)*100
    end

    D=D==1;

end



%% Area 
[ys,xs,zs]=size(D);
NT = D;
NTmax=max(NT,[],3);
se = strel('disk',Dil);
NTc=imclose(NTmax,se);
% image(NTc*1000),pause(.5)
NTp=bwperim(NTc);
image(NTp*1000),pause(.5)

%% Connect
NTcon=NTp*0;
NTpu=NTp;
res=1;
while sum(NTpu(:))>0
    [y x]= find(NTpu);
    Dists=axDist([y x],[y(1) x(1)]);
    ydif=y-y(1); xdif=x-x(1);
    for e = 1:length(x)
        for rep = 1:Dists(e)/res
            sy=y(1)+(ydif(e)/(Dists(e)/res)) * rep;
            sx=x(1)+(xdif(e)/(Dists(e)/res)) * rep;
            sy=round(sy); sx=round(sx);
            NTcon(sy,sx)=1;
            Off=abs(sy-y(e))+abs(sx-x(e));
            if Off
                NTpu(sy,sx)=0;
            end
        end
    end
    image(NTcon*50+NTp*50+NTpu*100),pause(.01)
    NTpu(y(1),x(1))=0;
end

%% Fill holes
NTfill=imclose(NTcon,se);
NThole=bwlabel(~NTfill);
Bord=NThole*0;
Bord(1:size(Bord,1),1)=1;
Bord(1:size(Bord,1),size(Bord,2))=1;
Bord(1,1:size(Bord,2))=1;
Bord(size(Bord,1),1:size(Bord,2))=1;

for i = 1: size(NThole)
    Overlap=Bord & NThole==i;
    if ~sum(Overlap(:))
        NTfill(NThole==i)=1;
    end
end
image(NTfill*1000),pause(Pt)



%% Get Dat
PixelArea=sum(NTfill(:));

[y x] = find(NTfill>0);
V=[y x];
[Co,Sc]=princomp(V);  %find principal components of object
minVfirst=find(Sc(:,1)==min(Sc(:,1)),1); %find position of min extreme position
maxVfirst=find(Sc(:,1)==max(Sc(:,1)),1); %find position of max extreme position
minVsecond=find(Sc(:,2)==min(Sc(:,2)),1); %find position of min extreme position
maxVsecond=find(Sc(:,2)==max(Sc(:,2)),1); %find position of max extreme position


NTpc=NTfill*0;
NTpc(V(minVfirst,1),V(minVfirst,2))=1;
NTpc(V(maxVfirst,1),V(maxVfirst,2))=1;
NTpc(V(minVsecond,1),V(minVsecond,2))=1;
NTpc(V(maxVsecond,1),V(maxVsecond,2))=1;

image(NTfill*40+NTpc*200),pause(Pt)

Length=max(Sc(:,1))-min(Sc(:,1));
Width=max(Sc(:,2))-min(Sc(:,2));

Convex.PixelArea=PixelArea;
Convex.Area = PixelArea * (xyum)^2; % axon territory in micron^2
Convex.Length=Length * xyum;
Convex.Width=Width *xyum;

save([TPN 'data/Area'], 'Convex') 



% NTcol=zeros(ys,xs,3,'uint8');
% NTcol(:,:,1)=double(NeuMax)*255/double(max(NeuMax(:)))*2;
% NTcol(:,:,2)=NTmax*0;
% NTcol(:,:,3)=NTperim*300;
% 
% image(NTcol),pause(.1)






