%% Extract surround
TPN=GetMyDir;

load([TPN 'Post.mat'])
load([TPN 'Dend.mat'])
load([TPN 'Colo.mat'])
Chans={'Post','Dend','Colo'};


[ys xs zs] = size(Post);

load([TPN 'ImageInfo.mat'])
load([TPN 'Use.mat'])

xyum=ImageInfo{4};
zum=ImageInfo{5};

VoxP=Use.DPos;
VoxP(:,1:2)=VoxP(:,1:2)/xyum;
VoxP(:,3)=VoxP(:,3)/zum;
VoxP=round(VoxP);

Side=5;  %distance in um to look
Sidxy=round(Side/xyum);
Sidz=round(Side/zum);
VoxP(:,1:2)=VoxP(:,1:2)+Sidxy;
VoxP(:,3)=VoxP(:,3)+Sidz;
BuffS=zeros(ys+Sidxy*2,xs+Sidxy*2,zs+Sidz*2,'uint8');
BuffS(Sidxy+1:ys+Sidxy,Sidxy+1:xs+Sidxy,Sidz+1:zs+Sidz)=1;
BuffI=zeros(ys+Sidxy*2,xs+Sidxy*2,zs+Sidz*2);


for c = 1 : 3

BuffI(Sidxy+1:ys+Sidxy,Sidxy+1:xs+Sidxy,Sidz+1:zs+Sidz)=eval(Chans{c});

Samp=zeros(Sidxy*2+1,Sidxy*2+1,Sidz*2+1);
Sum=Samp;
for i = 1:size(VoxP,1)
    Samp=Samp+double(BuffI(VoxP(i,1)-Sidxy:VoxP(i,1)+Sidxy,VoxP(i,2)-Sidxy:VoxP(i,2)+Sidxy,...
        VoxP(i,3)-Sidz:VoxP(i,3)+Sidz));
    Sum=Sum+double(BuffS(VoxP(i,1)-Sidxy:VoxP(i,1)+Sidxy,VoxP(i,2)-Sidxy:VoxP(i,2)+Sidxy,...
        VoxP(i,3)-Sidz:VoxP(i,3)+Sidz));
end

Ave=Samp./Sum;


Near{c}=Ave;
end

%% Draw

SampP=Near{1};
SampD=Near{2};
SampC=Near{3};
mSampP=max(SampP,[],3);
mSampD=max(SampD,[],3);
mSampC=max(SampC,[],3);

mSampA=mSampP*255/double(max(mSampP(:)));
mSampA(:,:,2)=mSampD*255/double(max(mSampD(:)));
mSampA(:,:,3)=mSampC*255/double(max(mSampC(:)));

image(uint8(mSampA))
sAve=max(Ave,[],3);
image(sAve*255/max(sAve(:)))

surf(sAve)
colormap(jet)
shading interp







