%%Read OIFs and assemble into single matrix that can be written as 2D
%%color tifs

IPN=GetMyDir;
IPNd=dir(IPN); IPNd=IPNd(3:size(IPNd,1));

%%Identify Tiffs
TiffNames={};
for i = 1: size(IPNd,1)
    siz=length(IPNd(i).name);
    if IPNd(i).name(siz-2:siz)== 'tif'
        TiffNames(length(TiffNames)+1,1)={IPNd(i).name}
    end
end


%%Identify image dims
LastName=cell2mat(TiffNames(length(TiffNames)));
siz=length(LastName);
Channels=str2num(LastName(siz-10:siz-8));
Planes=str2num(LastName(siz-6:siz-4));
It=imread([IPN LastName]);
[ys xs]=size(It);


%%Read image
I =zeros(ys,xs,3,Planes);
for i = 1:size(TiffNames,1)
    Name=cell2mat(TiffNames(i));
    Channel=str2num(Name(siz-10:siz-8));
    Plane=str2num(Name(siz-6:siz-4));
    I(:,:,Channel,Plane)=imread([IPN Name]);
end


%%Scale to 8bit
for i = 1:3
    minc=min(min(min(I(:,:,i,:))));
    I(:,:,i,:)=I(:,:,i,:)-minc;
    maxc=max(max(max(I(:,:,i,:))));
    I(:,:,i,:)=double(I(:,:,i,:))*256/maxc;    
end
I=uint8(I);
Imax=max(I,[],4);
image(Imax*3);


