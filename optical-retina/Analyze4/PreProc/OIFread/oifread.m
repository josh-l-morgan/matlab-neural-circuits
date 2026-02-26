function[I] = oifread(IPN)
%%Read OIFs and assemble into single matrix that can be written as 2D
%%color tifs.  Input is the directory that includes the tiffs

IPNd=dir(IPN); IPNd=IPNd(3:size(IPNd,1));

%%Identify Tiffs
TiffNames={};
for i = 1: size(IPNd,1)
    siz=length(IPNd(i).name);
    if IPNd(i).name(siz-2:siz)== 'tif'
        if IPNd(i).name(siz-8:siz-7)~='-R'
            TiffNames(length(TiffNames)+1,1)={IPNd(i).name};
        end
    end
end


%%Identify image dims
LastName=cell2mat(TiffNames(length(TiffNames)));
siz=length(LastName);

Cid=find(LastName=='C',1,'last');
Channels=str2num(LastName(Cid+1:Cid+3));
if isempty(Channels) | Channels==0, Channels=1; end

Pid=find(LastName=='Z',1,'last');
Planes=str2num(LastName(Pid+1:Pid+3));
if isempty(Planes) | Planes == 0; Planes=1; end

% Channels=str2num(LastName(siz-10:siz-8));
% Planes=str2num(LastName(siz-6:siz-4));

It=imread([IPN LastName]);
[ys xs]=size(It);


%%Read image
I =zeros(ys,xs,3,Planes);
for i = 1:size(TiffNames,1)
    Name=cell2mat(TiffNames(i));
    
    Cid=find(Name=='C',1,'last');
    Channel=str2num(Name(Cid+1:Cid+3));
    if isempty(Channel) | Channel == 1, Channel=1; end
    
    Pid=find(Name=='Z',1,'last');
    Plane=str2num(Name(Pid+1:Pid+3));
    if isempty(Plane) | Plane==0, Plane=1; end
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


