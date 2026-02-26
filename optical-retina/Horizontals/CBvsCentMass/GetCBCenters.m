%%Retrieves the centers of masses from individually stored tifs 
%%and relates them to reference point
clear all
TPN = GetMyDir;

TPNi= [TPN 'I\'];
TPNt= [TPN 'Ter\'];

%% Run Image (cb)
TPNd=dir(TPNi);TPNd=TPNd(3:length(TPNd));
Names={};
for i = 1:length(TPNd)
    nam=TPNd(i).name;
    ln=length(nam);
    if nam(ln-3:ln)=='.tif'
        Names{length(Names)+1}=nam;
    end
end

Cent=zeros(length(Names),2);
clear Ia
for i = 1: length(Names)
    Is=imread([TPNi char(Names{i})]);
    Is=sum(Is,3);
    
    %% Find CB
    SE = strel('disk',3);
    If=Is>0;
    for e = 1:1000
        Ib=If;
        If=imerode(If,SE);
        %image(If*200),pause(.01)
        
    if max(If(:))==0
        [y x]=find(Ib);
        CentCB(i,:)=[mean(y) mean(x)];
        break
    end

    end
    image(double(Is)*.4+Ib*100),pause(.01)
    
    if ~exist('Ia2'),Ia2=Is*0; end
    Ia2(Is>0)=Ia2(Is>0)+1;
    [y x z] = size(Is);
    siz= [y x];
    [y x]=find(Is>0);
    CentR(i,:)=[mean(y) mean(x)];

end


Reference=[   463.5809  466.4296 ];

DistsCB=dist(CentCB,Reference);

%% Run TErs

TPNd=dir(TPNt);TPNd=TPNd(3:length(TPNd));
Names={};
for i = 1:length(TPNd)
    nam=TPNd(i).name;
    ln=length(nam);
    if nam(ln-3:ln)=='.tif'
        Names{length(Names)+1}=nam;
    end
end


Cent=zeros(length(Names),2);
clear Ia
for i = 1: length(Names)
    Is=imread([TPNt char(Names{i})]);
    if ~exist('Ia'),Ia=Is*0; end
    Ia(Is>0)=Ia(Is>0)+1;
    [y x z] = size(Is);
    siz= [y x];
    [y x]=find(Is>0);
    CentTer(i,:)=[mean(y) mean(x)];
    image(Is),pause(.01)

    %% Store area
    Area(i,1)=length(find(Is>0));    

end

DistsTer=dist(CentTer,Reference);


Vals=Ia(Ia>0);
[y x] = find(Ia);
Scaled= [y x] .* double([Vals Vals]);
CentAll=sum(Scaled,1)/sum(Vals);

ER=sqrt(Area/pi); % Find Effective Radius

DistsCB2Ter=sqrt((CentTer(:,1)-CentCB(:,1)).^2 + (CentTer(:,2)-CentCB(:,2)).^2);

OffSet=DistsCB2Ter./ER

Direct=(DistsCB-DistsTer)./ER;


%%


Ic=zeros(size(Ia,1),size(Ia,2),3,'uint8');
Ic(:,:,3)=uint8(Ia2)*100+uint8(Ia)*50;
for i = 1:size(Cent,1)
Ic(round(CentTer(i,1)),round(CentTer(i,2)),2)=255;
Ic(round(CentCB(i,1)),round(CentCB(i,2)),1)=255;
end

image(Ic),pause(.01)


DatTitle={ 'ER' 'OffSet' 'DistsCB' 'DistsTer' 'Direct'};
Dat=[ER OffSet DistsCB DistsTer Direct];
DatTitle(2:size(Dat,1)+1,:)=mat2cell(Dat,ones(size(Dat,1),1),ones(size(Dat,2),1))


























