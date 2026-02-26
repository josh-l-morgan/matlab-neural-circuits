%% Define first time point unablated territories minus given cell
%% Select folder containing ablation experiment with the first time point
%% in PreUn\I

clear all
colormap gray(256)
GPN=GetMyDir
GPNf=[GPN 'PreUn\I'];
GPNd=dir(GPNf);
GPNd=GPNd(3:size(GPNd,1));
for i = 1: size(GPNd,1)
    AblatedLocation=[GPNf '\' GPNd(i).name];
    AbRegi=imread(AblatedLocation);
    AbRegi=AbRegi>0;
    AbRegi=max(AbRegi,[],3);
    AbReg(:,:,i)=AbRegi;
end
image(sum(AbReg,3)*60)

Names={};

%% Load ablated
APNf=[GPN 'Combined\Ter'];
APNd=dir(APNf);
AblatedLocation=[APNf '\' APNd(3).name];
Ab=imread(AblatedLocation);
Ab=Ab>0;
Ab=max(Ab,[],3);


%% Run all
for k = 1 : size(AbReg,3)
if ~exist([GPN 'UnAbOthers']), mkdir([GPN 'UnAbOthers']),end    
Names(k,1)={GPNd(k).name};    
I=AbReg(:,:,(1:size(AbReg,3))~=k);
%'Raw Image'

image(sum(I,3)*100),pause(.1)
It=sum(I,3);
It=It>0;  %define binary image
%'Binary Image'
image(It*1000),pause(.1)
%% Filter Results
    
    %% make distance filter

    CloseRad=40;
    FiltRad=2;
    Disk=fspecial('disk',CloseRad*2);


%% Find territory
    Disk2=fspecial('disk',FiltRad);%min(5,min(FiltRad,1)));
    TerFilt=imfilter(uint8(It*1000),Disk2,'same');
    TerFilt=TerFilt>0;
    %{
    %%Eliminate disconnected regions
    DFlab=bwlabel(TerFilt);
    for i = 1:max(DFlab(:))
        lSize(i)=size(find(DFlab==i),1); %#ok<AGROW>
    end
    TerFilt=DFlab*0;
    TerFilt(DFlab==find(lSize==max(lSize)))=1;
    %}
    %'TerFilt'
    image(TerFilt*200),pause(.1)
    clear DFlab lSize
    
    %Close
    Csize=round(CloseRad);
    SE=strel('disk',Csize);
    Buf=Csize*2;
    [tys txs] = size(TerFilt);
    BufT=zeros(tys+2*Buf,txs+2*Buf);
    BufT(Buf+1:Buf+tys,Buf+1:Buf+txs)=TerFilt;
    BufC=imclose(BufT,SE);
    TerC=BufC(Buf+1:Buf+tys,Buf+1:Buf+txs);
    %'TerC = closed image'
    image((TerFilt+TerC)*100), pause(.1)
    clear Buf BufT BufC
    
    %Remove holes
    TerFill=TerC;
    Frame=TerC*0; Frame(1:tys,1)=1; Frame(1:tys,txs)=1; Frame(1,1:txs)=1; Frame(tys,1:txs)=1;
    TerHole=bwlabel(~TerC);
    for i=1:max(TerHole(:))
        Surround=sum(sum((TerHole==i) .* Frame));
        AblateTer=sum((TerHole==i) .* Ab);
        if ~Surround & ~AblateTer
            TerFill(TerHole==i)=1;
        end        
    end
    %'TerFill'
    image(TerFill*100),pause(.1)
    clear Frame TerHole Surround 
    
    %Get Perimeter
    TerPerim=bwperim(TerFill,8);
    %'TerPerim'
    image((TerFill+TerPerim)*100),pause(.1)

    %DrawRays
    [tpy tpx]=find(TerPerim);
    TerRay=TerFill*0;

    %% find center of mass
    [y x]=find(It);
    Cent=[mean(y) mean(x)];
    
    for i = 1: size(tpy,1)
        ydif=Cent(1)-tpy(i); xdif=Cent(2)-tpx(i);
        Long=sqrt(ydif^2 + xdif^2);
        ystep=ydif/Long;
        xstep=xdif/Long;

        for l = 1:.5:Long
            TerRay(round(tpy(i)+ystep*l),round(tpx(i)+xstep*l))=1;
        end
    end
    %'TerRay'
    image(TerRay*100),pause(.1)
    
    
    Territory=TerFill;
    %'Territory'
    image(Territory*1000),pause(.1)
    imwrite(Territory,[GPN 'UnAbOthers\Not_' GPNd(k).name ],'Compression','none')
    
    Ters(:,:,k)=Territory>0;

end












     