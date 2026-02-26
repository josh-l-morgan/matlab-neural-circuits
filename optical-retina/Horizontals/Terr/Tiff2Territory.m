clear all
KPN=GetMyDir
TPN=[KPN 'I\']
Idir=dir(TPN); Idir=Idir(3:size(Idir));
colormap gray(255)



Names={};
for k = 1 : size(Idir,1)
if ~exist([KPN 'Ter']), mkdir([KPN 'Ter']),end    
Names(k,1)={Idir(k).name};    
I=imread([TPN Idir(k).name]);
%'Raw Image'

image(I),pause(.1)
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
        if ~Surround
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
    
    imwrite(Territory,[KPN 'Ter\Ter_' Idir(k).name ],'Compression','none')
    
    Ters(:,:,k)=Territory>0;

end

clear Area Overlap Dat

SumTers=sum(Ters,3);
image(SumTers*50),pause(.01)

Dat(1,1)={'Area'};
Dat(1,2)={'Overlap'};
Dat(1,3)={'Private'};

for i = 1:size(Ters,3)
   Area(i)=sum(sum(Ters(:,:,i))) ;
   
   Overlap(i)=sum(sum(SumTers(Ters(:,:,i)>0))) - Area(i);
   Private(i)=sum(sum(SumTers(Ters(:,:,i)>0)==1)); 
   
   Dat(i+1,1)={Area(i)};
   Dat(i+1,2)={Overlap(i)};
   Dat(i+1,3)={Private(i)};
   
end
Dat(2:size(Names,1)+1,4)=Names;
Dat
xlswrite([KPN 'Dat.xls'],Dat)













     