%%inputs a label field and spits out a 2D description of the arbor of each
%%cell
%%Assumes that chosen directory contains folder 'I' that contains 
%2D tiffs that correspond to masks of each cell

clear all
TPN=GetMyDir
TPNi=[TPN 'I\'];
Idir=dir(TPNi); Idir=Idir(3:size(Idir));
colormap gray(255)

Names={};
for i = 1:length(Idir)
    nam=Idir(i).name;
    LN=length(nam);
    if nam(LN-3:LN)=='.tif'
        Names{length(Names)+1}=nam;
    end
end

Planes=length(Names);
for i = 1:Planes
    I(:,:,i)=imread([TPNi Names{i}]);
end

Cells=unique(I(:));
Cells=Cells(Cells>0);
[ys xs zs]=size(I);
siz=[ys xs zs];
numCells=length(Cells);

Masks=zeros(ys,xs,length(Cells));
for c = 1: numCells
 D2cell=zeros(siz(1:2));
 [y x z] =  ind2sub(siz,find(I==Cells(c)));
 D2ind=sub2ind(siz(1:2),y,x);
 D2cell(D2ind)=1;  %define binary image
 Masks(:,:,c)=D2cell;
end

save([TPN 'Masks.mat'],'Masks')

%% Find all cell territories
for c = 1 : numCells
    It=Masks(:,:,c);
%'Binary Image'
image(It*1000),pause(.1)
%%Filter Results
    
    %% make distance filter
    Rsc=.15/.23;
    
    CloseRad=40*Rsc;
    FiltRad=2;
    Disk=fspecial('disk',CloseRad*2);


%%Find territory
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
      
    Ters(:,:,c)=uint8(Territory>0)*Cells(c);

end

save([TPN 'Ters.mat'],'Ters')


%%Retrieves the centers of masses from individually stored tifs
%%and relates them to reference point


%% Get all directories
% Got=0; TPNa={};
% while 1
%     Got=GetMyDir;
%     if strcmp(Got(2),'\'),break,end
%     TPNa{length(TPNa)+1}=Got;
% end

%% Remember

xyum=0.153;


%% calculate pixel by dist
Look = 250;
LookDs=.5:1:Look*xyum+.5;
fsize=Look*2+2;
field=zeros(fsize,fsize);
mid = [Look+1 Look+1];
[yp xp]=ind2sub([fsize fsize],find(field==0));
DPos=dist([yp xp],mid);
DPos=DPos*xyum;
Pixs=hist(DPos,LookDs);
%Pixs=Pixs(1:length(Pixs)-1);

CountTer=0; CountCB = 0;

    CentCB=zeros(numCells,2);

    for i = 1: numCells
        
        Is=Masks(:,:,i);
        
        %% Find CB
        SE = strel('disk',3);
        If=Is>0;
        for e = 1:1000
            Ib=If;
            If=imerode(If,SE);
            %image(If*200),pause(.01)

            If(1,:)=0; If(size(If,1),:)=0;
            If(:,1)=0; If(:,size(If,2))=0;
            if max(If(:))==0
                [y x]=find(Ib);
                CentCB(i,:)=[mean(y) mean(x)];
                break
            end
        end
        image(double(Is)*.4+Ib*100),pause(.01)


        %%  DRP
        CB=round(CentCB(i,:));
        [y x] = find(Is>0);
        Dists=dist([y x], CB);
        Dists=Dists*xyum;
        CountCB=CountCB+1;
        MaskHist(CountCB,:)=hist(Dists,LookDs);

    end


    Cent=zeros(length(Names),2);
    clear Ia CentTer
    for i = 1: numCells
        Is=Ters(:,:,i);
        if ~exist('Ia'),Ia=Is*0; end
        Ia(Is>0)=Ia(Is>0)+1;
        [y x z] = size(Is);
        siz= [y x];
        [y x]=find(Is>0);
        CentTer(i,:)=[mean(y) mean(x)];
        image(Is),pause(.01)


        %%  DRP
        CB=round(CentCB(i,:));
        [y x] = find(Is>0);
        Dists=dist([y x], CB);
        Dists=Dists*xyum;
        CountTer=CountTer+1;
        TerHist(CountTer,:)=hist(Dists,LookDs);
        Area(CountTer,1)=length(y)*xyum^2;


    end




ER=sqrt(Area/pi); % Find Effective Radius
DenTer=mean(TerHist,1)./Pixs;
DenMask=mean(MaskHist,1)./Pixs;


bar(DenTer,'b'), hold on
bar(DenMask,'r')
hold off

ArborProfile.ER=ER;
ArborProfile.DenTer=DenTer;
ArborProfile.DenMask=DenMask;

save([TPN 'ArborProfile.mat'],'ArborProfile')
load([TPN 'ArborProfile.mat'])


CellStats=double([]);
CellStats(:,1)=Cells;  %
CellStats(:,4:5)=CentTer;  %in pixels
CellStats(:,2:3)=CentCB;













     