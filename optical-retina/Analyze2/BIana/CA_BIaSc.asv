%%Uses OFONNa instead of OFFON so that there is no top or bottom limit on
%%cell arbors


'Get Cells'
load(['UseCells.mat'])
[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);
GetBIs=find(BIs);
BIAge=Age(BIs);
'Cells Gotten'

KPN = '\\128.208.64.36\wonglab\Josh\Analyzed\Output2\';

for k = 1:size(GetBIs,1)
    TPN =  char(UseCells(GetBIs(k))); 
    Slash=find(TPN=='\')';
    CellName=char(TPN(Slash(size(Slash,1)-1)+1:Slash(size(Slash,1))-1));
    
     clear CA BiDotDist DotShow DotFilt DendFilt DotMap DendMap DPos 
     clear Cell Red Green Blue TerFill TerFilt TerPerim TerRay Territory
     clear Mids sMids NN sNN
    if exist([TPN 'find\SG.mat'])
      
    
%% load info
'Load Data'
clear Use DPos Cent Mids Length Top Bottom
load([TPN 'UseSc.mat'])
load([TPN 'ONOFFa.mat'])

xyum=Use.Res;

cmap=(jet(256));
cmap(1,:)=0;
colormap(cmap)

set(gcf,'Colormap',cmap)


for B=1:2

DPos=Use.DPos(ONOFF.DotLayer==B,:);
Cent=Use.Cent;

Mids=Use.Mids(ONOFF.MidLayer==B,:);
Length=Use.Length(ONOFF.MidLayer==B,:);

NN=Use.NN(ONOFF.DotLayer==B,:);

Idir=dir([TPN 'mask']);
IMtest=imread([TPN 'mask\' Idir(3).name]);
[yM xM]=size(IMtest);
yM=fix(yM*xyum)+1;
xM=fix(xM*xyum)+1;
Msize=[yM xM];

%% create Nearest Node list (nearest node for each dot and distance)
'Sort Data'


% Get image properties

%fix stray data
    DPos(DPos<1)=1;
    Mids(Mids<1)=1;
    NN(NN<1)=1;
    DPos(DPos(:,1)>yM,1)=yM;
    Mids(Mids(:,1)>yM,1)=yM;
    NN(NN(:,1)>yM,1)=yM;
    DPos(DPos(:,2)>xM,2)=xM;
    Mids(Mids(:,2)>xM,2)=xM;
    NN(NN(:,2)>xM,2)=xM;
    
  

%%
  
    clear sMids sNN

    sMids=Mids;
    sNN=NN;
        
    sMids=round(sMids);
    sNN=round(sNN);

    %% Draw maps
    'Draw Maps'

    %%Draw Dot map
    DotMap=zeros(Msize(1),Msize(2));
    for i=1:size(sNN,1)
        DotMap(sNN(i,1),sNN(i,2))=DotMap(sNN(i,1),sNN(i,2))+1;
    end
    
    %%Draw Dend map
    DendMap=zeros(Msize(1),Msize(2));
    for i=1:size(sMids,1)
        DendMap(sMids(i,1),sMids(i,2))=DendMap(sMids(i,1),sMids(i,2))+Length(i);
    end

%% Filter Results
    
    %% make distance filter

    AreaS=10;
    Disk=fspecial('disk',AreaS);
    DotFilt=imfilter(DotMap,Disk,'same');
    DendFilt=imfilter(DendMap,Disk,'same');
    BlankFilt=imfilter(double(DendMap>-1000),Disk,'same');
    DotFilt=DotFilt./BlankFilt;
    DendFilt=DendFilt./BlankFilt;

    DendSkel=DendFilt;
    DendSkel(DendMap>0)=0;
    image(DendSkel*1000)

    %DDSkel=DotFilt./DendFilt;
    %DDSkel(DendMap>0)=0;
    %image(DDSkel*500)

    DotSkel=DotFilt;
    DotSkel(DendMap>0)=0;
    image(DotSkel*2000)


%% Find territory
    image(DendMap*100), pause(.3)

    Disk2=fspecial('disk',5);%min(5,min(AreaS/2,1)));
    TerFilt=imfilter(DendMap,Disk2,'same');

    DFlab=bwlabel(TerFilt);
    for i = 1:max(DFlab(:))
        lSize(i)=size(find(DFlab==i),1); %#ok<AGROW>
    end
    TerFilt=DFlab*0;
    TerFilt(DFlab==find(lSize==max(lSize)))=1;
    image(TerFilt*200),pause(.3)
    clear DFlab lSize
    
    %Close
    Csize=round(AreaS/2);
    SE=strel('disk',Csize);
    Buf=Csize*2;
    [tys txs] = size(TerFilt);
    BufT=zeros(tys+2*Buf,txs+2*Buf);
    BufT(Buf+1:Buf+tys,Buf+1:Buf+txs)=TerFilt;
    BufC=imclose(BufT,SE);
    TerC=BufC(Buf+1:Buf+tys,Buf+1:Buf+txs);
    image((TerFilt+TerC)*100), pause(.3)
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
    image(TerFill*100),pause(.01)
    clear Frame TerHole Surround 
    
    %Get Perimeter
    TerPerim=bwperim(TerFill,8);
    image((TerFill+TerPerim)*100),pause(.3)

    %DrawRays
    [tpy tpx]=find(TerPerim);
    TerRay=TerFill*0;

    for i = 1: size(tpy,1)
        ydif=Cent(1)-tpy(i); xdif=Cent(2)-tpx(i);
        Long=sqrt(ydif^2 + xdif^2);
        ystep=ydif/Long;
        xstep=xdif/Long;

        for l = 1:.5:Long
            TerRay(round(tpy(i)+ystep*l),round(tpx(i)+xstep*l))=1;
        end
    end
    %image(TerRay*100),pause(.01)
    
    Territory=TerFill;
    DotDist=DotSkel.*Territory;
    image(DotDist*2000)
    
%%   Collect Data

 
    CA.Arbor(B).DotMap=DotMap;
    CA.Arbor(B).DendMap=DendMap;
    CA.Arbor(B).Territory=Territory;
    CA.Arbor(B).DotDist=DotFilt.*Territory;
    CA.Arbor(B).DendDist=DendFilt.*Territory;
    CA.Arbor(B).TerPerim=TerPerim;
    CA.Arbor(B).MeanDotOArea=median(DotFilt(Territory>0));
    

    CA.mDtA(B)=CA.Arbor(B).MeanDotOArea;
    CA.mDdA(B)=median(CA.Arbor(B).DendDist(CA.Arbor(B).Territory>0));
    CA.mDDA(B)=CA.mDtA(B)/CA.mDdA(B);   
    CA.Area(B)=sum(CA.Arbor(B).Territory(:));
    
    CA.WholeA.DD(B)=sum(DotMap(:))/sum(DendMap(:));
end

end % for B = 1: 2 run both arbors


%%Write images
clear Cell

   ScaleRG=max(max(CA.Arbor(1).DotDist(:)),max(CA.Arbor(2).DotDist(:)));
   ScaleRG=255/ScaleRG;
   Red=CA.Arbor(1).DotDist*ScaleRG;
   Red=Red+CA.Arbor(1).TerPerim*75;
   Green=CA.Arbor(2).DotDist*ScaleRG;
   Green=Green+CA.Arbor(2).TerPerim*75;
   Blue=~(CA.Arbor(2).Territory & CA.Arbor(1).Territory);
   BiDotDist(:,:,1)=uint8(Red);
   BiDotDist(:,:,2)=uint8(Green);
   BiDotDist(:,:,3)=uint8(Blue)*0;
   CA.Im.BiDotDist=BiDotDist;
   image(BiDotDist),pause(.3)
   
   if exist('Cell')
    cName=['P' Cell.Age '_' Cell.Type '_' CellName '.tif']
    Name=[KPN 'Images\BiCom\Bi' cName];
    else
        Name=[KPN 'Images\BiCom\Bi' CellName '.tif'];
    end
    if isempty(find(Name=='?'));
        imwrite(BiDotDist,Name,'Compression','none')
    end

   
    save([TPN 'CAbiaSc.mat'],'CA')
    
    for a = 1:2
      [ysA(a) xsA(a)]=size(CA.Arbor(a).DotMap);
    end
      ys=max(ysA); xs=max(xsA);
      CA4=zeros(ys*2,xs*4);
      ys=212; xs=212;
      for B = 1:2
          ystart=(B-1)*ys+1;
          DotDist=CA.Arbor(B).DotDist;
          DendDist=CA.Arbor(B).DendDist;
          DotMap=CA.Arbor(B).DotMap;
          DendMap=CA.Arbor(B).DendMap;
          DDMap=DendMap*3000000;
          DDMap(DotMap>0)=100;
          
          
          CA4(ystart:ystart+ysA(B)-1,1:xsA(B))=DDMap;
          CA4(ystart:ystart+ysA(B)-1,xs+1:xs+xsA(B))=DendDist*300;
          CA4(ystart:ystart+ysA(B)-1,xs*2+1:xs*2+xsA(B))=DotDist*1500;
          CA4(ystart:ystart+ysA(B)-1,xs*3+1:xs*3+xsA(B))=DotDist./DendDist*300;
          image(CA4),pause(.3)
      end
      
      
      if exist([TPN 'Cell.mat'])
            load([TPN 'Cell.mat'])
            cName=['P' Cell.Age '_' Cell.Type '_' CellName '.tif']
            Name=[KPN 'Images\BiCom\C4\' cName];
        else
            Name=[KPN 'Images\BiCom\C4\' CellName '.tif'];
        end
        if isempty(find(Name=='?'));
            imwrite(CA4,cmap,Name,'Compression','none')
        end
      
    

pause(.3)
end % Run all cells

     