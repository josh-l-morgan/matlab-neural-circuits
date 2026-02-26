colormap gray(100)
dur=10000;
Type = 5;

fsize=201;
field=zeros(fsize)+100;
Change=ones(fsize)> 0;
Change(round(fsize/2),round(fsize/2))=0;
Cycle=30;
Blank=20;
Bright=90;
field=field*0+Bright/2;

if Type == 1

Rate=.1;
for i = 1:dur
    image(sin(i*Rate)*50+50),pause(.01)
end

elseif Type == 2

for i = 1:dur
   val=(mod(i,Cycle)>Blank)*Bright;
   field(Change)=val;
   image(field),pause(.01)
    
end

elseif Type == 3
    %%Drift
    Rad=round(fsize/10);
    Sig=Rad/5;
    Rfield=rand(fsize)*100;
    SE=fspecial('gaussian',Rad,Sig);
    ChangeSpeed=10;
    for i = 1: dur
    
        Cfield=randn(fsize)*ChangeSpeed;
        Rfield=Rfield+Cfield;
        Ffield=imfilter(Rfield,SE,'same');
        image(Ffield),pause(.01)
    
    end
elseif Type ==4
    %%color Drift
    Rad=round(fsize/10);
    Sig=Rad/5;
    Rfield=rand(fsize)*100;
    Bfield=rand(fsize)*100;
    Gfield=rand(fsize)*100;
    SE=fspecial('gaussian',Rad,Sig);
    ChangeSpeed=10;
    Afield=zeros(fsize,fsize,3,'uint8');
    for i = 1: dur
    
        Cfield=randn(fsize)*ChangeSpeed;
        Rfield=Rfield+Cfield;
        RFfield=imfilter(Rfield,SE,'same');
        
        Cfield=randn(fsize)*ChangeSpeed;
        Bfield=Bfield+Cfield;
        BFfield=imfilter(Bfield,SE,'same');
        
        Cfield=randn(fsize)*ChangeSpeed;
        Gfield=Gfield+Cfield;
        GFfield=imfilter(Gfield,SE,'same');
        
        Afield(:,:,1)=RFfield; Afield(:,:,2)=GFfield; 
        Afield(:,:,3)=BFfield;
        image(Afield*2.55),pause(.01)
    
    end
elseif Type==5
    %%Shift
    Rad=round(fsize/10);
    Sig=Rad/5;
    Rfield=rand(fsize)*100;
    SE=fspecial('gaussian',Rad,Sig);
    ChangeSpeed=3;
    SatStd=2.5;  %Number of standard deviations to saturate
    Ffield=imfilter(Rfield,SE,'same','circular');
    meanSc=mean(Ffield(:));
    stdSc=std(Ffield(:));
    for i = 1: dur
    
        Cfield=randn(fsize)*ChangeSpeed;
        Rfield=Rfield+Cfield;%Cfield/10+Cfield.*abs(Rfield-50)/50;
        Rfield(Rfield<0)=0;
        Rfield(Rfield>100)=100;
        Ffield=imfilter(Rfield,SE,'same','circular');
        subplot(4,1,1:3)
        image((Ffield-meanSc)*stdSc*SatStd+meanSc)
        subplot(4,1,4)
        hist(Ffield(:),0:1:100);
        pause(.01)
    end
end


