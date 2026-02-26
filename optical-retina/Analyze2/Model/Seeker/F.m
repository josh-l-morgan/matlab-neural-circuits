
fsize=30
MnVal=.2;

Wall=zeros(fsize,fsize);
Wall(1:fsize,[1 fsize])=1;
Wall([1 fsize],1:fsize)=1;
field=rand(fsize,fsize)>.5;
TargF=rand(fsize,fsize)>0.95;
TargF((Wall |field)>0)=0;

F.mem=field*0;  %initialize memory
F.T=0;
F.F=0;
F.M=0;


[F.y F.x]= find(~(Wall | field),1);
Ff=field*0;
Ff(F.y,F.x)=1;


Show=uint8(TargF*200);
Show(:,:,2)=Ff*200;
Show(:,:,3)=Wall*10000 + field*150;


image(Show),pause(.1)

reps = 10;
side=[0 1 ; 0 -1; -1 0; 1 0];
for rep = 1:reps
   rep, pause
    %%look 
   neary=side(:,1)+F.y;
   nearx=side(:,2)+F.x;
   inear=sub2ind([fsize fsize],neary,nearx);
   Wn=Wall(inear);
   Tn=TargF(inear);
   Fn=field(inear);
   Mn=F.mem(inear);
   Val=Tn-Mn*MnVal-Fn-Wn*1000;
   ToGo=find(Val==max(Val),1);
   Goy=neary(ToGo); Gox=nearx(ToGo);
   
   
   %%move
   F.y=Goy; F.x=Gox;
   F.T = F.T +Tn(ToGo);
   F.F=F.F+Fn(ToGo);
   F.M=F.M+Mn(ToGo);
   
   %% change spot
   TargF(Goy,Gox)=0;
   field(Goy,Gox)=0;
   F.mem(Goy,Gox)=F.mem(Goy,Gox)+1;
   
   %% Show
    Show=uint8(TargF*200);
    Ff=field*0;
    Ff(F.y,F.x)=1;
    Show(:,:,2)=Ff*200;
    Show(:,:,3)=Wall*10000 + field*150;
    
    image(Show), pause(.01)

end