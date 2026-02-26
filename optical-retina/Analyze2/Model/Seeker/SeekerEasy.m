%% Set the weights to find the food breaking the fewest walls


%%Choice determined by
%%Val=Tn * TnVal - Mn*MnVal - Fn *FnVal-Wn*1000 + Sn * SnVal +Dn* DnVal;
MnVal = 1; % Value of avoiding retracing steps
FnVal = 5; % Value of avoiding breaking walls
DnVal = 10; % Value of following sound (through walls)

SnVal = 100; % Value of following nose (no wall breaking)
TnVal = 1000  % Value of grabbing target

%%shape field
OpenSpace= .6; % 0-1

while 1 == 1
fsize=30
'Building Field'
Wall=zeros(fsize,fsize);
Wall(1:fsize,[1 fsize])=1;
Wall([1 fsize],1:fsize)=1;
TargF=uint8(rand(fsize,fsize)>0.95);

Filt = [
    0.0751    0.1238    0.0751;
    0.1238    0.2042    0.1238;
    0.0751    0.1238    0.0751]

field=uint8(rand(fsize,fsize)>OpenSpace)*100;
field=filter2(Filt,uint8(field*100),'same');
field=field>median(field(:));
image(field*200)

TargF((Wall |field)>0)=0;

%% Map sound and smell
[ty tx] = find(TargF);
Smells=zeros(size(TargF,1),size(TargF,2),size(ty,1));
side=[0 1 ; 0 -1; -1 0; 1 0];

%%Smell
for t = 1 : size(ty,1) % run all targets
    
    starty=ty(t); startx=tx(t);
    Path = ~(field | Wall );
    for s = 1:100;  %run 100 itterations for each target
        nexty=[]; nextx = [];
        
        for n = 1: size(starty,1) %run each pixel in wave
           
            neary=side(:,1)+starty(n);
            nearx=side(:,2)+startx(n);
            inear=sub2ind([fsize fsize],neary,nearx);
            Pn=Path(inear);
            Path(inear)=0;
            Loc=sub2ind([fsize fsize size(ty,1)], neary(Pn), nearx(Pn),ones(sum(Pn),1)*t)
            Smells(Loc)=1/s;
            nexty=[nexty ;neary(Pn)];
            nextx=[nextx ;nearx(Pn)];
        end
        starty=nexty; startx = nextx;
        %image(Smells(:,:,t)*100),pause(.01)
        if isempty(nexty),break,end
    end
    
end

%%Sound
Dists = Smells*0;
   for y = 1:fsize
       for x = 1:fsize
            Dists(y,x,:)=sqrt((y-ty).^2+(x-tx).^2);
       end
   end
Dists=1/Dists;
   
Smell=sum(Smells,3);
Dist=sum(Dists,3);

image(Smell*300),pause(.1)
image(Dist*10),pause(.1)

%% 
[ty tx] = find(TargF);
for t = 1:size(ty,1)
    TargF(ty(t),tx(t))=t;
end
useT=ones(size(ty,1),1)>0;


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

%smell

image(Show),pause(.1)

reps = 10000;
side=[0 1 ; 0 -1; -1 0; 1 0];
'Play'
for rep = 1:reps
    %%look 
   neary=side(:,1)+F.y;
   nearx=side(:,2)+F.x;
   inear=sub2ind([fsize fsize],neary,nearx);
   Wn=Wall(inear);
   Tn=TargF(inear)>0;
   Fn=field(inear);
   Mn=F.mem(inear);
   Dn=Dist(inear);
   Sn=Smell(inear);
   Val=Tn * TnVal - Mn*MnVal - Fn *FnVal-Wn*1000 + Sn * SnVal +Dn* DnVal;
   ToGo=find(Val==max(Val),1);
   Goy=neary(ToGo); Gox=nearx(ToGo);
   
   
   %%move
   F.y=Goy; F.x=Gox;
   F.T = F.T +Tn(ToGo);
   F.F=F.F+Fn(ToGo);
   F.M=F.M+Mn(ToGo);
   
   %% change spot
   if TargF(Goy,Gox)
    useT(TargF(Goy,Gox))=0;
    F.mem=F.mem*0;
   end
   Smell=sum(Smells(:,:,useT),3);
   Dist=sum(Dists(:,:,useT),3);
   %image(Dist*20),pause(.01);   
   
   TargF(Goy,Gox)=0;
   field(Goy,Gox)=0;
   F.mem(Goy,Gox)=F.mem(Goy,Gox)+1;
   
   %% Show
    Show=uint8(TargF*200);
    Ff=field*0;
    Ff(F.y,F.x)=1;
    Show(:,:,2)=Ff*200;
    Show(:,:,1)=Show(:,:,1)+uint8(Ff)*200;
    Show(:,:,2)=Show(:,:,2)+uint8(F.mem)*30;
    Show(:,:,3)=Wall*10000 + field*150;
    
    image(Show), pause(.01)
    if ~sum(useT), break,end
    
end

F
WallsBroken=F.F
Time = rep

end
