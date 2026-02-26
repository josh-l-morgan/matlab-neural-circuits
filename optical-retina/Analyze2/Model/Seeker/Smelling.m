function[Smells,Dists]=Smelling(TargF, field, Wall);

fsize=size(TargF,1)
[ty tx] = find(TargF);
Smells=zeros(size(TargF,1),size(TargF,2),size(ty,1));
side=[0 1 ; 0 -1; -1 0; 1 0];


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
            Smells(neary(Pn),nearx(Pn),t)=1/s;
            nexty=[nexty ;neary(Pn)];
            nextx=[nextx ;nearx(Pn)];
        end
        starty=nexty; startx = nextx;
        %image(Smells(:,:,t)*100),pause(.01)
        if isempty(nexty),break,end
    end
    
end
    
Dists = Smells*0;
   for y = 1:fsize
       for x = 1:fsize
            Dists(y,x,:)=sqrt((y-ty).^2+(x-tx).^2);
       end
   end
Dists=1/Dists;
   