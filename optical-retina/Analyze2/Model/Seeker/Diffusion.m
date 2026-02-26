
[ty tx] = find(TargF);
Diffs=zeros(size(TargF,1),size(TargF,2),size(ty,1));
side=[0 1 ; 0 -1; -1 0; 1 0];


for t = 1 : size(ty,1) % run all targets
    
    starty=ty(t); startx=tx(t);
    Path = ~(field | Wall );
    Diffs(starty(1),startx(1),t)=1000;
    for s = 1:fsize*2;  %run 100 itterations for each target
        
        for n = 1: size(starty,1) %run each pixel in wave
           
            neary=side(:,1)+starty(n);
            nearx=side(:,2)+startx(n);
            inear=sub2ind([fsize fsize],neary,nearx);
            Pn=Path(inear);
            Wn=field(inear);
            %Path(inear)=0;
            Conc=Diffs(starty(n),startx(n),t)
            Spread=(sum(Pn)+1)*10 + sum(Wn);
            Diffs(starty(n),startx(n),t)=10*(Conc/Spread);
            Loc=sub2ind([fsize fsize size(ty,1)], neary(Pn), nearx(Pn),ones(sum(Pn),1)*t);
            Diffs(Loc)= Diffs(Loc)+10*(Conc/Spread);
            Loc=sub2ind([fsize fsize size(ty,1)], neary(Wn), nearx(Wn),ones(sum(Wn),1)*t);
            Diffs(Loc)= Diffs(Loc)+(Conc/Spread);
            
            
        end
        [starty startx] = find(Diffs(:,:,t));
        image(Diffs(:,:,t)*5000/max(max(Diffs(:,:,t)))),pause(.01)
    end
    
end
    
