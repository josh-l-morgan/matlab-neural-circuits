%%Vector, the friendly vector

Vbody=ones(5);
[Vy Vx]=find(Vbody);
field=zeros(20);

for i = 1:100
    
    
    for n=1:size(Vy,1)
       
       moves=[];
       for y = -1:1:1, for x = -1 : 1 :1
            d=sqrt((Vy-(Vy(n)+y)).^2+(Vx-(Vx(n)+x)).^2);
            ok=sum(d<1.5)>1;
            block=sum(d==0);
            if ok & ~block
                moves=[moves ; y x];
                
            end
            %d,y,x,ok,block,moves,pause
            
            
            
       end, end
            


            if ~isempty(moves)
                pick=fix(rand*size(moves,1))+1
                Vy(n)=Vy(n)+ moves(pick,1);
                Vx(n)=Vx(n)+ moves(pick,2);            
            end
   
    end
        %%Show
            miny=min(Vy); minx=min(Vx);
            Showy=Vy-miny+10;
            Showx=Vx-minx+10;
            fsizey=max(Showy)+10;
            fsizex=max(Showx)+10;
            field=zeros(fsizey,fsizex);
            field(sub2ind([fsizey fsizex],Showy,Showx))=100;
            %field(Showy(n),Showx(n))=300;
            image(field),pause(.1)


end