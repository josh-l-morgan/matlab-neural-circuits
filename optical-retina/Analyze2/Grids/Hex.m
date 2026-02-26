fsize=100;
Diam=50;

field=zeros(fsize);

Grid=hexes([fsize fsize],Diam);

Col=rand(size(Grid,1),1)*255;


for y = 1:fsize
    for x=1:fsize
        Pos=[y x];
        Dist=dist(Grid,Pos);
        Pick=find(Dist==min(Dist),1);
        field(y,x)=Col(Pick);
    end
end

image(field),pause(.01)