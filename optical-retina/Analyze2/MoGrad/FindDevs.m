function[Use] = FindDevs(Step)
%%find gaussian SDs the correspond to series of inner/outer ratios
%%Step = incriment between Rdend ratios

fsize=100;

Rad=fsize/2; % Maximum radius of arbor



%% get dists
Dmap=zeros(fsize,fsize);
    for y = 1:fsize
        for x= 1 : fsize
            Dmap(y,x) = sqrt((y-fsize/2)^2 + (x-fsize/2)^2);      
        end
    end
%image(Dmap(:,:)),pause(.04) 

%% Run all StanDevs
StanDevR=.4:.01:2;
StanDev=StanDevR*Rad; %1 Stadard deviation of dendritic gaussian
clear Rdend
for i = 1:size(StanDev,2)
    
    
%% f(x) = ae^(-(x-b)^2/(2c^2))
a = 1;
b = 0;
c = StanDev(i);
e = 2.71828183;

showg=a * e.^(-1*((0:Rad)-b).^2./(2*c^2));
%plot(showg),pause(.4)

Gmap=a * e.^(-1*(Dmap-b).^2./(2*c^2));
%image(Gmap*100)

InG=mean(Gmap(Dmap<(Rad/2)));
OutG=mean(Gmap(Dmap>(Rad/2) & Dmap<Rad));


Rdend(i)=InG/OutG;

end
scatter(StanDevR,Rdend)

c=1;
UseDend=Rdend(1);
Use=StanDevR(1);
for i = 2: size(StanDevR,2)
    Test=UseDend(c);
    if Rdend(i) < Test- Step
    c=c+1;
        Use(c)=StanDevR(i);
        UseDend(c)=Rdend(i);
    end
       
end

Use