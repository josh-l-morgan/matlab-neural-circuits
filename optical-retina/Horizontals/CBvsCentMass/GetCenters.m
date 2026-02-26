%%Retrieves the centers of masses from individually stored tifs 
%%and relates them to reference point

TPN = GetMyDir;

TPNd=dir(TPN);TPNd=TPNd(3:length(TPNd));
Names={};
for i = 1:length(TPNd)
    nam=TPNd(i).name;
    ln=length(nam);
    if nam(ln-3:ln)=='.tif'
        Names{length(Names)+1}=nam;
    end
end

Cent=zeros(length(Names),2);
clear Ia
for i = 1: length(Names)
    Is=imread([TPN char(Names{i})]);
    if ~exist('Ia'),Ia=Is*0; end
    Ia(Is>0)=Ia(Is>0)+1;
    [y x z] = size(Is);
    siz= [y x];
    [y x]=find(Is>0);
    Cent(i,:)=[mean(y) mean(x)];
    image(Is),pause(.01)
end

Vals=Ia(Ia>0);
[y x] = find(Ia);
Scaled= [y x] .* double([Vals Vals]);
CentAll=sum(Scaled,1)/sum(Vals);


Reference=[463.5809476	466.4295663];


Dists=dist(Cent,Reference);


image(Ia*100),pause(.01)
