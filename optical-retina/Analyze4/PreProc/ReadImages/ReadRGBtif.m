%% read in RGB tiffs

TPN=GetMyDir;

TPNd=dir(TPN); TPNd=TPNd(3:length(TPNd));
zs=length(TPNd);

It=imread([TPN TPNd(1).name]);
Class=class(It);
[ys xs cs]= size(It);
I = zeros(ys, xs, cs, zs, Class);
I(:,:,:,1)=It; clear It;

for i =2:zs
    I(:,:,:,i)=imread([TPN TPNd(i).name]);
end

