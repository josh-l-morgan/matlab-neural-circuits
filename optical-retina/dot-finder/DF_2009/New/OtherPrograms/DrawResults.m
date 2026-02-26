TPN = GetMyDir;
TPNd = dir(TPN); TPNd=TPNd(3:length(TPNd));

for i = 1:length(TPNd)
   nam = TPNd(i).name
   if strcmp(nam(length(nam)-3:length(nam)),'.mat')
        load([TPN nam])
   end
end


Hit= Classify.Hit;
Miss = Classify.Miss;
Unknown = 1:length(Classify.only);
Unknown = setdiff(Unknown,Hit);
Unknown = setdiff(Unknown,Miss);
Found = Unknown(result>0);
Eliminated=Unknown(result==0);

ImSize = Dots.ImSize;

I = zeros(ImSize,'uint8');

%% Draw Hits
for i = 1:length(Hit)
   I(Dots.Vox(Hit(i)).Ind) = 255;   
end

%% Draw Found
for i = 1:length(Found)
   I(Dots.Vox(Found(i)).Ind) = 100;       
end

imwriteNp(TPN,I,'Found');

I = I * 0;

%% Draw Hits
for i = 1:length(Miss)
   I(Dots.Vox(Miss(i)).Ind) = 255;   
end

%% Draw Found
for i = 1:length(Eliminated)
   I(Dots.Vox(Eliminated(i)).Ind) = 100;       
end

imwriteNp(TPN,I,'Eliminated')


%% Draw Missed Hits
for i = 1:length(Hit)
   I(Dots.Vox(Hit(i)).Ind) = 100;       
end

plot(result)




