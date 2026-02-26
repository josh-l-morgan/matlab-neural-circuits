clear all

TPN = GetMyDir;
inams = GetPics(TPN);
overlap = 0.05;

inf = imfinfo([TPN inams{1}]);

for i = 1: length(inams)
   nam = inams{i};
   ms = find(nam == 'm');
   dots = find(nam == '.');
   mos(i) = str2num(nam(ms(end)+1:dots(end)-1));
end


mkdir([TPN 'maxs'])
for m = 1:max(mos)
    sprintf('running tile %d of %d',m,max(mos))
   planes = find(mos==m); 
   maxM = double(imread([TPN inams{planes(1)}]));
   for p = 2:length(planes)
       I = double(imread([TPN inams{planes(p)}]));
       maxM(:,:,2) = maxM(:,:,2) + I(:,:,2);
       maxM(:,:,1) = max(maxM(:,:,1),I(:,:,1)); 
       maxM(:,:,3) = max(maxM(:,:,3),I(:,:,3)); 
       %image(uint8(maxM)),pause
   end
   maxM(:,:,2) = maxM(:,:,2)/length(planes);
   imwrite(uint8(maxM),[TPN 'maxs\max' num2str(m) '.tif'],'Compression','none');
end

['finished ' TPN]

don = rand(1000,1)*.2;
sound(don,1000)

