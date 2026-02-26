clear all
colormap gray(256)
TPN = GetMyDir;
dTPN = dir(TPN);
dTPN = dTPN(3:end);

iNam = [];
for i = 1:length(dTPN)
    nam = dTPN(i).name; 
    if length(nam)>4
        if strcmp(nam(end-3:end), '.tif' )
        iNam{length(iNam)+1} = nam;
        end
    end
end

for i = 1 : length(iNam)
   rawI{i} = imread([TPN iNam{i}]); 
   I(i,:,:) = double(rawI{i});
end

%% Scale 
backGround = mean(mean(I(:,1:3,:),1),2);
sI = I;
sI(:,:,1) = sI(:,:,1)-  backGround(1);
sI(:,:,2) = sI(:,:,2) - backGround(2);
for i = 1:size(sI,1)
   for c = 1:2
       sI(i,:,c)=sI(i,:,c)/mean(sI(i,:,c));
   end
end
image(uint8(sI * 50)),pause(.1)

%%consolidate
msI = squeeze(mean(sI,1));
image(uint8(msI * 50))
bar(msI)
colormap([1 0 0; 0 1 0])

