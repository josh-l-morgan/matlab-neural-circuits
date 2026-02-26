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



%% Start Reading
for i = 1:min(10,length(iNam));
    i
I = double(imread([TPN iNam{i}]));

for s = 1:20
    subI = abs(I(s+1:end,:)-I(1:end-s,:));
    dif(1,s) = mean(subI(:));
end
difs(i,:) = dif;
subplot(1,2,1)
plot(dif/max(dif)),pause(.01)
%plot((dif-dif(1))/max(dif-dif(1))),pause(.01)
comp2to(i) = (dif(2)-dif(1))/(dif(end)-dif(1));
hold on
end
hold off
ylim([0 1])
subplot(1,2,2)
hist(comp2to)

iNam{comp2to<.2}



save([TPN 'difs.mat'],'difs')