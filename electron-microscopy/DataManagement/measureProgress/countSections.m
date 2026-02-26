%% Speed/Productivity check


SPN = GetMyDir

dSPN = dir(SPN), dSPN = dSPN(3:end);
%%
sNum = 0;
for i = 1:length(dSPN)
    
    nam = dSPN(i).name;
    if sum(regexp(nam,'_Montage'))
        sNum = sNum + 1;
        secInfo(sNum) = dSPN(i);
    end
    
end

secDirDates = cat(1,secInfo.datenum);
secDirDates = secDirDates - min(secDirDates);
histSec = hist(secDirDates,0:1:max(secDirDates));

%% 
for i = 1:length(secInfo)
    sprintf('Checking folder %d of %d',i,length(secInfo))
   nam = secInfo(i).name;
   dSec = dir([SPN nam]);
   fdates(i) = min(cat(1,dSec.datenum));
   
    
end
    
%%
ddates = (fdates-min(fdates));
hdates = hist(ddates,0:1/8:max(ddates));
bar(hdates)

TBperday = hdates * 10/1000*8;
bar(TBperday)
% 
% maxTB = 20000000/1000000000000 * 60 * 60 * 24;
% 
% capacity = TBperday/maxTB;
% bar(capacity)

capacity = TBperday

%%Running average
weekCap = capacity * 0;
for i = 1:length(capacity)
    startC = max(1,i-(15*8));
    stopC = min(length(capacity),i+(15*8));
    weekCap(i) = mean(capacity(startC:stopC));
    
end
hold on
plot(weekCap,'r','Linewidth',5)
xlabel(datestr(min(fdates)))
hold off






    