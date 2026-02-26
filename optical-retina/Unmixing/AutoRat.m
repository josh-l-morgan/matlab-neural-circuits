 
colormap gray(256)

med1=median(channel1(:));
med2=median(channel2(:));
std1=uint16(std(single(channel1(:))));
std2=uint16(std(single(channel2(:))));

Thresh1=med1+3*std1;
Thresh2=med2+3*std2

Cht1=channel1>Thresh1;
Cht2=channel2>Thresh2;
Cht12=Cht1 | Cht2;
Cht12 = Cht12 & (channel1 < max(channel1(:))-max(channel1(:))/10) & (channel2 < max(channel2(:))-max(channel2(:))/10);

image(sum(Cht1,3)*20);
image(sum(Cht2,3)*20);
image(sum(Cht12,3)*30);

BackC1=mean(channel1((channel1<(med1+std1))&(channel2<(med2+std2))));
BackC2=mean(channel2((channel1<(med1+std1))&(channel2<(med2+std2))));

Dat=[channel1(Cht12) channel2(Cht12)];

Rat=single(channel1)./(single(channel1)+single(channel2));
[IDX,C] = kmeans(Rat,2);

ChxSi1 = Dat(IDX==1,1);
ChySi1 = Dat(IDX==1,2);

ChxSi2 = Dat(IDX==2,1);
ChySi2 = Dat(IDX==2,2);

b1 = robustfit(ChxSi1,ChySi1);
b2 = robustfit(ChxSi2,ChySi2);

% jlm1 = median(ChySi1)/median(ChxSi1);
% jlm2 = median(ChySi2)/median(ChxSi2);

Sub=1:length(channel1(:));
Sub=mod(Sub,300)==0;
scatter(channel1(Cht12(:) & Sub'),channel2(Cht12(:) & Sub'))

