


%{
%input 1=time(h) 2=length(um) 3=startNum 4=gain 5=loss
c = [] ;% enter control data
e = []; % enter experimental data
%}

cDD=c(:,3)./c(:,2);
eDD = e(:,3)./e(:,2);

cGain = c(:,4)./c(:,2)./c(:,1) *1000;
eGain = e(:,4)./e(:,2)./e(:,1) *1000;

cLoss = c(:,5)./c(:,3)./c(:,1) *1000;
eLoss = e(:,5)./e(:,3)./e(:,1) *1000;


cDat=[cDD  cGain cLoss];
eDat=[eDD  eGain eLoss];

mean(cLoss)

%% Run from sums
cSum=sum(c,1);
eSum=sum(e,1);

cSumDD = cSum(3)./cSum(2);
eSumDD = eSum(3)./eSum(2);

cSumGain = cSum(4)./cSum(2)./cSum(1) *1000;
eSumGain = eSum(4)./eSum(2)./eSum(1) *1000;

cSumLoss = cSum(5)./cSum(3)./cSum(1) *1000;
eSumLoss = eSum(5)./eSum(3)./eSum(1) *1000;

cDatS= [cSumDD  cSumGain  cSumLoss];
eDatS= [eSumDD  eSumGain  eSumLoss];

%% Plot
cN=size(cDat,1);
eN=size(eDat,1);

for i = 1: 3
    
    subplot(1,3,i)
    hold off
    cE=std(cDat(:,i))/sqrt(cN);
    eE=std(eDat(:,i))/sqrt(eN);
    eM=mean(eDat(:,i));
    cM=mean(cDat(:,i));
    eS=eDatS(i);
    cS=cDatS(i);

    errorbar([cM eM],[cE eE],'k','LineWidth',2,...
        'LineStyle','none','Marker','x','MarkerSize',20)
    hold on
    scatter(ones(cN,1),cDat(:,i),'b')
    scatter(ones(eN,1)*2,eDat(:,i),'b')

    %scatter(1.1,cS,'r')
    %scatter(2.1,eS,'r')
    
    y=ylim;
    ylim([0 y(2)])

end














