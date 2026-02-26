
%% Collect and format data

clear Name Data dat DatTime
[DFN DPN] = uigetfile('.xls','Get Time Point 1','C:\Documents and Settings\a\My Documents\work\Huckfeldt\Fixed\')
Data(:,:,1)=xlsread([DPN DFN])
Name(1,1)={[DPN DFN]};
[DFN DPN] = uigetfile('.xls','Get Time Point 2','C:\Documents and Settings\a\My Documents\work\Huckfeldt\Fixed\')
Data(:,:,2)=xlsread([DPN DFN])
Name(2,1)={[DPN DFN]};
[DFN DPN] = uigetfile('.xls','Get Time Point 3','C:\Documents and Settings\a\My Documents\work\Huckfeldt\Fixed\')
Data(:,:,3)=xlsread([DPN DFN])
Name(3,1)={[DPN DFN]};

DataMax=max(Data,[],3);
GoodDat=DataMax>0 & DataMax<1;

for i = 1: 3
    dat=Data(:,:,i);
    DatTime(:,i)=dat(GoodDat>0);
end



slash=find(DPN == '\');
GPN=DPN(1:slash(size(slash,2)-1));

xlswrite([GPN 'All' DFN],DatTime)