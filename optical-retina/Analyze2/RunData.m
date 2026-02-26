%% Run all data in a folder

DPN=uigetdir; %get directory containing all data folders
DPN=[DPN '\'];
d=dir(DPN); %make directory listing
d=d(3:size(d,1));
c=0;

Report=cell(1,4);
Report(1,1)={'Folder'};
Report(1,2)={'D Status'};
Report(1,3)={'Dot Status'};
Report(1,4)={'Dend Status'};
c=1;

for i=1:size(d,1), 
    if d(i).isdir
       c=c+1;
       Report(c,1)={d(i).name};
     if exist([DPN d(i).name '\data'])
        TPN=[DPN  d(i).name];
        Result=anaFix(TPN)
        Report(c,2:4)=Result(1,2:4);
     end %is there a data file?
        
    end %is it a directory
    progress=i/size(d,1)*100
end %run all files in DPN


