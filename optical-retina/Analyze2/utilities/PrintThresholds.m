TPN=GetMyDir

Tdir=dir(TPN); Tdir=Tdir(3:size(Tdir,1));
ListThresh=cell(size(Tdir,1),2);
for i = 1: size(Tdir,1)
   GetThresh=[TPN '\' Tdir(i).name '\data\Threshold.mat'];
   ListThresh(i,1)=cellstr(Tdir(i).name);
   if exist(GetThresh)
       load(GetThresh)
       ListThresh(i,2)=num2cell(Thresh);
      
   else
       ListThresh(i,2)=cellstr('None Found');
   end
end