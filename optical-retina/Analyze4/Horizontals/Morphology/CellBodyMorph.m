%%Retrieves the centers of masses from individually stored tifs
%%and relates them to reference point
clear all


%% Remember
TPN = GetMyDir
TPNi=[TPN 'I\'];
xyum=0.153;
xyum = 0.23;

Idir=dir(TPNi); Idir=Idir(3:length(Idir));

Names={};
for i = 1:length(Idir)
    nam=Idir(i).name;
    LN = length(nam);
    if nam(LN-3:LN)=='.tif'
        Names{length(Names)+1}=nam;
    end
end

clear I
for i = 1:length(Names)
    I(:,:,i)=imread([TPNi Names{i}]);
end


Cells=unique(I(:));
Cells=Cells(Cells>0);
[ys xs zs]=size(I);
siz=[ys xs zs];
numCells=length(Cells);

%% Get Cell info


for i = 1: numCells %% Make 2D pics
    c=Cells(i);
   Ic(:,:,i)=sum(I==c,3); 
end
image(sum(Ic,3)*30),pause(.01)

%%Get Area
for i = 1: numCells
    Area(i)=length(find(Ic(:,:,i)>0));
end
Area=Area*xyum^2;

%%Get Equivilant Radius A=pi R^2
ER=sqrt(Area/pi);
Stats(ER)


%{

[I numCells] = bwlabeln(I);
for i = 1: numCells
    Area(i)=length(find(I==i));
end
Area=Area*xyum^2;

%%Get Equivilant Radius A=pi R^2
ER=sqrt(Area/pi);
Stats(ER)



%}












