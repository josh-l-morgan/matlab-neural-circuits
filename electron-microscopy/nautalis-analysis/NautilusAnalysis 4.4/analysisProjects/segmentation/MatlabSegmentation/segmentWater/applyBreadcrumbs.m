
clear all
TPN = GetMyDir;  %Get Tif

labName = 'fused2\'; 
load([TPN 'binCrumb.mat'])

obs = binCrumb.obs;
clear binCrumb

%% Transform crumbs
% 150 - z, x * 2, y * 2

mins = [ 10 ^10 10^10 10^10];
maxs = [1 1 1];
for o = 1:length(obs)
   ob = obs{o};
   ob(:,1) = 151 - ob(:,1);
   ob(:,2:3) = ob(:,2:3) * 2;
   obs{o} = ob;
   
   mins = min(min(ob,[],1),mins);
   maxs = max(max(ob,[],1),maxs);
end

%% turn into planes

planes = cell(maxs(1),1);

for o = 1:length(obs)
   ob = obs{o};
   for i = 1:size(ob,1)
      plane = planes{ob(i,1)}; 
      plane = [plane; o ob(i,2:3)];
      planes{ob(i,1)} = plane;
   end
end


%% Grab cell
grab = 114;

lnams = getPics([TPN labName]);
inams = getPics([TPN]);
rasOb = zeros(2048,2048,160,'uint16');
myCol = hsv(256);
myCol(1,:) = 0;
colormap(myCol)

for i = 1:length(planes);
    
    sprintf('Running plane %d of %d',i,length(planes))
    %%Get data
    ras2D = rasOb(:,:,i);
    lI = imread([TPN labName lnams{i}]);
    plane = planes{i};
    planeObs = unique(plane(:,1));
    for g = 1: length(planeObs)
        pos = plane(plane(:,1) == g,:);
        posi = sub2ind(size(lI),pos(:,3),pos(:,2));
    %     ras2D(posi) = 1000;
    %     ras2D = imdilate(ras2D,strel('disk',10));
        ids = lI(posi);
        ids = unique(ids(ids>0));
        for l = 1:length(ids)
           ras2D(lI == ids(l)) = g; 
        end
    end
    rasOb(:,:,i) = ras2D;
    
    image(mod(ras2D,256)),pause(.01)
    
end
    
    
    image(sum(rasOb,3))
    



%% Display with Label

lnams = getPics([TPN labName]);
inams = getPics([TPN]);
for i = 1;%:length(planes);
    %%Get data
    lI = imread([TPN labName lnams{i}]);
    %I = imread([TPN inams{i}]);
    I3 = cat(3,I,I,I);
    plane = planes{i};
    
    myCol = hsv(double(max(lI(:)))+1)* 100;
    [r rix] = sort(rand(size(myCol,1)),1);
    myCol= myCol(rix,:);
    myCol = cat(1,[0 0 0],myCol);
    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
    skipCol = 5 + round(rand*10);

    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
    image(I3  + colO),pause(.01)
    xax = get(gca,'xlim'); yax = get(gca,'ylim');
    
    hold on
    scatter(plane(:,2),plane(:,3),'.','LineWidth',10)
    hold off
    
end
    
    
    
    
    
    


