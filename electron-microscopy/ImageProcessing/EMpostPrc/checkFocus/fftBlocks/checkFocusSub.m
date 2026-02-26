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

myCol = {'r','b'};

%% Start Reading
for i = 1:2%:length(iNam);
    i
    iDiv = 10;
    I = double(imread([TPN iNam{i}]));
    [ys xs] = size(I);
    xI = fix(1:xs/iDiv:xs-xs/iDiv)+1;
    xI(xI>xs) = xs;
    yI = fix(1:ys/iDiv:ys-ys/iDiv)+1;
    yI(yI>ys) = ys;
    for y = 1:length(yI)-1
        for x = 1:length(xI)-1
            %subplot(1,2,1)
            Isub = I(yI(y):yI(y+1) ,xI(x):xI(x+1));
           % image(Isub),pause(.01)
            for s = 1:20
                subI = abs(Isub(s+1:end,:)-Isub(1:end-s,:));
                dif(1,s) = mean(subI(:));
            end
            difs(i,:) = dif;
            %subplot(1,2,2)
            %plot(dif),pause(0.01)
            %plot(dif/max(dif)),pause(.01)
            %plot((dif-dif(1))/max(dif-dif(1))),pause(.01)
            comp2to(y,x) = (dif(2)-dif(1))/(dif(3)-dif(1));
            comp2to(y,x) = dif(1)/dif(5);
            hold on
        end
      

    end
    
    
  subplot(2,2,2*(i-1) + 1)
  image(I)
  subplot(2,2,2*(i-1) + 2)
  image(comp2to*200)
    
end

hold off
  
  
        %subplot(1,2,2)
        %hist(comp2to(:,:,i))

save([TPN 'difs.mat'],'difs')