
KPN=GetMyDir;

Kdir=dir(KPN); Kdir=Kdir(3:size(Kdir,1));


for k = 1:size(Kdir,1) %run all files in directory

    TPN = [ KPN Kdir(k).name] %get image name
    if strcmp('tif',TPN(size(TPN,2)-2:size(TPN,2))) %if a tif
        I=imread([TPN]); %read image
        image(I),pause(.01) %Show raw
        Im=max(I,[],3); %combind channels
        It=Im>median(Im(:)); %Threshold out background
        Il=bwlabel(It); %find all dots

        Dnum=max(Il(:));

        %% find centers

        for d= 1:Dnum
            [y x] = find(Il==d);
            ym=mean(y); xm=mean(x);
            DPos(d,:)=[ym xm];
        end

        Dists=zeros(Dnum);
        for d = 1:Dnum
            Dists(d,:)=sort(dist(DPos, DPos(d,:)));
        end

        mDists=mean(Dists,1);
        Near=Dists(:,2);
        hNear=hist(Near,2.5:5:max(Near(:)));
        
        Name=Kdir(k).name(1:size(Kdir(k).name,2)-4);
        xlswrite([KPN Name '.xls'],mDists','meanDists')
        xlswrite([KPN Name '.xls'],Near,'NearestNeighbor')
        

    end %if a tif
end %run all