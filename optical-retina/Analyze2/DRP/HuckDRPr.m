
KPN=GetMyDir;

Kdir=dir(KPN); Kdir=Kdir(3:size(Kdir,1));

colormap gray(255)

for k = 1:size(Kdir,1) %run all files in directory

    DPN = [ KPN Kdir(k).name] %get image name
    if strcmp('tif',DPN(size(DPN,2)-2:size(DPN,2))) %if a tif

        TPN=[DPN(1:size(DPN,2)-4) '\']
        mkdir(TPN)

        I=imread([DPN]); %read image
        image(I),pause(.01) %Show raw
        Im=max(I,[],3); %combind channels
        It=Im>median(Im(:)); %Threshold out background
        Il=bwlabel(It); %find all dots

        Dnum=max(Il(:));

        %% find centers
        DPos=zeros(Dnum,2);
        for d= 1:Dnum
            [y x] = find(Il==d);
            ym=mean(y); xm=mean(x);
            DPos(d,:)=[ym xm];
        end
        %% Find Dists to centers
        Dists=zeros(Dnum);
        for d = 1:Dnum
            Dists(d,:)=sort(dist(DPos, DPos(d,:)));
        end

        mDists=mean(Dists,1);
        Near=Dists(:,2);
        hNear=hist(Near,2.5:5:max(Near(:)));

        Name=Kdir(k).name(1:size(Kdir(k).name,2)-4);
        xlswrite([TPN Name '.xls'],mDists','meanDists')
        xlswrite([TPN Name '.xls'],Near,'NearestNeighbor')



        %% Density Recovery Profile
        %Calculate distance to areas
        A=I*0+1;
        [Ay Ax]=find(A);
        APos=[Ay Ax];
        clear Ay Ax
        maxDist=max(Dists(:));
        DUsed=((0:fix(maxDist))+.5)';
        Cnum=zeros(Dnum,size(DUsed,1),1);
        Anum=Cnum;
        
        for c = 1:Dnum
            clear DistsA
            DistsA=dist(APos, DPos(c,:));
            %% Map Density by distance Area
            Cnum(c,:)=hist(Dists(c,:),DUsed);
            Anum(c,:)=hist(DistsA,DUsed);
            PercentDoneMapping=c/Dnum*100
        end
        Csum=sum(Cnum,1)';
        Asum=sum(Anum,1)';
        CoA=Csum./Asum;
        plot(CoA),pause(.01)
        bar(Asum)
        bar(Csum)
        bar(CoA)
        bar(CoA(2:200))
        %ylim([-max(CoA(:))/10 max(CoA(:))])

%         save([TPN 'Cnum.mat'],'Cnum')
%         save([TPN 'Anum.mat'],'Anum')
%         xlswrite([TPN Name '.xls'],Cnum,'Cnum')
%         xlswrite([TPN Name '.xls'],Anum,'Anum')

        %% Bin
        Bin=6;
        for d = 1 : size(Asum,1)
            Cbin(d)=sum(Csum(max(d-Bin/2,2):min(d+Bin/2,size(Csum,1))));
            Abin(d)=sum(Asum(max(d-Bin/2,2):min(d+Bin/2,size(Asum,1))));
        end
        CAbin=Cbin./Abin;
        plot(CAbin)
        bar(Abin)
        bar(Cbin)
        subplot(2,1,1)
        image(It*1000)
        subplot(2,1,2)
        bar(CAbin(1:200)),pause(.01)
%         xlswrite([TPN Name '.xls'],CAbin,['Density Binned ' num2str(Bin)])
        %%

    end %if a tif
end %run all










