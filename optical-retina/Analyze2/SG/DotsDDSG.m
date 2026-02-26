%function[] = anaDD(TPN)
% Dot Dend Does something useful with puncta and dendritic data from FD
%%and Skel please

'Starting analysis of Dots and Dendrites'


%% Get Path to Big Centroid
%Get directory name
KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:size(Kdir,1)
    clear ADotDend ADots ALength APDotDend APDots APLength AllSegCut AverageDensity
    clear BC Back BoarderFind Cname DPN DPNd DendBD DendBDmin DendBinDepth DotBD
    clear DendBDmin DendBinDepth DotBD DotBinDepth DotDendAt DotPerDendDepth Dots PeakFind
    clear Results SegLength SegMP TotalDots TotalLength VDotDend VDots VLength VPDotDend VPDots
    clear VPLength Valley a ase b bin bin2 d de maxDepth minDepth peaks r top vcent vdm vleft vr vright
    clear x xyum y yxum zum
    TPN = [KPN '\' Kdir(k).name '\'];
    DPN = [TPN 'I\']   
    DPNd=[TPN 'data\'];
    
    de(k) = exist([TPN 'find\SG.mat']);
    ase(k) = exist([DPNd 'AllSegCut.mat']);
    if de(k) & ase(k) %Get results if info available
    
    yxum=0.103;
    xyum=yxum;
    zum=0.3;

    Cname=TPN;
    %choose source data based on availability
    if exist([TPN 'Dots.mat']) & exist([TPN 'data\AllSegCut.mat']) 
        load([TPN 'Dots.mat']) 
        load([TPN 'find\SG.mat'])
        load([TPN 'data\AllSegCut.mat'])
        load([TPN 'ShiftInfo.mat']) 
        
        BC=Dots.Pos; %Extract Dot Position
        BC=BC(SG.passF,:); 
        %BC(:,3)=BC(:,3)+ShiftInfo.ShiftDots; %shift dots accordingly
        BC(:,1:2)=BC(:,1:2)* 0.103;
        BC(:,3)=BC(:,3)* 0.3;
        


    TotalDots=sum(SG.passF);
    save([DPNd 'TotalDots.mat'],'TotalDots')

    %% Assign Dend properties
    'Finding Dendrite Segment Properties'

    %%Find Segment Lengths
    SegLength=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 +(AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 +(AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);
    save([DPNd 'SegLength.mat'],'SegLength')
    TotalLength=sum(SegLength);
    save([DPNd 'TotalLength.mat'],'TotalLength')
    AverageDensity=TotalDots/TotalLength;
    save([DPNd 'AverageDensity'], 'AverageDensity')

    %%Find Segment Midpoints
    SegMP=(AllSegCut(:,:,1)+AllSegCut(:,:,2))/2;

    %% Bin with Depth
    'Binning Data with Depth'

    %%Depth
    bin=zum;%zum;
    step=zum;
    bin2=7; %size to average
    filt=ones(bin2,1);

    maxDepth=fix(max(max(SegMP(:,3)), max(BC(:,3))))+1;
    minDepth=max(fix(min(min(SegMP(:,3)), min(BC(:,3)))),1);
    minDepth = minDepth-mod(minDepth,zum); %Set min depth to factor of zum
    DotBinDepth=zeros(round((maxDepth-minDepth)/step+1),1);
    DendBinDepth=zeros(round((maxDepth-minDepth)/step+1),1);

    b=0;
    for d=minDepth:step:maxDepth
        b=b+1;
        DotBinDepth(b)=sum(BC(:,3)>=(d-bin/2) & BC(:,3) <= (d+bin/2));
        DendBinDepth(b)=sum(SegLength(SegMP(:,3)>=(d-bin/2) & SegMP(:,3) <= (d+bin/2)));
    end %end d, run all depth bins

    %%filter to smooth
    DotBD=filter(filt,bin2,DotBinDepth);
    DendBD=filter(filt,bin2,DendBinDepth);


    DotPerDendDepth=DotBD./DendBD;

    subplot(3,5,1:4)
    plot(DendBD,'r')
    subplot(3,5,6:9)
    plot(DotBD,'g')
    subplot(3,5,11:14)
    plot(DotPerDendDepth)
    pause(.01)
    %}

    %% Identify Strata
    PeakFind='Manual'
    BoarderFind='full width half maxima';
    %Have strata been previously identified
    clear P right left
    if exist([DPNd 'Results.mat']), 
        load([DPNd 'Results.mat'])
        for a=1:size(Results.Arbor,2)
            left(a)=Results.Arbor(a).Top/zum;
            right(a)=Results.Arbor(a).Bottom/zum;
            P(a,1)=Results.Arbor(a).Peak/zum
        end
    else
        %%manual Identify Strat

        if exist([DPNd 'Results.mat']), 
            load([DPNd 'Results.mat'])
            Results.Arbor.Peak
        end
        'Select peaks and click return'
        [x y]=ginput
        P=round(x);
        peaks=DendBinDepth*0;
        if size(x,1)>1
            'Identify Valley '
            [x y]=ginput;
            Valley=round(x);
        else, Valley=0; %stray value for Valley
        end

            DendBDmin=imregionalmin(DendBD);
        %%find peak edges using full width half maxima and reginal minimus
        for i =1:size(P,1)
            top=DendBD(P(i));
            left(i)=1; right(i)=size(DendBD,1);
            for l=1:P(i)-1
                if (DendBD(P(i)-l)< top/2) | sum((P(i)-l)==Valley)
                    left(i)=P(i)-l; break;   end %search for left side
            end %end searching for left side
            for r=1:size(DendBD,1)-P(i)
                if (DendBD(P(i)+r)< top/2) | sum((P(i)+r)==Valley)
                    right(i)=P(i)+r; break;   end %search for left side
            end %end searching for left side
        end %find peak width
    end % get edges of peaks

    %}

    %%find arbor stats
    APLength=DendBD*0;
    APDots=DendBD*0;
    APDotDend=DendBD*0;
    for i=1:size(P,1)
        ALength(i)=sum(DendBinDepth(left(i):right(i)));
        APLength(left(i):right(i))=ALength(i);
        ADots(i)=sum(DotBinDepth(left(i):right(i)));
        APDots(left(i):right(i))=ADots(i);
        ADotDend(i)=ADots(i)/ALength(i);
        APDotDend(left(i):right(i))=ADotDend(i);
    end %run all arbors


    %%Look at Valley

    if size(P,1)>1 %if bi
        vdm=DendBD*0;
        vr=2; %valley width = vr *2 +1
        for i=min(right):max(left) %run valley
            vdm(round(i))=mean(DotPerDendDepth(i-vr:i+vr));
        end
        vcent=find(vdm==min(vdm(min(right):max(left))));

        VPLength=DendBD*0;
        VPDots=DendBD*0;
        VPDotDend=DendBD*0;

        %%Find edges of valley
        vleft=vcent-vr;
        vright=vcent+vr;
        %%Find stats
        VLength=sum(DendBinDepth(vleft:vright));
        VPLength(vleft:vright)=VLength;
        VDots=sum(DotBinDepth(vleft:vright));
        VPDots(vleft:vright)=VDots;
        VDotDend=VDots/VLength;
        VPDotDend(vleft:vright)=VDotDend;
    end %if bi



    %%Draw stats
    subplot(3,5,5)
    plot(APLength,'r')
    subplot(3,5,10)
    plot(APDots,'g')
    subplot(3,5,15)
    plot(APDotDend,'b')
    title(Cname)
    if exist('VPDotDend','var')
        hold on
        plot(VPDotDend,'b')
        hold off
    end
    %}

    if isdir([TPN 'images'])==0, mkdir([TPN 'images']); end %create directory to store steps
    saveas(gcf,[TPN 'images/DD'],'ai') %save figure as illustrator file in images
    save([DPNd 'DepthFigBin.mat'], 'bin2') %save size of figure bin 

    pause(2)


    %% Enter Data to Structured Array
        if exist([TPN 'data\Results.mat'])
            load([TPN 'data\Results.mat'])
        end
        %Get information
        %Results.type=input('What is the Cell type?  ', 's')
        %Results.age=input('What is the Cell age? ')
        %Results.notes=input('Enter notes here ==> ', 's')


        Results.location=Cname;
        Back=find(Cname=='\');
        Results.name=Cname(Back(size(Back,2)-1)+1:Back(size(Back,2))-1)
        Results.ImageInfo.xyum=xyum;
        Results.ImageInfo.zum=zum;
        Results.CellStats.TotalLength=TotalLength;
        Results.CellStats.TotalDots=TotalDots;
        Results.CellStats.AverageDensity=AverageDensity;
        Results.Depth.DotsBinDepth=DotBD;
        Results.Depth.DendBinDepth=DendBD;
        Results.Depth.DotPerDendDepth=DotPerDendDepth;
        Results.Depth.AverageingBinWidthInMicrons=bin2*zum;
        for i=1:size(P,1)
            Results.Arbor(i).PeakFind=PeakFind;
            Results.Arbor(i).BoarderFind=BoarderFind;
            Results.Arbor(i).Length=ALength(i);
            Results.Arbor(i).Dots=ADots(i);
            Results.Arbor(i).DotDend=ADotDend(i);
            Results.Arbor(i).Top=left(i)*zum;
            Results.Arbor(i).Peak=P(i)*zum;
            Results.Arbor(i).Bottom=right(i)*zum;
        end

        %save Results
        save([TPN 'data\Results.mat'],'Results')


    %{
    if exist('./Dat.mat')
        load('./Dat.mat')
        for i= 1:size(Dat,2)
            if strcmp(Cname,Dat(i).name); targ=i; break
            else targ=size(Dat,2)+1; end
        end
        Dat(targ)=Results;

        save('./Dat.mat','Dat')

    end %if Dat exists
    %}



    %% Draw Dot and Dend
    %{
    'Drawing Dots and Dendrites'

    Sc=(1/xyum)/2;
    DD=uint8(zeros(round(max(max(AllSegCut(:,1,:)))*Sc),round(max(max(AllSegCut(:,2,:)))),round(max(max(AllSegCut(:,3,:)))*Sc)));


    %%Draw Segments
    SkelRes=.1;
    for i=1:size(AllSegCut,1)
        Dist=sqrt((AllSegCut(i,1,1)-AllSegCut(i,1,2))^2 + (AllSegCut(i,2,1)-AllSegCut(i,2,2))^2 + (AllSegCut(i,3,1)-AllSegCut(i,3,2))^2); %find distance
        Length(i)=Dist;
          devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
        for d=1:devs+1
            sy=AllSegCut(i,1,1)+((AllSegCut(i,1,2)-AllSegCut(i,1,1))/devs)*(d-1);
            sx=AllSegCut(i,2,1)+((AllSegCut(i,2,2)-AllSegCut(i,2,1))/devs)*(d-1);
            sz=AllSegCut(i,3,1)+((AllSegCut(i,3,2)-AllSegCut(i,3,1))/devs)*(d-1);
            DD(round(sy*Sc)+1,round(sx*Sc)+1,round(sz*Sc)+1)=1; %draw Skel
        end
    end
    clear Dist

    %%DrawNodes
    for i=1:size(AllSegCut,1)
        DD(round(AllSegCut(i,1,1)*Sc)+1,round(AllSegCut(i,2,1)*Sc)+1,round(AllSegCut(i,3,1)*Sc)+1)=2;
        DD(round(AllSegCut(i,1,2)*Sc)+1,round(AllSegCut(i,2,2)*Sc)+1,round(AllSegCut(i,3,2)*Sc)+1)=2;
    end

    %%Draw Dots
    for i=1:size(BC,1)
        DD(round(BC(i,1)*Sc)+1,round(BC(i,2)*Sc)+1,round(BC(i,3)*Sc)+1)=3;
    end

    save([TPN 'pics\DD.mat'])

    %}

    %% Finish

    [TPN(size(TPN,2)-6:size(TPN,2)-1)]
    DotDendAt=uint16(clock)
    save([TPN 'data/DotDendAt.mat'],'DotDendAt')
    else
        'Data Not Available'
    end
    'Done DotDend'
    end

end % run all cells
