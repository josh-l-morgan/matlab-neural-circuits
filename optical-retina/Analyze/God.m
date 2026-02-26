%% Look at Data from all Cells while working out power issues

clear all

%% Read Data from Raid into Master List
figure
%%Map Raid folders
Raid='\\128.208.64.36\wonglab\Josh\Analyzed\Output';
RaidDir=dir(Raid);
RD=struct2cell(RaidDir);
RD=RD(1,:); RD=RD(3:size(RD,2));

%%Compile Cell, Status, and Results into master list
for i=1:size(RD,2)
    if exist([Raid '\' char(RD(i)) '\Cell.mat'])
    load([Raid '\' char(RD(i)) '\Cell.mat'])
    Master(i).Cell=Cell; end
    if exist([Raid '\' char(RD(i)) '\Status.mat'])
     load([Raid '\' char(RD(i)) '\Status.mat'])
    Master(i).Status=Status; end
    if exist([Raid '\' char(RD(i)) '\data\Results.mat'])
    load([Raid '\' char(RD(i)) '\data\Results.mat'])
    Master(i).Results=Results; end
end



%% Remove Cells that are not Ready

%%pick out cells
for i=1:size(Master,2), c(i)=~isempty(Master(i).Cell); end
Master=Master(c); clear c

for i = 1 : size(Master,2)
    AC=struct2cell(Master(i).Cell)';
    AS=Master(i).Status.Results;
    AllCells(i,1)=AC(1);
    AllCells(i,2)=AC(3);
    AllCells(i,3)=AC(4);
    AllCells(i,4)={AS};
    AllCells(i,5)=AC(2);
end



%%pick out results
for i=1:size(Master,2), c(i)=~isempty(Master(i).Results); end
Master=Master(c); clear c



%% Make lists of Cell types and Ages

ONs=logical(zeros(1,size(Master,2)));
OFFs=ONs; BIs=ONs; Others=ONs; Unknowns=ONs;

for i = 1:size(Master,2)

    %% read Cell type and mark appropriate matrix
    switch Master(i).Cell.Type
        case 'ON'
            ONs(i)=1;
        case 'OFF'
            OFFs(i)=1;
        case 'BI'
            BIs(i)=1;
        case 'Other'
            Others(i)=1;
        otherwise
            Unknowns(i)=1;
    end %End Cell type switch
    
    %% make list of Cell Class (Chalupa)
    Class(i)=str2double(Master(i).Cell.Class);
    

    %% Put cell ages in Ages list
    Ages(i)=str2double(Master(i).Cell.Age);
    %% mark unknons as zero
    if isnan(Ages(i)), Ages(i)=0; end;

end

%% Create Goupings
G(1,1,:)=ONs; G(2,1,:)=OFFs; G(3,1,:)=BIs; G(4,1,:)=Others;
G(5,1,:)=Unknowns;
G(1,2,:)=Ages==5; G(2,2,:)=Ages==7; G(3,2,:)=Ages==12; G(4,2,:)=Ages>30;
Gname(:,1)={'ONs'; 'OFFs'; 'BIs'; 'Others'; 'Unknowns'};
Gname(:,2)={'5';'7';'12';'30+';'none'};


%%Extract Data
clear AverageDensity Arbor2DD Arbor21DD
for i = 1: size(Master,2)
   AverageDensity(i)=Master(i).Results.CellStats.AverageDensity;
   Arbor1DD(i)=Master(i).Results.Arbor(1).DotDend;
   if size(Master(i).Results.Arbor,2)>1
       Arbor2DD(i) = Master(i).Results.Arbor(2).DotDend;
       MeanDD(i)=mean([Arbor1DD(i) Arbor2DD(i)]);
   else
       Arbor2DD(i)=0;
       MeanDD(i)=Arbor1DD(i);
   end
   if isfield(Master(i).Results.CellStats,'NN')
       nn=Master(i).Results.CellStats.NN;
   else,    nn(i)=0 ;   end
   NN{i}=nn;
   NNm(i)=mean(nn);
   NNu2(i)=sum(nn<=2)/sum(nn>0);
end




%% Display Results


%%Group By Age
subplot(4,1,1)
[n x]= hist(MeanDD(Ages==5),0:.05:.4); %% Show All P5
bar(x,n),title('P5','FontSize',12);
subplot(4,1,2)
[n x]= hist(MeanDD(Ages==7),0:.05:.4);
bar(x,n),title('P7','FontSize',12);
subplot(4,1,3)
[n x]= hist(MeanDD(Ages==12),0:.05:.4);
bar(x,n),title('P12','FontSize',12);
subplot(4,1,4)
[n x]= hist(MeanDD(Ages>30),0:.05:.4);
bar(x,n),title('One Month +','FontSize',12);


%%Compair ON and OFF arbors of BI
subplot(1,1,1)
ONvsOFF=Arbor1DD(BIs)-Arbor2DD(BIs);
[n x] = hist(ONvsOFF)
bar(x,n),title('Upper-Lower arbor Dot Density in Bistratified cells')

%% Show all Groupings
for t=1:4 % run all types
    for a=1:4 %run all ages
        TypeC=G(t,1,:); %Get type
        AgeC=G(a,2,:); %Get age
        subplot(4,4,(t-1)*4+a)
       [n x]= hist(MeanDD(TypeC & AgeC),0:.05:.4); %% Show All P5
       bar(x,n),title([Gname(a,2) Gname(t,1)],'FontSize',8);
       %%Calculate Ns
       Ns((t-1)*4+a,1:3)=[Gname(a,2) Gname(t,1) sum(TypeC & AgeC)];
    end
end



%% Test for Significance

[p h]=ranksum(MeanDD(Ages==5),MeanDD(Ages==7));
ProbThat_P5isP7=p

[p h]=ranksum(MeanDD(Ages==7),MeanDD(Ages==12));
ProbThat_P7isP12=p

[p h]=ranksum(MeanDD(Ages==12),MeanDD(Ages>30));
ProbThat_P12isP30=p

[p h]=ranksum(MeanDD(Ages==7),MeanDD(Ages>7));
ProbThat_P7isP12_P30=p

%%Probability that ON and OFF arbors have different densities
ONvsOFF=Arbor1DD(BIs)-Arbor2DD(BIs);
[h p]=ttest(ONvsOFF);
ProbONvsOFF=p



%% plot ages together
figure
bin= .05

%subplot(4,1,1)
[n x]= hist(MeanDD(Ages==5),0:bin:.4); %% Show All P5
plot(x,n,'r'),title('P5','FontSize',12);
hold on
%subplot(4,1,2)
[n x]= hist(MeanDD(Ages==7),0:bin:.4);
plot(x,n,'m'),title('P7','FontSize',12);
%subplot(4,1,3)
[n x]= hist(MeanDD(Ages==12),0:bin:.4);
plot(x,n,'b'),title('P12','FontSize',12);
%subplot(4,1,4)
[n x]= hist(MeanDD(Ages>30),0:bin:.4);
plot(x,n,'g'),title('One Month +','FontSize',12);
hold off



%% Bar ages together
figure
st=.05; % Start point
sp=.35; % Stop Point
bin = .1; clear n
%subplot(4,1,1)
[n(:,1) x]= hist(MeanDD(Ages==5),st:bin:sp); %% Show All P5
[n(:,2) x]= hist(MeanDD(Ages==7),st:bin:sp);
[n(:,3) x]= hist(MeanDD(Ages==12),st:bin:sp);
[n(:,4) x]= hist(MeanDD(Ages>30),st:bin:sp);
%subplot(2,1,1)
bar(x,n,1.5), title('Number of cells in .1DD bins')

np(:,1)= n(:,1)/sum(n(:,1))*100;
np(:,2)= n(:,2)/sum(n(:,2))*100;
np(:,3)= n(:,3)/sum(n(:,3))*100;
np(:,4)= n(:,4)/sum(n(:,4))*100;
%subplot(2,1,2)
bar(x,np,1.5), title('Percent of cells in .1DD bins')
legend(['P5 mean density   = ' num2str(mean(MeanDD(Ages==5)))],...
       ['P7 mean density   = ' num2str(mean(MeanDD(Ages==7)))],...
       ['P12 mean density  = ' num2str(mean(MeanDD(Ages==12)))],...
       ['P30+ mean density = ' num2str(mean(MeanDD(Ages>30)))])

%% Plot increase with age
figure

scatter(Ages,MeanDD)
hold on

%%Find standard Error as StandardDeviation / sqrt(n)
E(1)=std(MeanDD(Ages==5))/sqrt(sum(Ages==5));
E(2)=std(MeanDD(Ages==7))/sqrt(sum(Ages==7));
E(3)=std(MeanDD(Ages==12))/sqrt(sum(Ages==12));
E(4)=std(MeanDD(Ages>30))/sqrt(sum(Ages>30));


errorbar([5 7 12 35],[mean(MeanDD(Ages==5)) mean(MeanDD(Ages==7)) ...
    mean(MeanDD(Ages==12)) mean(MeanDD(Ages>30))],E,'r')

hold off


%% Plot density with age and type
figure
ONtime=[mean(MeanDD(Ages==5)) mean(MeanDD(Ages==7 & ONs)) ...
        mean(MeanDD(Ages==12 & ONs)) mean(MeanDD(Ages>30 & ONs))];
    
OFFtime=[mean(MeanDD(Ages==5 )) mean(MeanDD(Ages==7 & OFFs)) ...
        mean(MeanDD(Ages==12 & OFFs)) mean(MeanDD(Ages>30 & OFFs))];
BItime=[mean(MeanDD(Ages==5 )) mean(MeanDD(Ages==7 & BIs)) ...
        mean(MeanDD(Ages==12 & BIs)) mean(MeanDD(Ages>30 & BIs))];
BIONtime=[mean(MeanDD(Ages==5)) mean(Arbor1DD(Ages==7 & BIs)) ...
        mean(Arbor1DD(Ages==12 & BIs)) mean(Arbor1DD(Ages>30 & BIs))];
BIOFFtime=[mean(MeanDD(Ages==5)) mean(Arbor2DD(Ages==7 & BIs)) ...
        mean(Arbor2DD(Ages==12 & BIs)) mean(Arbor2DD(Ages>30 & BIs))];



plot([5 7 12 35], ONtime, 'b')
hold on
plot([5 7 12 35], OFFtime, 'r')
plot([5 7 12 35], BItime, 'g')
plot([5 7 12 35], BIONtime, 'c')
plot([5 7 12 35], BIOFFtime, 'm')

scatter(Ages(Ages==5),MeanDD(Ages==5),'k')
scatter(Ages(ONs),MeanDD(ONs),'b')
scatter(Ages(OFFs),MeanDD(OFFs),'r')
scatter(Ages(BIs),Arbor1DD(BIs),'c')
scatter(Ages(BIs),Arbor2DD(BIs),'m')

legend('ON','OFF','BI','BI ON arbor', 'BI OFF arbor')



hold off

%% plot BI arbors pair wise Density
figure

%%All
bar([Arbor1DD(BIs); Arbor2DD(BIs)]')
title('All BI arbors Dot Density')
%%P5
subplot(4,1,1)
bar(MeanDD(Ages==5))
title('All P5 arbors Dot Density')
%%P7
subplot(4,1,2)
bar([Arbor1DD(BIs & Ages==7); Arbor2DD(BIs & Ages==7)]',1)
title('P7 BI arbors Dot Density')
%%P12
subplot(4,1,3)
bar([Arbor1DD(BIs & Ages==12); Arbor2DD(BIs & Ages==12)]',1)
title('P12 BI arbors Dot Density')
%%P30+
subplot(4,1,4)
bar([Arbor1DD(BIs & Ages>30); Arbor2DD(BIs & Ages>30)]',1)
title('P30+ BI arbors Dot Density')

%%Test Statistics
%%Is ON more Dense then OFF?
[h BiPairTT_P7]=ttest(Arbor1DD(BIs & Ages==7),Arbor2DD(BIs & Ages==7))
[h BiPairTT_P12]=ttest(Arbor1DD(BIs & Ages==12),Arbor2DD(BIs & Ages==12))
[h BiPairTT_P35]=ttest(Arbor1DD(BIs & Ages>30),Arbor2DD(BIs & Ages>30))



%% plot BI arbors pari wise with Dot Number
figure

bis=find(BIs);
for i = 1 : size(bis,2) %extract Dot number for each arbor
   DotsA1(bis(i))=Master(bis(i)).Results.Arbor(1).Dots;
   DotsA2(bis(i))=Master(bis(i)).Results.Arbor(2).Dots;
end

P5s=find(Ages==5);
for i=1: size(P5s,2)
   DotsP5(i)=Master(P5s(i)).Results.Arbor(1).Dots; 
end


bar([DotsA1(BIs);DotsA2(BIs)]',1)
title('All BI arbors Dot Number')
%%P5
subplot(4,1,1)
bar(DotsP5)
title('All P5 arbors Dot Number')
%%P7
subplot(4,1,2)
bar([DotsA1(BIs & Ages==7);DotsA2(BIs & Ages==7)]',1)
title('P7 BI arbors Dot Number')
%%P12
subplot(4,1,3)
bar([DotsA1(BIs & Ages==12); DotsA2(BIs & Ages==12)]',1)
title('P12 BI arbors Dot Number')
%%P30+
subplot(4,1,4)
bar([DotsA1(BIs & Ages>30); DotsA2(BIs & Ages>30)]',1)
title('P30+ BI arbors Dot Number')

%% Plot Dend Length of BIs
figure

bis=find(BIs);
for i = 1 : size(bis,2) %extract Dot number for each arbor
   DendA1(bis(i))=Master(bis(i)).Results.Arbor(1).Length;
   DendA2(bis(i))=Master(bis(i)).Results.Arbor(2).Length;
end

P5s=find(Ages==5);
for i=1: size(P5s,2)
   DendP5(i)=Master(P5s(i)).Results.Arbor(1).Length; 
end


bar([DendA1(BIs);DendA2(BIs)]',1)
title('All BI arbors Dend Length')
%%P5
subplot(4,1,1)
bar(DendP5)
title('All P5 arbors Dend Length')
%%P7
subplot(4,1,2)
bar([DendA1(BIs & Ages==7);DendA2(BIs & Ages==7)]',1)
title('P7 BI arbors Dend Length')
%%P12
subplot(4,1,3)
bar([DendA1(BIs & Ages==12); DendA2(BIs & Ages==12)]',1)
title('P12 BI arbors Dend Length')
%%P30+
subplot(4,1,4)
bar([DendA1(BIs & Ages>30); DendA2(BIs & Ages>30)]',1)
title('P30+ BI arbors Dend Length')



%% Save Results
save('C:\Notes\Notebook\matlab\Analyze\Data\Gout.mat')




