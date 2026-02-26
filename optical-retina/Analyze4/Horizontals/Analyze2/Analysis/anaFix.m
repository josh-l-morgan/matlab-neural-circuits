function[Result]=anaFix

%%A number of errors were made in some versions of the analysis. 
%%This program searches for those errors and corrects them
%%The first error is that many stacks were read in with starting with two
%%blank planes.  To fix this all data must be shifted by two planes. 

%%The second problem is a single pixel xy shift in the the block buffering
%%this problem also requires the a single pixel shift of all data. 

%%The third problem is D (Thresholded Dendrite) should be in data

'Fixing Data'

global TPN
%%Create directory in which to enter fixed data
Fixed='dataFix';
if ~exist([TPN  Fixed]), mkdir([TPN  Fixed]); end
DataFixTPN=TPN
Result=cell(1,4);
Result(1,:)={0};
Sta=find(TPN=='\');
Sta=Sta(size(Sta,2));
Result(1)={TPN(Sta+1:size(TPN,2))};

%% Move D (also necessary for shift)
if ~exist([TPN 'data\D.mat']) %if D isnt where it should be
    if exist([TPN 'temp\D.mat']) %pull it from temp if possible
        load([TPN 'temp\D.mat'])
        save([TPN 'data\D.mat'],'D')
        clear D
        Result(1,2)={'D moved'};
    elseif exist([TPN 'pics\D'])
        Result(1,2)={'D created'};
        Dd=dir([TPN 'pics\D']); Dd=Dd(3:size(Dd));
        clear D
        for Dp=1:size(Dd,1)
            D(:,:,Dp)=imread([TPN 'pics\D\' Dd(Dp).name]);
        end
        save([TPN 'data\D.mat'],'D'), clear D
      else
        Result(1,2)={'No D'};
    end
else
    Result(1,2)={'D present'};
end




%% Correct the two pixel z shift present in many data files

%%Check number of real image planes
d=dir(TPN);
d=d(3:size(d,1));
Found=0;  %how many did you find
for i=1:size(d,1)
    if ~strcmp(d(i).name,'images')
        Files=dir([TPN  d(i).name]);
        if size(Files,1)>2
        Tag=Files(3).name;
        %look for tifs
        if strcmp('.tif',Tag(1,max(1,size(Tag,2)-3):size(Tag,2)))
            OK=1; %check if all the same
            for f=4:size(Files,1)
              if ~Files(f-1).bytes==Files(f).bytes, 
                    OK=0;
                end
            end %check all files
            if OK, Found=Found+1; DPN=[TPN  d(i).name];end %name path
        end
        end
    end
end

if exist('DPN') %if you found an image directory
    
ImageZsize=size(dir(DPN),1)-2;
zum=.3;

%%Fix Dot Data
if exist([TPN 'data\BigFilled.mat']); %if Dots have been processed
    load([TPN 'data\BigFilled.mat']);
    if size(BigFilled,3)-ImageZsize==2 %check if it has the double offset
        %%Fix BigFilled
        BigFilled=BigFilled(:,:,3:size(BigFilled,3)); %Fixit
        save([TPN  Fixed '\BigFilled.mat'],'BigFilled');
        clear BigFilled
        %%Fix BigCentroid
        load([TPN 'data\BigCentroid.mat']);
        BigCentroid=BigCentroid(:,:,3:size(BigCentroid,3)); %Fixit
        save([TPN  Fixed '\BigCentroid.mat'],'BigCentroid');
        clear BigCentroid
        Result(1,3)={'Dots Fixed'};
    else
        Result(1,3)={'Dots OK'};
    end%if offset
else
    Result(1,3)={'No Dots'};
    
end %if Dots have been processed
clear BigFilled

%%Fix Dend Data
if exist([TPN 'data\D.mat'])
        load([TPN 'data\D.mat'])
        if size(D,3)-ImageZsize==2; %%if shifted by 2
            D=D(:,:,3:size(D,3)); %fix D
            save([TPN  Fixed '\D.mat'])
            clear D
            %fix AllSeg
            if exist([TPN 'data\AllSeg.mat'])
                load([TPN 'data\AllSeg.mat'])
                AllSeg(:,3,:)=AllSeg(:,3,:)-2*zum; %shift all Z segments
                save([TPN Fixed '\AllSeg.mat'],'AllSeg')
                clear AllSeg
            end
            Result(1,4)={'Dend Fixed'};
        else
            Result(1,4)={'Dend OK'};
        end %%If D is shifted by 2
else
    Result(1,4)={'No Dend'};
end %if Dendrites have been processed

end %if an image directory was found

'Done Fixing Data'


