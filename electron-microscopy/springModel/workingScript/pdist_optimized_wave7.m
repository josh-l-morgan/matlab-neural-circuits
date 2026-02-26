%% This is a scratch pad for the calculation of the forces from the
%distances. It relies heavily on the 'pdist' function in order to calculate
%the pairwise differences in the positional matrix.

%The code currently is for a square matrix of set A to set B, it will not
%work for a matrix of interconnected nodes without slight modification to
%the next few sections.

%% Here is where the connectivity matrix is calculated

%The number of nodes needs to be explicitly stated.
nodeNum=60;
%Declare the square matrix containing set A -> set B connectivity.
load('testMat.csv');
firstFrame=testMat;
%making a list of the synapses from the square matrix.
fullsynlist=zeros((nodeNum/2)^2,3);
%This is importing the data from Josh
% load('saveCon.mat');
% firstFrame=saveCon(:,:,1);
%This takes the syn values from the square matrix and puts them into the
%list format.
for f=1:(nodeNum/2)
    for g=1:(nodeNum/2)
        fullsynlist(((f-1)*(nodeNum/2))+g,:)=[f g firstFrame(f,g)];
    end
end
%Adding number to half the nodes for unique indices
fullsynlist(:,2)=fullsynlist(:,2)+(nodeNum/2);
%get a list of the non-zero connections for drawing lines later.
conListReal=fullsynlist(fullsynlist(:,3)~=0,:);
%Get a trimmed list that only includes the IDs of the syn parts.
fullsynlist_trm=(fullsynlist(:,1:2));
%create the nodeList that will be used for the rest of the script
nodeList = unique(fullsynlist_trm(:,:));
%get the total number of nodes
nodeNum = length(nodeList);
%get the unique connection rows (pre post). 
uniqueConn = unique(fullsynlist_trm,'rows');
%make a pairwise list of possible synapses. All of them.
nodeListCombArray = nchoosek(nodeList,2);
%Comment one of the following two lines, depending on whether the node IDs
%are the same as the indices.
% nodeListCombArrayInd = nchoosek(1:length(nodeList),2);
nodeListCombArrayInd = nodeListCombArray;

% %% This is a slow way to make a connectivity matrix
% % Be careful with this. It destroys the directionality of synapses.
% fullsynlist_trm_ord=fullsynlist_trm;
% %go through and sort each of the synapses so that the smaller index is
% %first.
% for s=1:length(fullsynlist_trm_ord)
%     fullsynlist_trm_ord(s,:)=sort(fullsynlist_trm_ord(s,:));
% end
% 
% [synList,b,c]=unique(fullsynlist_trm_ord,'rows');
% f=histc(c,1:numel(b));
% synHist=[synList f];

% %% This part is completely mad.
% 
% %Get the number of synapses for each of those
% %WARNING: THIS LOOP TAKES ABOUT 10 MINUTES TO RUN. CONSIDER PARALLELIZATION
% %OR CHANGING THE APPROACH.
% connDist=[nodeListCombArray nodeListCombArrayInd zeros(length(nodeListCombArray),1)];
% tic
% for n=1:length(fullsynlist_trm_ord)
%     %order the cells for matching the array. WARNING: THIS DESTROYS DIRECTION OF
%     %SYNAPTIC CONTACT!!
%     matchSyn = sort(fullsynlist_trm_ord(n,:));
%     connDist(ismember(connDist(:,1:2),matchSyn,'rows'),5)=(connDist(ismember(connDist(:,1:2),matchSyn,'rows'),5)+1);
% %     connDist(n,3)=sum(ismember(fullsynlist_trm,nodeListCombArray(n,:),'rows'));
% end
% connDist(:,2)=connDist(:,2)-(nodeNum/2);
% 
% toc

%% This is particular to Josh's data format. I need to put the conn strengths into the connDist
%create the con strength array.
connDistDist=zeros(length(nodeListCombArray),1);
%fill up the array with the strengths of the connections
tic
for f=1:(nodeNum/2)
    %The second synapse
    for g=1:(nodeNum/2)
        connDistDist( (((nodeNum-2)*(f-1))-sum(1:(f-2))+(nodeNum/2)+g-1) ,1)=firstFrame(f,g);        
    end
%     disp(toc);
end
toc

%% Set up the dimensions and such

%set the number of dimensions
dimensions = 4;
%set up some random positions for everything
posMat=rand(length(nodeList),dimensions).*50;

%% This is where the static and dynamic nodes are defined

nodeListStatic=[];
nodeListNoForce=[];

%% Prams for the loop

%Note on constants: (Cr/Cc)^(1/3)=optimum distance.
repConst = 1000;
conConst = 1;
Rreducer = 100;
Creducer = 100;
speedGov = 100;

plotMod = 1;
frameCount = 2000;
hold on;
h=plot([1;1],[2;2]);
i=1;
jpglist={};
plotBool=1;
writeBool=1;
lineBool=1;
figSaveBool=0;
logfile='timeLog.txt';
debug=0;

colormap(jet);
testColor=ind2rgb((conListReal(:,3).*10),jet);

log=fopen(logfile,'w');
fprintf(log,'Time for each step in wave2\n');
fclose(log);

%% Loop
for frame=1:frameCount
    tic
    fprintf('%i\n',frame);
    %get a new matrix started to put the new positions into
    oldPosMat=posMat;
    %get the distances and squared distances
    dist=pdist(oldPosMat,'euclidean')';
    distSq=pdist(oldPosMat,'squaredeuclidean')';
    distCb=dist.*distSq;

    parfor ndim=1:dimensions
        
        %Get the positions, slice for the parfor
        currDimPos=oldPosMat(:,ndim);
        %calculate the differences in positions for all of them.
        % Note: this is a potential slow spot because I power and compare
        % and all that.
        currDispArray = (-1).*((-1).^(currDimPos(nodeListCombArrayInd(:,1),:)<currDimPos(nodeListCombArrayInd(:,2),:)).*(pdist(currDimPos,'euclidean')'));
        %Get the repulsive forces
        currRepForceArray = (repConst.*(-currDispArray)./distCb);
        currRepForceArray(currRepForceArray>100,:)=100;
        currRepForceArray(currRepForceArray<(-100),:)=(-100);
        %Get the connective forces
        currAttForceArray = currDispArray.*connDistDist.*conConst;
        currAttForceArray(currAttForceArray>100,:)=100;
        currAttForceArray(currAttForceArray<(-100),:)=(-100);
        %Get the net force
%         currNetForceArray=currRepForceArray-currAttForceArray;
        deltaDimArray=currRepForceArray+currAttForceArray;
        %speedGoc is a stand-in for actual dynamic oscillation and
        %freak-out control.
        deltaDimArray=deltaDimArray./speedGov;
        %Josh's wizardry to get the lists back to a square.
        sqrInd = sub2ind([nodeNum nodeNum],nodeListCombArrayInd(:,1),nodeListCombArrayInd(:,2));
        con = zeros(nodeNum,nodeNum); 
        con(sqrInd) = deltaDimArray(:,:); %connDistDist(:,:);
        %gotta invert to get the other part of the matrix
        con=con-con';
        %get the sum of the position changes.
        deltaDimPos=sum(con,2);
        
        %% Diagnostics
        if debug==1
            repSquare=zeros(nodeNum,nodeNum);
            conSquare=zeros(nodeNum,nodeNum);
            netSquare=zeros(nodeNum,nodeNum);
            repSquare(sqrInd)=currRepForceArray(:,:);
            conSquare(sqrInd)=currAttForceArray(:,:);
            netSquare(sqrInd)=currNetForceArray(:,:);
            repSquare=repSquare-repSquare';
            conSquare=conSquare-conSquare';
            netSquare=netSquare-netSquare';

            HeatMap(repSquare)
            HeatMap(conSquare)
            HeatMap(netSquare)
        end
        
        
        %% write the output to the posMat
        
        posMat(:,ndim)=currDimPos+deltaDimPos;
        

    end
    
%     log=fopen(logfile,'a');
%     dimTimeList=reshape(dimTime,[],1);
%     for k=1:length(dimTimeList)
%         fprintf(log,'%s\r\n',string(dimTimeList{k}));
%     end
%     fclose(log);

    if mod(frame,plotMod)==0 && plotBool==1
        hold on;
        cla(h);
        tic
%         colormap(hsv);
        h=scatter3(posMat(:,1),posMat(:,2),posMat(:,3),500,'.k'); %posMat(:,4),
%         h=scatter3(posMat(7,1),posMat(7,2),posMat(7,3),1000,'.r');
        if lineBool==1
            for l=1:length(conListReal(:,1))
                pre=conListReal(l,1);
                post=conListReal(l,2);
                %shortened this to just use the indices that I have. Need
                %to make some from IDs if ID~=index.
                prepos=posMat(pre,:);
                postpos=posMat(post,:);
                prePostPos=[prepos ; postpos];
                h=line(prePostPos(:,1),prePostPos(:,2),prePostPos(:,3),'Color',testColor(l,:),'LineWidth',2);  %conListReal(l,3),
            end
        end
        drawnow;
        axis tight;
        camorbit(2,0,'data',[0,1,0])
        if figSaveBool==1
            figname=sprintf('figFile_frame_%i',frame);
            savefig(figname);
        end
%         axis off;
        if writeBool==1
            jpgname = sprintf('new_plotter_1_%i.png',frame);
            print(jpgname,'-dpng');
            jpglist{i}=jpgname;
            i=i+1;
        end
%         log=fopen(logfile,'a');
%         fprintf(log,'frame %i - plotPos: %s\r\n',frame,num2str(toc,'%.3f'));
%         fclose(log);
    end
    toc
end

if writeBool==1
    writerObj = VideoWriter('testout.avi');
    writerObj.FrameRate=15;
    open(writerObj);
    % pnglist2=string(pnglist);
    for K = 1 : length(jpglist)
      filename = jpglist{K};
    %   fnstr=cellstr(filename);
      thisimage = imread(filename);
      writeVideo(writerObj, thisimage);
    end
    close(writerObj);
end

%% This is for speed testing

'Just pdist'
tic
    nodeListCombArray = nchoosek(nodeList,2);
    nodeListCombArrayInd = nchoosek(1:length(nodeList),2);    
    dist=pdist(oldPosMat,'euclidean')';
    distSq=pdist(oldPosMat,'squaredeuclidean')';
    distCb=dist.*distSq;
toc

%Including the indices takes about 0.03 seconds extra
'Pdist plus indices'
tic
    nodeListCombArray = nchoosek(nodeList,2);
    nodeListCombArrayInd = nchoosek(1:length(nodeList),2);
    disArray = [nodeListCombArrayInd nodeListCombArray pdist(posMat,'euclidean')'];
    disArraySq = [nodeListCombArrayInd nodeListCombArray pdist(posMat,'squaredeuclidean')'];
    disArrayCb = [nodeListCombArrayInd nodeListCombArray disArray(:,5).*disArraySq(:,5)];
toc






