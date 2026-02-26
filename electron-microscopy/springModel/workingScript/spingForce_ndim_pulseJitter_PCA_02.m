%% This is a scratch pad for the calculation of the forces from the
%distances. It relies heavily on the 'pdist' function in order to calculate
%the pairwise differences in the positional matrix.

%The code currently is for a square matrix of set A to set B, it will not
%work for a matrix of interconnected nodes without slight modification to
%the next few sections.

%% Here is where the connectivity matrix is calculated

%The number of nodes needs to be explicitly stated.
nodeNum=120;
%Declare the square matrix containing set A -> set B connectivity.
% load('testMatBig.csv');
% firstFrame=testMatBig;
%making a list of the synapses from the square matrix.
fullsynlist=zeros((nodeNum/2)^2,3);
%This is importing the data from Josh
% load('saveCon.mat');
% firstFrame=saveCon(:,:,59);
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
dimensions = 6;
%set up some random positions for everything
posMat=rand(length(nodeList),dimensions).*100;

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
jitScale = 5;
jitMod = 5;

%parameters for the repPulse
repPulseBool = 1;
repPulseMult = 100;
repPulseMod = 225;

%plotting paramaters
plotMod = 5;
PCAplotMod = 25;
frameCount = 1000;
plotBool=1;
writeBool=1;
lineBool=0;
turnDeg=0;
%remove all connections under this percentile
lineDistCut=50;
%draw lines for every Nth neuron
%can use things like "round((nodeNum/2)/10,0)"
lineNeurPerc=1;
%save the fig files and how often?
figSaveBool=0;
figMod=250;

%does the simulation pause every N frames and spin in a 360 to show the
%figure? Does it draw the lines during this time?
spinBool=1;
spinDeg=9;
spinLineBool=1;
%spinMod should be a multiple of the plotMod
spinMod=200;

%start the graph
% figure('Position',[100 100 800 600]);
h=figure('Position',[100 100 1280 480]);
% g=gca();
subplot(1,3,1);
g=scatter3([1;2],[1;2],[1;2],500,'.k');
subplot(1,3,2);
g=scatter3([1;2],[1;2],[1;2],500,'.k');
hold on;
axis off;
axis vis3d;
subplot(1,3,3);
plot(1,1);
drawnow;

%initiate the other things
i=1;
jpglist={};
logfile='timeLog.txt';
debug=0;

conListBig=conListReal(conListReal(:,3)>prctile(conListReal(:,3),lineDistCut),:);
colormap(parula);
testColor=ind2rgb((conListBig(:,3).*(1/(round(max(conListBig(:,3)),5)))),jet);

%copypasta
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(conListBig(:,3)),max(conListBig(:,3)),L),1:L,conListBig(:,3)));
testColor = reshape(C(Gs,:),[size(Gs) 3]);


% log=fopen(logfile,'w');
% fprintf(log,'Time for each step in wave2\n');
% fclose(log);

%set up the video file for recording
if writeBool==1
    vidObj=VideoWriter('testing4.avi');
    vidObj.FrameRate=15;
    open(vidObj);
end

%% Loop

for frame=1:frameCount
%     tic;
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
        if repPulseBool==1 && mod(frame,repPulseMod)==0
            repConstTemp=repConst*repPulseMult;
        else
            repConstTemp=repConst;
        end
        currRepForceArray = (repConstTemp.*(-currDispArray)./distCb);
        currRepForceArray(currRepForceArray>500,:)=500;
        currRepForceArray(currRepForceArray<(-500),:)=(-500);
        %Get the connective forces
        currAttForceArray = currDispArray.*connDistDist.*conConst;
        currAttForceArray(currAttForceArray>500,:)=500;
        currAttForceArray(currAttForceArray<(-500),:)=(-500);
        %Get the net force
        currNetForceArray=currRepForceArray-currAttForceArray;
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
        
        % Diagnostics
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
        
        
        % write the output to the posMat
        
        posMat(:,ndim)=currDimPos+deltaDimPos;
        

    end
    
    %Here is the jitter
    
    if mod(frame,jitMod)==0
%         tic
        jitMat=rand(nodeNum,dimensions).*jitScale;
        posMat=posMat+jitMat;
%         toc
    end
    
%     log=fopen(logfile,'a');
%     fprintf(log,'frame %i - calcsDone: %s\r\n',frame,num2str(toc,'%.3f'));
%     fclose(log);

    if mod(frame,plotMod)==0 && plotBool==1
        
%         tic
%         colormap(hsv);
        clf();
        subplot(1,3,2);
        [ccoeff,cscore,clatent,ctsquared,cexplained,cmu] = pca(posMat(:,2:end));
        PCApos=cscore(:,1:3).*repmat(((cexplained(1:3,:)'./100)), [length(cscore(:,1)) 1]);
        g=scatter3(PCApos(:,1),PCApos(:,2),PCApos(:,3),500,'.');
        maxRange=max(range(PCApos(:,:)));
%         limits=[min(PCApos(:,1)) min(PCApos(:,1))+maxRange min(PCApos(:,2)) min(PCApos(:,2))+maxRange min(PCApos(:,3)) min(PCApos(:,3))+maxRange];  
        limits=[mean(PCApos(:,1))-(maxRange/2) mean(PCApos(:,1))+(maxRange/2) mean(PCApos(:,2))-(maxRange/2) mean(PCApos(:,2))+(maxRange/2) mean(PCApos(:,3))-(maxRange/2) mean(PCApos(:,3))+(maxRange/2)];
        axis(limits);
        axis off;
%         axis tight;
%         g=scatter3(posMat(:,1),posMat(:,2),posMat(:,3),500,'.k'); %posMat(:,4),
        if lineBool==1 || (mod(frame,spinMod)==0 && spinLineBool==1)
            for l=1:(nodeNum/2)
                if mod(l,lineNeurPerc)==0
                    preID=nodeList(l,:);
                    synList=conListBig(conListBig(:,1)==preID,:);
                    plotCol=testColor(conListBig(:,1)==preID,:,:);
                    for s=1:length(synList(:,1))
                        pre=(synList(s,1));
                        post=(synList(s,2));
                        prepos=PCApos(pre,1:3);
                        postpos=PCApos(post,1:3);
                        prePostPos=[prepos ; postpos];
                        g=line(prePostPos(:,1),prePostPos(:,2),prePostPos(:,3),'Color',plotCol(s,:),'LineWidth',1);  %conListReal(l,3),
                        
                    end
                end
                        %post=conListBig(l,2);
                %shortened this to just use the indices that I have. Need
                %to make some from IDs if ID~=index.
            end
        end
        
        
        
        subplot(1,3,1);
        g=scatter3(posMat(:,1),posMat(:,2),posMat(:,3),500,'.');
        axis off;
        axis tight;
%         g=scatter3(posMat(:,1),posMat(:,2),posMat(:,3),500,'.k'); %posMat(:,4),
        if lineBool==1 || (mod(frame,spinMod)==0 && spinLineBool==1)
            for l=1:(nodeNum/2)
                if mod(l,lineNeurPerc)==0
                    preID=nodeList(l,:);
                    synList=conListBig(conListBig(:,1)==preID,:);
                    plotCol=testColor(conListBig(:,1)==preID,:,:);
                    for s=1:length(synList(:,1))
                        pre=(synList(s,1));
                        post=(synList(s,2));
                        prepos=posMat(pre,1:3);
                        postpos=posMat(post,1:3);
                        prePostPos=[prepos ; postpos];
                        g=line(prePostPos(:,1),prePostPos(:,2),prePostPos(:,3),'Color',plotCol(s,:),'LineWidth',1);  %conListReal(l,3),
                        
                    end
                end
                        %post=conListBig(l,2);
                %shortened this to just use the indices that I have. Need
                %to make some from IDs if ID~=index.
            end
        end
           
        
        
        drawnow;
        camorbit(turnDeg,0,'data',[0,1,0])
        axis off;

        subplot(1,3,3);
        pareto(cexplained);

        
        if figSaveBool==1 && mod(frame,figMod)==0
            figname=sprintf('figFile_frame_%i',frame);
            csvwrite(sprintf('figPos_frame_%i',frame),posMat);
            savefig(figname);
            matName=sprintf('posMat_frame_%i',frame);
            save(matName,'posMat');
        end         
        if writeBool==1
            writeVideo(vidObj,getframe(gcf));
        end
        
        if spinBool==1 && mod(frame,spinMod)==0
            for p=1:2
                subplot(1,3,p)
                axis vis3d;
                figCol=gcf;
                figCol.Color=[0.80 0.80 1];
                for sp=1:(720/spinDeg)
                    camorbit(spinDeg,spinDeg*0.5,'data',[0,1,0])
                    drawnow;
                    writeVideo(vidObj,getframe(gcf));
                end
                figCol.Color=[0.94 0.94 0.94];
            end
        end
%         log=fopen(logfile,'a');
%         fprintf(log,'frame %i - plotsDone: %s\r\n',frame,num2str(toc,'%.3f'));
%         fclose(log);
    end
end

close(vidObj);

if writeBool==1000
    writerObj = VideoWriter('testout2.avi');
    writerObj.FrameRate=15;
    open(writerObj);
    % pnglist2=string(pnglist);
    for K = 1 : length(jpglist)
      filename = jpglist{K};
    %   fnstr=cellstr(filename);
      thisimage = imread(filename);
      writeVideo(writerObj, thisimage);
    end
    
    close(g);
    close(writerObj);
end

%% This is for speed testing

% 'Just pdist'
% tic
%     nodeListCombArray = nchoosek(nodeList,2);
%     nodeListCombArrayInd = nchoosek(1:length(nodeList),2);    
%     dist=pdist(oldPosMat,'euclidean')';
%     distSq=pdist(oldPosMat,'squaredeuclidean')';
%     distCb=dist.*distSq;
% toc
% 
% %Including the indices takes about 0.03 seconds extra
% 'Pdist plus indices'
% tic
%     nodeListCombArray = nchoosek(nodeList,2);
%     nodeListCombArrayInd = nchoosek(1:length(nodeList),2);
%     disArray = [nodeListCombArrayInd nodeListCombArray pdist(posMat,'euclidean')'];
%     disArraySq = [nodeListCombArrayInd nodeListCombArray pdist(posMat,'squaredeuclidean')'];
%     disArrayCb = [nodeListCombArrayInd nodeListCombArray disArray(:,5).*disArraySq(:,5)];
% toc
% 
% 




