



    nnDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\nn\';
    expDir = [nnDir expName '\'];
    cellDir = sprintf('%scid%d\\',expDir,runCid);
    expInfoName = sprintf('%sexpInfo.mat',cellDir);




%% Read in the result
load(expInfoName);
numExp = expInfo.numExp;
%synVs = zeros(expNum,numSyn,1);
clear synVs
clear maxVs
clear synIs
%allVs = zeros(numExp,size(nnRes.allV,2),size(nnRes.allV,1));
%maxVs = zeros(expNum,numNodes);
c = 0;
for e = 1:numExp
    expFileName = sprintf('%sex%d_con%02.0f.mat',cellDir,e,con);
    if exist(expFileName,'file')
        c = c+1;
        % disp(sprintf('found %s',expFileName))
        load(expFileName)
        synVs(c,:,:) =  nnRes.synV;
        synIs{c} = nnRes.synI;
        maxVs(e,:) = nnRes.maxV;
        %allVs(c,:,:) = nnRes.allV';
    end
end
numRun = c;
nTime = nnRes.time;






%% Syn with Polarity
if showVideo

    synVs(e,:,:) = synV(uSynIdx,:);
    synMax = max(synVs,[],3);
    synDev = synMax;
    synDev = synDev - median(synDev(:));
    synDev(synDev<0) = 0;
    sumDev = sum(synDev,3);
    allProp = sumDev;

    meanV = squeeze(mean(mean(synVs,1),2));
    difV = abs(meanV-median(meanV))>.01;
    showV = find(difV,1);
    markerType = {'o','o','d'}
    experimentFigureString = {'off', 'on'}

    cla
    hold on
    axis 'equal'
    axis 'off'
    set(fig,'clipping','off')
    hAx = gca;
    set(hAx, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])

    view(00,00)
    s1 = scatter3(tr.X,tr.Y,tr.Z,1,'.','markeredgecolor',[.6 .6 .6]);
    % s1 = scatter3(tr.X,tr.Y,tr.Z,10,'o','markerfacecolor','flat','markerfacealpha',.2,...
    %     'markeredgecolor','none');
    s4 = scatter3(tr.X(nnRes.syn2.nodes),tr.Y(nnRes.syn2.nodes),tr.Z(nnRes.syn2.nodes),...
        50,markerType{3},'k','linewidth',2);

    s3 = scatter3(tr.X(synIDs),tr.Y(synIDs),tr.Z(synIDs),50,'o','markerfacecolor','flat','markerfacealpha',1,...
        'markeredgecolor','none');

    recordMovie = 0;
    movieDir = sprintf('D:\\WorkDocs\\Publications\\VG3\\Revisions\\movies\\allVG3\\cid%d\\',runCid);
    if ~exist(movieDir,'dir'),mkdir(movieDir);end
    c = 0;

    for r = 1


        for i = 1:numRun

            a = activate{i};
            eV = squeeze(synVs(i,:,:));
            eSum = sum(abs(eV-eV(1,1)),1);
            eStart = find(eSum>0,1)
            %eV = squeeze(allVs(i,:,:));
            sV = eV + 90;
            sV = ceil(sV*255/max(sV(:)));
            sV(sV<1)=1;
            cmap = jet(256);
            tV = sV(:,1);
            vColMap = cmap(tV,:);
            s2 = scatter3(tr.X(a),tr.Y(a),tr.Z(a),100,markerType{i},'k','linewidth',2);
            title = experimentFigureString{i};


            for t = eStart:size(sV,2)
                t
                tV = sV(:,t);
                vColMap = cmap(tV,:);
                set(s3,'cdata',vColMap)
                %set(s1,'cdata',vColMap)

                % delete(s3)
                % s3 = scatter3(tr.X(synIDs),tr.Y(synIDs),tr.Z(synIDs),100,'.','cdata',vColMap);
                pause(.01)

                if recordMovie
                    c = c+1;
                    printName  = sprintf('%s%05.0f.png',movieDir,c);
                    print(fig,printName,'-dpng')

                end
            end

            s2.delete;

        end

    end
end
