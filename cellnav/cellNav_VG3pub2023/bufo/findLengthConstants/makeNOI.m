global glob


SPN = [glob.datDir 'Analysis\Data\preproc\'];
load([SPN 'COI.mat']);
clear NOI

cids = COI.vgcCids;
NOI.cids = cids;
%% load data
makeMPNcnv
load('MPN.mat')
swcDir = [WPN 'swc\'];
smDir = [WPN 'SMs\'];
swcDir = [WPN 'swc\'];

for c = 1:length(cids)
    cid = cids(c);
    disp(sprintf('Making NOI for cid%d (%d of %d)',cid,c,length(cids)))

    NOI.good(c) = 1;

    smxName = sprintf('smx_cid%d.mat',cid);
    if exist([smDir smxName],'file')
        smx = load([smDir smxName],'syn','nep','syn2Skel');
        NOI.d{c} = smx.syn2Skel.skelTopoDist;
    else
        smName = sprintf('sm_cid%d.mat',cid);
        if exist([smDir smName])
            load([smDir smName]);
            NOI.d{c} = sm.syn2Skel.skelTopoDist;
        else
            NOI.good(c) = 0;
        end
    end
    
    if NOI.good(c)
        nrnName = sprintf('nrn_cid%d.mat',cid);
        if exist([swcDir nrnName],'file')
            nrn = load([swcDir nrnName]);
            W = nrn.maxVolt;
            swc2arbor = smx.nep.swcS.swc2arborID;
            W = W(:,swc2arbor);
            NOI.isNeu(c) = 1;
            NOI.W{c} = W;
            NOI.p{c}.ex = nrn.ex;
        else
            lc = 20;
            NOI.isNeu(c) = 0;
            NOI.W{c} = exp(-NOI.d{c}/lc); % Apply length constant
            NOI.p{c}.lc = lc;
        end
    end
end


save([WPN 'NOI.mat'],'NOI')



if 0
    %% Pick out variables
    swc = smx.nep.swcS;
    pred = swc.pred;

    %% map synapse to edge
    synPos = smx.syn.pos;
    pos =  smx.nep.pos;
    edges = smx.nep.edges;
    swcPos = swc.pos;
    synSecPos = nrn.synSecPos;
    close = sm.syn2Skel.closest;
    synClosePos = pos(close,:);



    d = smx.syn2Skel.skelTopoDist;
    exAx = gca;
    newID = swc.arbor2swcID;
    d2 = d;
    d2(:,newID) = d;

    pos2 = pos;
    pos2(newID,:) = pos;




    %%%%%%%%%%%%%%%%%%%%%%%%%%% Optional display and testing

    %% Estimate length constant from previous topology
    W2 = W;
    for i = 1:size(W,1)
        w = W(i,:);
        w = w- min(w);
        w = w/max(w);
        W2(i,:) = w;
    end

    scatter(d2(:),W2(:),'.')
    maxD = max(d2(:));
    binW = 1;
    bins = 0:.1:maxD;
    vW = bins * 0;
    for b = 1:length(bins)
        isD = find((d2(:)>(bins(b)-binW)) & (d2(:)<(bins(b)+binW)));
        vW(b) = median(W2(isD));
    end
    plot(bins,vW)
    passL = find(vW <= 0.37,1)
    estL = bins(passL);

    %% Reshape
    if 0
        nrnPos = nrn.secPos;
        nrnPos = nrnPos(:,[2 1 3]);
        meanPos = mean(pos,1);
        nrnSynPos = swcPos(sm2nrn.synCloseEdge,:);
        %nrnSynPos = nrnPos(sm2nrn.synCloseEdge,:);
    end


    %% Find length constant





    if 0
        subplot(1,1,1)
        cla
        exAx = gca

        for s = 1:size(d,1)

            prop2 = d(s,:);
            %propCol2  = (1-exp(-prop2/10));
            propCol2 = 1-prop2;
            propCol2 = propCol2-min(propCol2);
            propCol2 = propCol2*99/max(propCol2)+1;
            propCol2 = propCol2;
            propCol2(propCol2>100) = 100;
            propCol2(propCol2<1) = 1;

            s
            prop = W(s,:);
            propCol  = (1-exp(-prop/10));
            propCol = propCol-min(propCol);
            propCol = propCol*99/max(propCol)+1;
            propCol = propCol;
            propCol(propCol>100) = 100;
            propCol(propCol<1) = 1;

            %%scatSkel2 = scatter3(exAx,pos(:,1),pos(:,2),pos(:,3),'o','k','filled');
            hold off

            scatSkel2 = scatter3(exAx,pos(:,1),pos(:,2),pos(:,3),30,'o','filled','markerfacealph',1);
            set(scatSkel2,'CData',propCol2);%cmap(round(propCol),:);
            hold on
            synPos
            scatSkel5 = scatter3(exAx,synPos(s,1),synPos(s,2),synPos(s,3),530,'p','filled','markerfacealph',1);
            scatSkel5 = scatter3(exAx,synPos(:,1),synPos(:,2),synPos(:,3),330,'k','x','linewidth',3);

            scatSkel3 = scatter3(exAx,swcPos(:,1)+1,swcPos(:,2),swcPos(:,3),30,'d','filled','markerfacealph',1);
            set(scatSkel3,'CData',propCol);%cmap(round(propCol),:);
            hold on

            pause

        end
    end



    if 0
        for s2 = 1:size(W,1)
            startN = ceil(rand* size(W,2));
            for ss = 1;%size(W,2)
                s = startN;
                s1 = s2;
                %% Mark position on arbor
                tA = subplot(2,1,1);
                scatSkel = scatter(tA,swcPos(:,1),swcPos(:,2),10,'k','o','filled','markerfacealph',1);
                hold on
                scatSkel2 = scatter(tA,swcPos(s,1),swcPos(s,2),150,'r','d','filled','markerfacealph',1);
                scatSkel2 = scatter(tA,pos2(s,1),pos2(s,2),150,'g','s','filled','markerfacealph',1);
                scatSynS = scatter3(tA,synClosePos(s1,1),synClosePos(s1,2),synClosePos(1,3),400,'p','k','filled');
                hold off
                %% Show correlation between distance and electrotonic
                cA = subplot(2,1,2);
                scatter(d2(s1,:),W(s2,:),'.')
                hold on
                scatter(cA,d2(s1,s),W(s2,s),'o','filled','r')
                scatter(cA,d2(1,close(s1)),W(1,close(s2)),'o','filled','g')
                hold off
                startN = pred(startN)+1;
                if startN == 1
                    startN = ceil(rand* size(W,2));
                end
                drawnow
                pause(.1)
            end
        end
    end

    %% test skeleton
    if 0
        cla
        pred = swc.pred+1; %adjust index from python to matlab
        nodes = swc.nodes' + 1;
        last = find(pred<1);
        fillin = pred*0>0;
        fillin(last) = 1;

        xE = [pos(2:end,2) pos(pred(2:end),2)];
        yE = [pos(2:end,1) pos(pred(2:end),1)];
        zE = [pos(2:end,3) pos(pred(2:end),3)];

        plot3(yE,xE,zE)

        cla
        hold on
        for i = 2:length(pred)
            plot3([swcPos(i,1) swcPos(pred(i),1)], [swcPos(i,2) swcPos(pred(i),2)], [swcPos(i,3) swcPos(pred(i),3)])
        end


        cla
        hold on
        for i = 1:size(edges,1)
            plot3([pos(edges(i,1),1) pos(edges(i,2),1)], ...
                [pos(edges(i,1),2) pos(edges(i,2),2)], ...
                [pos(edges(i,1),3) pos(edges(i,2),3)])
        end



        while ~isempty(last)
            next = [];
            for n = 1:length(last)
                newNext = find(pred==last(n));
                next = cat(1,next,newNext);
                if 0
                    for i = 1:length(newNext)
                        plotNext = [swcPos(last(n),:);swcPos(newNext(i),:)]
                        plot3(plotNext(:,1),plotNext(:,2),plotNext(:,3))
                        drawnow
                        hold on
                    end
                end
            end
            last = next
            fillin(last) = 1;
            %     scatSkel4 = scatter3(exAx,swcPos(fillin,1),swcPos(fillin,2),...
            %         swcPos(fillin,3),50,'.','k');
            %     drawnow
        end

        predNotNode = setdiff(unique(pred),unique(nodes))
        predBiggerThanNode = find(pred>=nodes)
        extraNodes = length(nodes) - length(unique(nodes))


    end

    %% match nodes
    if 0

        nrnMatch = zeros(size(swcPos,1),1);
        nrnDist = nrnMatch;
        for i = 1:size(swcPos,1)

            dists = sqrt((nrnPos(:,1)-swcPos(i,1)).^2 + ....
                (nrnPos(:,2)-swcPos(i,2)).^2 + ....
                (nrnPos(:,3)-swcPos(i,3)).^2);
            minDist = min(dists);
            nrnDist(i) = minDist;
            nrnMatch = find(dists==minDist,1);

        end
        lostSwc = find(nrnDist>.1)
        lostSwc = setdiff(lostSwc,1)


        %% length rad
        cla,exAx = gca;
        L = sqrt((swcPos(2:end,1)-swcPos(pred(2:end),1)).^2 + ...
            (swcPos(2:end,2)-swcPos(pred(2:end),2)).^2+ ...
            (swcPos(2:end,3)-swcPos(pred(2:end),3)).^2);

        R = swc.rad(2:end);

        scatter(exAx,L,R,'k','.')
        hold on
        scatter(exAx,L(lostSwc-1),R(lostSwc-1),'r','o')
        plot([0 1],[0 1])
        hold off

        %%  show edges


        cla,exAx = gca;
        e1 = 2:size(swcPos,1);
        e2 = pred(2:end);
        plot3(exAx,[swcPos(e1,1) swcPos(e2,1)]', [swcPos(e1,2) swcPos(e2,2)]',...
            [swcPos(e1,3) swcPos(e2,3)]');
        set(gca,'Clipping','off')
        hold on
        scatSkel5 = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),100,'.','k');
        scatSkel5 = scatter3(exAx,swcPos(lostSwc,1),swcPos(lostSwc,2),swcPos(lostSwc,3),100,'+','k');
        hold off



    end


    %% show nrnOrder
    if 0
        cla,exAx = gca;
        hold on
        for i = 1:size(nrnPos,1)

            scatSkel = scatter3(exAx,nrnPos(i,1),nrnPos(i,2),nrnPos(i,3),50,'r','o','filled','markerfacealph',.3);
            try     scatSkel5 = scatter3(exAx,swcPos(i,1),swcPos(i,2),swcPos(i,3),100,'+','k');
            end
            drawnow
            pause
        end

        hold off
    end

    %% Show weights
    if 0
        cla,exAx = gca;
        hold on
        cmap = colorcube(1000);
        for i = 1:size(W,1);
            subplot(size(W,1),1,i),cla

            plot(W(i,:),'color',cmap(ceil(rand*900),:))
        end
        hold off
    end


    if 0
        %% Show result
        cmap = jet(100);
        colormap = cmap
        subplot(1,1,1)
        globSC.ax = gca
        exAx = globSC.ax;
        set(gcf,'color',[.5 .5 .5])
        for s = 1:size(W,1)

            prop = W(s,:);
            propCol  = prop;
            propCol = propCol-min(propCol);
            propCol = propCol*99/max(propCol)+1;
            propCol = propCol;
            propCol(propCol>100) = 100;
            propCol(propCol<1) = 1;

            prop2 = d(s,:);
            propCol2  = (1-exp(-prop2/10));
            propCol2 = propCol2-min(propCol2);
            propCol2 = propCol2*99/max(propCol2)+1;
            propCol2 = propCol2;
            propCol2(propCol2>100) = 100;
            propCol2(propCol2<1) = 1;




            %%scatSkel2 = scatter3(exAx,pos(:,1),pos(:,2),pos(:,3),'o','k','filled');
            hold off
            if 0
                scatSkel2 = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),150,'o','markerfacealph',1);
                set(scatSkel2,'CData',propCol2);%cmap(round(propCol),:);
            else
                scatSkel2 = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),50,'o','k','markerfacealph',1);

            end
            hold on
            set(gca,'Clipping','off')
            view(0,0)

            %%Draw voltage
            if 0
                scatSkel = scatter3(exAx,pos(:,1),pos(:,2),pos(:,3),50,'o','filled','markerfacealph',1);
            else
                scatSkel = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),50,'o','filled','markerfacealph',1);
            end
            set(scatSkel,'CData',propCol);%cmap(round(propCol),:);

            biggest = find(prop == max(prop),1);
            scatSkel5 = scatter3(exAx,swcPos(biggest,1),swcPos(biggest,2),swcPos(biggest,3),...
                500,'o','filled','markerfacealph',.8);
            %scatSkel4 = scatter3(exAx,swcPos(:,1),swcPos(:,2),swcPos(:,3),50,'x','k');
            %scatSkel5 = scatter3(exAx,swcPos(lostSwc,1),swcPos(lostSwc,2),swcPos(lostSwc,3),100,'+','k');

            %%Draw synapses
            %scatSynAll = scatter3(exAx,synPos(:,1),synPos(:,2),synPos(:,3),150,'o','g','filled','markerfacealph',.3);
            %scatSynSN = scatter3(exAx,nrnSynPos(:,1),nrnSynPos(:,2),nrnSynPos(:,3),200,'d','r','filled','markerfacealph',.3);
            try scatSynS = scatter3(exAx,synPos(s,1),synPos(s,2),synPos(s,3),400,'d','k','filled','markerfacealph',.8); end
            try scatSynS = scatter3(exAx,synClosePos(s,1),synClosePos(s,2),synClosePos(s,3),400,'p','k','filled','markerfacealph',.8); end
            %try scatSynSN = scatter3(exAx,nrnSynPos(s,1),nrnSynPos(s,2),nrnSynPos(s,3),600,'d','k','filled','markerfacealph',.8); end
            %     try scatSynSecPos = scatter3(exAx,synSecPos(s,1),synSecPos(s,2),synSecPos(s,3),600,'d','k',...
            %             'markerfacealph',.8,'linewidth',3); end
            %
            title(sprintf('Length constant = %0.3f',estL))

            %scatSynS = scatter3(exAx,swcPos(1,1),swcPos(1,2),swcPos(1,3),200,'o','r','filled','markerfacealpha',.5);
            drawnow
            s
            hold off


            pause(1)

        end

    end

end % run each cid








