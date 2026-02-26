function[] = prepareNeps(cids,app)


load('MPN.mat')
load([WPN 'tis.mat']);

nepDir = [WPN 'neps\'];
smDir = [WPN 'SMs\'];
swcDir = [WPN 'swc\'];

for c = 1 : length(cids)
    cid = cids(c);
    str = sprintf('running cell %d, (%d of %d)',...
        cid,c,length(cids));
    app.textOut.Value = str; drawnow
    
    smxName = sprintf('smx_cid%d.mat',cid);
    smx = load([smDir smxName],'syn','nep','syn2Skel');
    clear sm2nrn
    
    sm2nrn.synCloseNode = smx.syn2Skel.closest;
    sm2nrn.syn = smx.syn;
    sm2nrn.nep = smx.nep;
    sm2nrn.synCloseEdge = []%smx.syn2Skel.closest;
    
    if 1 %% map to SWC instead
        %%map synapse to edge
        synPos = smx.syn.pos;
        %pos2 =  smx.nep.pos;
        pos = smx.nep.swcS.pos;
        
        closeEdge = zeros(size(synPos,1),1);
        for s = 1:size(synPos,1)
            dists = sqrt((synPos(s,1)-pos(:,1)).^2 + ...
                (synPos(s,2)-pos(:,2)).^2 + (synPos(s,3)-pos(:,3)).^2);
            closeEdge(s) = find(dists==min(dists),1);
        end
        
        if 0
            for s = 1:size(synPos,1)
                scatter(pos(:,1),pos(:,2),'.','k')
                hold on
                scatter(synPos(s,1),synPos(s,2),'o','filled','g','markerfacealpha',.5)
                scatter(pos(closeEdge(s),1),pos(closeEdge(s),2),'o','k')
                hold off
                pause(.01)
            end
        end
        
        sm2nrn.synCloseEdge = closeEdge;
    end
    
    sm2nrnFile = sprintf('sm2nrn_cid%d.mat',cid);
    save([swcDir sm2nrnFile],'sm2nrn');
    
end

app.textOut.Value = 'Finished preparing NEPs'; drawnow






