function[sm] = getSynSM(sm)

load('MPN.mat')
load([WPN 'tis.mat'])

cid = sm.cid;

isTarg = find((tis.syn.edges(:,1)==cid | tis.syn.edges(:,2) == cid));

sm.syn.synID = isTarg;
sm.syn.edges = tis.syn.edges(isTarg,:);
sm.syn.isIn = find(sm.syn.edges(:,1) == cid);
sm.syn.isOut = find(sm.syn.edges(:,2) == cid);
sm.syn.pre = tis.syn.pre(isTarg);
sm.syn.post = tis.syn.post(isTarg);
sm.syn.obID = tis.syn.obID(isTarg);
sm.syn.synType = tis.syn.synType(isTarg);
sm.syn.preClass = tis.syn.preClass(isTarg);
sm.syn.postClass = tis.syn.postClass(isTarg);
sm.syn.pos = tis.syn.pos(isTarg,:);
sm.syn.synPosRaw = tis.syn.synPosRaw(isTarg,:);
sm.syn.synPosDS = tis.syn.synPosDS(isTarg,:);
sm.syn.order = tis.syn.order;




%% Fetch syn props
isPre = []; isPost = [];
for i = 1: length(tis.syn.synProp)
   sp = tis.syn.synProp{i};
   if isfield(sp,'preID')
      if sum(sp.preID == cid)
          isPre =[isPre i];
      end
      if sum(sp.postID == cid)
          isPost = [isPost i];
      end
      
   end
    
end

sm.syn.syPropsIn = tis.syn.synProp(isPost);
sm.syn.synPropsOut = tis.syn.synProp(isPre);
sm.syn.synProp = tis.syn.synProp(unique([isPre isPost]));




