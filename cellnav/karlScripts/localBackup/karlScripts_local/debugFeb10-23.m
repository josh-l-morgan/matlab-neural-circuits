oldSM2=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs_230201\sm_cid2.mat');
newSM2=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\sm_cid2.mat');

size(oldSM2.sm.arbor.nodes.nodes)
size(newSM2.sm.arbor.nodes.nodes)

size(oldSM2.sm.nep.fv.vertices)
size(newSM2.sm.nep.fv.vertices)

find(any(oldSM2.sm.arbor.edges==23456,2))
find(any(newSM2.sm.arbor.edges==23456,2))

length(unique(oldSM2.sm.arbor.edges,'rows'))
length(unique(newSM2.sm.arbor.edges,'rows'))

length(unique(oldSM2.sm.arbor.edges(:)))
length(unique(newSM2.sm.arbor.edges(:)))

figure(); 
hold on; 
pOld=patch(oldSM2.sm.nep.fv); 
pNew=patch(newSM2.sm.nep.fv);
pOld.FaceAlpha=0.5; 
pNew.FaceAlpha=0.5;
pOld.FaceColor=[1 0 1];
pNew.FaceColor=[0 1 1];
pNew.EdgeColor='none';
pOld.EdgeColor='none';

figure(); 
hold on; 
scatter3(oldSM2.sm.nep.pos(:,1),oldSM2.sm.nep.pos(:,2),oldSM2.sm.nep.pos(:,3),'m.'); 
scatter3(newSM2.sm.nep.pos(:,1),newSM2.sm.nep.pos(:,2),newSM2.sm.nep.pos(:,3),'c.'); 

scatter3(oldSM2.sm.arbor.nodes.pos(:,1),oldSM2.sm.arbor.nodes.pos(:,2),oldSM2.sm.arbor.nodes.pos(:,3),'mo'); 
scatter3(newSM2.sm.arbor.nodes.pos(:,1),newSM2.sm.arbor.nodes.pos(:,2),newSM2.sm.arbor.nodes.pos(:,3),'co'); 


figure(); hold on
plot3(newSM2.sm.arbor.skel.bones(1).nodePos(:,1),newSM2.sm.arbor.skel.bones(1).nodePos(:,2), ...
    newSM2.sm.arbor.skel.bones(1).nodePos(:,3),'m.');
plot3(oldSM2.sm.arbor.skel.bones(1).nodePos(:,1),oldSM2.sm.arbor.skel.bones(1).nodePos(:,2), ...
    oldSM2.sm.arbor.skel.bones(1).nodePos(:,3),'c.');


