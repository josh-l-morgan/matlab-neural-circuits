function[] = sm2swc(sm)

%%
%%http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html
%%https://neuroinformatics.nl/HBP/morphology-viewer/


%%Origin node parent must be -1
%%Parent nodes must occur before children nodes

load('MPN')

nodes = sm.nep.nodes(:);
type = nodes * 0;
pos = sm.nep.pos;
rad = sm.nep.nodeRad(:);
edges = sm.nep.edges;
parent = type*0;
parent(sm.nep.edges(:,1)) = sm.nep.edges(:,2);


offset = min(pos,[],1)+1;
offset = mean(pos,1);

pos = pos-repmat(offset,[size(pos,1) 1]);


datMat = cat(2,nodes,type,pos,rad,parent);



%% Reorder Nodes so that parents are always first



counts = hist(parent,nodes)
sum(counts==0)















%% Format swc variable
swcStr = [];
line1 = '# sm to swc';

line2 = sprintf('\n# cell ID %d', sm.cid);
line3 = sprintf('\n# position of cell zero point  = %.3f %.3f %.3f ',...
    offset(1), offset(2), offset(3));



swcStr = [swcStr line1 line2 line3 '\n'];
for i = 1 : size(datMat,1)
    newLine = sprintf('\n%d %d %.3f %.3f %.3f %.3f %d',...
        nodes(i),type(i),pos(i,1),pos(i,2),pos(i,3),rad(i),parent(i));
    swcStr = [swcStr newLine];
end






%%Write SWC

swcDir = [WPN 'swc'];
if ~exist(swcDir,'dir'), mkdir(swcDir); end
fileName = sprintf('%s\\cid%d.swc',swcDir,sm.cid)
fid = fopen(fileName,'w+');
fprintf(fid, swcStr);
fclose(fid);




    




















    
    
    
    