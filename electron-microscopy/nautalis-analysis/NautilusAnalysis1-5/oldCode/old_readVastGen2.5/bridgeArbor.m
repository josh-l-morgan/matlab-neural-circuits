function[] = bridgeArbor(arbor,surfVox);



tips = [arbor.branches.tip];
tipSurf = arbor.nodes.node2vox(tips);
tipSeed = surfVox.seedPath.segID(tipSurf);
sourceSeed = surfVox.seedPath.seed;
parents = [arbor.branches.parent];

sourceList = find(tipSeed == sourceSeed);
breakList = find(tipSeed ~= sourceSeed);





for i = 1
    
end