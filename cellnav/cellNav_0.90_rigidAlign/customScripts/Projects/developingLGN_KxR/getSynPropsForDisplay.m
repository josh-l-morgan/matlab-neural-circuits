function[synProp] = getSynPropsForDisplay(parsedTypes, synGroup)



if strcmp(synGroup,'rgc') %RGC sheath vs giant
synProp(1).name = 'one';
synProp(1).is = (parsedTypes.ps.preSynVal(:,1)==1) & parsedTypes.tp.rgc;
synProp(1).color = 'none';
synProp(1).mark = 'v';
synProp(1).size = 40;
synProp(1).alpha = .7;
synProp(1).jit = [0 0 0];
synProp(1).edgeColor = [1 .1 1];


synProp(2).name = 'two';
synProp(2).is = (parsedTypes.ps.preSynVal(:,1)==2) & parsedTypes.tp.rgc;
synProp(2).color = 'none';
synProp(2).mark = 'v';
synProp(2).size = 120;
synProp(2).alpha = .7;
synProp(2).jit = [0 0 0];
synProp(2).edgeColor = [1 .1 1];


synProp(3).name = 'three';
synProp(3).is = (parsedTypes.ps.preSynVal(:,1)==3) & parsedTypes.tp.rgc;
synProp(3).color = 'none';
synProp(3).mark = 'v';
synProp(3).size = 300;
synProp(3).alpha = .7;
synProp(3).jit = [0 0 0];
synProp(3).edgeColor = [1 .1 1];

synProp(4).name = 'four';
synProp(4).is = (parsedTypes.ps.preSynVal(:,1)==4) & parsedTypes.tp.rgc;
synProp(4).color = 'none';
synProp(4).mark = 's';
synProp(4).size = 400;
synProp(4).alpha = .7;
synProp(4).jit = [0 0 0];
synProp(4).edgeColor = [1 .1 1];


synProp(5).name = 'four';
synProp(5).is = (parsedTypes.ps.isSpine & parsedTypes.tp.rgc);
synProp(5).color = [0 .3 1];
synProp(5).mark = 'o';
synProp(5).size = 20;
synProp(5).alpha = 1;
synProp(5).jit = [0 0 0];
synProp(5).edgeColor = 'none';

synProp(6).name = 'four';
synProp(6).is = (parsedTypes.ps.isSheath & parsedTypes.tp.rgc);
synProp(6).color = 'none';
synProp(6).mark = 'o';
synProp(6).size = 700;
synProp(6).alpha = .6;
synProp(6).jit = [0 0 0];
synProp(6).edgeColor = [1 1 0];

elseif strcmp(synGroup,'mito')

synProp(1).name = 'one';
synProp(1).is = ( parsedTypes.tp.rgc);
synProp(1).color = [0 1 0];
synProp(1).mark = 'o';
synProp(1).size = 200;
synProp(1).alpha = .5;
synProp(1).jit = [0 0 0];
synProp(1).edgeColor = [1 1 1];


synProp(2).name = 'two';
synProp(2).is = (parsedTypes.ps.darkMito);
synProp(2).color = [1 0 0];
synProp(2).mark = 'o';
synProp(2).size = 200;
synProp(2).alpha = .5;
synProp(2).jit = [0 .1 0];
synProp(2).edgeColor = [1 1 1];


synProp(3).name = 'three';
synProp(3).is = (parsedTypes.ps.noMito);
synProp(3).color = [.1 .3 1];
synProp(3).mark = 'o';
synProp(3).size = 200;
synProp(3).alpha = .5;
synProp(3).jit = [0 -.1 0];
synProp(3).edgeColor = [1 1 1];

synProp(4).name = 'four';
synProp(4).is = (parsedTypes.ps.isSpine );
synProp(4).color = [1 0 1];
synProp(4).mark = 'p';
synProp(4).size = 120;
synProp(4).alpha = 1;
synProp(4).jit = [.01 0 0];
synProp(4).edgeColor = 'none';

elseif strcmp(synGroup,'denseVec')

synProp(1).name = 'one';
synProp(1).is = ( parsedTypes.ps.isDenseVec & parsedTypes.ps.darkMito);
synProp(1).color = [0 1 0];
synProp(1).mark = 'o';
synProp(1).size = 200;
synProp(1).alpha = .5;
synProp(1).jit = [0 0 0];
synProp(1).edgeColor = [1 1 1];


synProp(2).name = 'two';
synProp(2).is = ( parsedTypes.ps.isDenseVec & parsedTypes.ps.noMito);
synProp(2).color = [1 0 0];
synProp(2).mark = 'o';
synProp(2).size = 200;
synProp(2).alpha = .5;
synProp(2).jit = [0 .1 0];
synProp(2).edgeColor = [1 1 1];

elseif strcmp(synGroup,'darkMito')

synProp(1).name = 'one';
synProp(1).is = (parsedTypes.ps.darkMito & parsedTypes.ps.notLarge);
synProp(1).color = [0 0 .7];
synProp(1).mark = '^';
synProp(1).size = 30;
synProp(1).alpha = .3;
synProp(1).jit = [0 0 0];
synProp(1).edgeColor = [1 1 1];


synProp(2).name = 'two';
synProp(2).is = (parsedTypes.ps.darkMito & parsedTypes.ps.notSmall);
synProp(2).color = [0 0 .7];
synProp(2).mark = '^';
synProp(2).size = 100;
synProp(2).alpha = .3;
synProp(2).jit = [0 0 0];
synProp(2).edgeColor = [1 1 1];


synProp(3).name = 'three';
synProp(3).is = (parsedTypes.ps.darkMito & parsedTypes.ps.giant);
synProp(3).color = [0 0 .7];
synProp(3).mark = '^';
synProp(3).size = 400;
synProp(3).alpha = .3;
synProp(3).jit = [0 0 0];
synProp(3).edgeColor = [1 1 1];

synProp(4).name = 'four';
synProp(4).is = (parsedTypes.ps.isSpine & parsedTypes.ps.darkMito);
synProp(4).color = [1 0 0];
synProp(4).mark = 'o';
synProp(4).size = 40;
synProp(4).alpha = .8;
synProp(4).jit = [.01 0 0];
synProp(4).edgeColor = [1 1 1];

synProp(5).name = 'five';
synProp(5).is = (parsedTypes.ps.isSheath & parsedTypes.ps.darkMito);
synProp(5).color = [0 1 0];
synProp(5).mark = 'o';
synProp(5).size = 400;
synProp(5).alpha = .2;
synProp(5).jit = [0 .01 0];
synProp(5).edgeColor = [1 1 1];

elseif strcmp(synGroup,'axonType')

synProp(1).name = 'one';
synProp(1).is = (parsedTypes.ps.lightMito);
synProp(1).color = 'none';
synProp(1).mark = 'v';
synProp(1).size = 200;
synProp(1).alpha = .5;
synProp(1).jit = [0 0 0];
synProp(1).edgeColor = [0 1 0];


synProp(2).name = 'two';
synProp(2).is = (parsedTypes.ps.darkMito);
synProp(2).color = 'none';
synProp(2).mark = 'v';
synProp(2).size = 200;
synProp(2).alpha = .5;
synProp(2).jit =  [0 0 0];
synProp(2).edgeColor = [1 0 0];


synProp(3).name = 'three';
synProp(3).is = (parsedTypes.ps.noMito);
synProp(3).color = 'none';
synProp(3).mark = 'v';
synProp(3).size = 200;
synProp(3).alpha = .5;
synProp(3).jit =  [0 0 0];
synProp(3).edgeColor = [.2 .3 1];

synProp(4).name = 'four';
synProp(4).is = (parsedTypes.ps.isSpine );
synProp(4).color = 'none';
synProp(4).mark = 'p';
synProp(4).size = 50;
synProp(4).alpha = .8;
synProp(4).jit = [0 0 0];
synProp(4).edgeColor = [1 1 1];

synProp(5).name = 'four';
synProp(5).is = (parsedTypes.ps.isDenseVec );
synProp(5).color = 'none';
synProp(5).mark = 's';
synProp(5).size = 200;
synProp(5).alpha = .8;
synProp(5).jit =  [0 0 0];
synProp(5).edgeColor = [1 1 1];

synProp(6).name = 'four';
synProp(6).is = (parsedTypes.ps.isSheath );
synProp(6).color = 'none';
synProp(6).mark = 'o';
synProp(6).size = 400;
synProp(6).alpha = .4;
synProp(6).jit = [0 0 0];
synProp(6).edgeColor = [1 1 0];

end
