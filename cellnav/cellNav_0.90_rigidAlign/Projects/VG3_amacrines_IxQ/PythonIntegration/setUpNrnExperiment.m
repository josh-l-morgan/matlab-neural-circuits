%% set up neuron experiment file to be run in python




load('MPN.mat')

swcDir = [WPN 'swc\'];

ex.cids = [1008];
ex.e = 0; % reversal potential
ex.gmax = 0.8;
ex.tau = 0.1;
ex.Ra = 100; %Axial resistance in Ohm * cm
ex.cm = 1; %# Capacitance in microFarad/cm2
ex.g_pas = 0.001; %Passive conductance in S/cm2
ex.e_pas = -65; 


exFile = sprintf('experiment.mat');

save([swcDir exFile],'ex');


