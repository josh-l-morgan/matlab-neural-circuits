function[] = cnvSetUpNrnExperiment(app)

%% set up neuron experiment file to be run in python
global globSC


makeVolMPNcnv
load('MPN.mat')
swcDir = [WPN 'swc\'];

ex.cids = globSC.pickCID; %[1008];
ex.e = str2num(app.reversalPotentialEditField.Value); %0; % reversal potential
ex.gmax = str2num(app.gmaxEditField.Value); %0.8;
ex.tau = str2num(app.tauEditField.Value); %0.1;
ex.Ra = str2num(app.RaEditField.Value); %100; %Axial resistance in Ohm * cm
ex.cm = str2num(app.cmEditField.Value); %1; %# Capacitance in microFarad/cm2
ex.g_pas = str2num(app.g_passEditField.Value); %0.001; %Passive conductance in S/cm2
ex.e_pas = str2num(app.e_passEditField.Value); %-65; 


exFile = sprintf('experiment.mat');

save([swcDir exFile],'ex');


