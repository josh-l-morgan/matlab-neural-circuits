clear all

load('MPN.mat')

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])
load([WPN 'tis.mat'])

fileName = ' G:\IxQ\Matlab\Analysis\stlCellLibrary\dSamp1_1084.stl';
stl = stlread(fileName);
fv.vertices = stl.Points;
fv.faces = stl.ConnectivityList;

[p] = renderFV(fv,[1 1 1],1);































