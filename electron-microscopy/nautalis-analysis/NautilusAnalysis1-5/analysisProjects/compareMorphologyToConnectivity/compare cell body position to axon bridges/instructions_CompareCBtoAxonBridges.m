%%The goal is to find the distribution of axonal links from a seed
%%thalamocortical cell to surrounding thalamocortical cells using a
%%combination of network tracing and labeling of surrounding nuclei.

%{

1) Run readCBstack to turn Vast export image stack into structure 
2) Run connectVScb - to map traced neurons onto map of nuclei in volume
3) Run compare CBtoBridge.m to run montecarlo estimating probability of
bridge clustering

%}