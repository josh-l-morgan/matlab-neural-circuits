%%Instructions for comparing skeleton overlap to connectivity.

%{

Some functions take from readVast / analyzeMorphology / skelOverlap

After cells are skeletonized,

run overlapSkel... to find out how much each cell overlaps


1) export cells and read into mat folder
2) skeletonizes cells
3) run filteredSkelOverlap07 to calculate overlap based on node distances
and to find node distance for each synapse
4)overlapPredictSynOfSub can produce a synapse prediction from overlap at a
desired resolution
5) predictPreferenceFromOverlap finds seed cell preference
6) conError finds the difference (prediction error) betwee a predicted
connectivity graph and the actual connectivity graph. 
