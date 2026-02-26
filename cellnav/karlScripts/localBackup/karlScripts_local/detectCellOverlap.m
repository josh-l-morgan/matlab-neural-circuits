%
%Plan
% take a list of N cids. Return an N x N matrix showing the size of overlap
% in DS voxels.
% Input: cidList, dsObj, obi or tis
% steps:
%  for each cell, get a list of object IDs
%  get the voxels for each cid
%  do a convex hull on a zprojection
%  compare via same code as overlap

