function ret=getinheritedfromids(data,sourceidlist)
%Uses the 24-column-data matrix as analyzed from VAST color files
%Gets a list of IDs of the displayed (collapsed) segments in VAST of the segments in sourceidlist

ret=data(sourceidlist(:),18);
