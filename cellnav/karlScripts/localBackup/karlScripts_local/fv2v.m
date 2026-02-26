function locVast=fv2v()
% takes the coordinates from FV space and spits out high res vast
d=datacursormode;
xyz = getfield(getCursorInfo(d),'Position');
locV=xyz([3 2 1]).*[250 250 25];
locOut=uint16(locV)
clipboard('copy',locOut) ;
