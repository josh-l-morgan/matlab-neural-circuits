function yesno = boundingboxoverlap(minx1,miny1,maxx1,maxy1,minx2,miny2,maxx2,maxy2)

yesno=1;
if (maxx2<minx1) || (maxx1<minx2) || (maxy1<miny2) || (maxy2<miny1) 
  yesno=0;
end;
