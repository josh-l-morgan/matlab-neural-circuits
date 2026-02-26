function loc=getCurLocMed()
d=datacursormode();
inf=getCursorInfo(d);
pos=inf.Position;
loc=pos([3 2 1]);
end