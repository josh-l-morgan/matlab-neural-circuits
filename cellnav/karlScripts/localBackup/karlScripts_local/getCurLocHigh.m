function loc=getCurLocHigh()
d=datacursormode();
inf=getCursorInfo(d);
pos=inf.Position;
loc=uint16(pos([3 2 1]).*[250 250 25]);
end