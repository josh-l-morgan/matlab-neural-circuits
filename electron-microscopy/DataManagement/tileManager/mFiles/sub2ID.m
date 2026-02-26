function[ID] = sub2ID(w, s, r, c)

IDstr = sprintf('1%03.0f%03.0f%03.0f%03.0f%',w,s,r,c);
ID = str2num(IDstr);






