function[Iout] = uint8keepRat(Iin);



maxI = max(Iin,[],3);
tooHigh = repmat(maxI>255,[1 1 3]);
fixI = repmat(255./maxI,[1 1 3]);
Iout = Iin;
Iout(tooHigh) = Iin(tooHigh).* fixI(tooHigh);
Iout = uint8(Iout);