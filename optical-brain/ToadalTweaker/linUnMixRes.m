function[Y X] = linUnMix01(A,B,unMix)


A = double(A);
B = double(B);

R = unMix.R;
Rxa = R.Rxa;
Rxb = R.Rxb;
Rya = R.Rya;
Ryb = R.Ryb;

outID = unMix.outID;

Y = (A * Rxb - B * Rxa)./(Rya * Rxb - Ryb * Rxa);
X = (A * Ryb - B * Rya)./(Rxa * Ryb - Rxb * Rya);
if Rxa==0
    Y = A;
end
if Rxb == 0;
    Y = B;
end
if Rya ==0
    X = A;
end
if Ryb == 0;
    X = B;
end
Ys = Y;%myNorms(Y,2);
Xs = X;%myNorms(X,2);

cc = corrcoef(Ys,Xs);

if outID == 1
    Y = Ys;
    X = Xs;
else
    Y = Xs;
    X = Ys;
end


