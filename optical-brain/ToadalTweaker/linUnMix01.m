function[W res] = linUnMix01(A,B,unMix)

global globNM
showProgress = 1;
%
% Y = (A * Rxb - B * Rxa)./(Rya * Rxb - Ryb * Rxa);
% X = (B * Ryb - A * Rya)./(Rxa * Ryb - Rxb * Rya);

A = double(A);
B = double(B);
if ~exist('unMix','var')
    outID = 1;
else
    outID = unMix.outID;
end

if showProgress
    pAx = {};
    pAx{1} = subplot(1,3,1,'Parent',globNM.app.PreviewplanePanel);
    pAx{2} = subplot(1,3,2,'Parent',globNM.app.PreviewplanePanel);
    pAx{3} = subplot(1,3,3,'Parent',globNM.app.PreviewplanePanel);

    image(pAx{1},A);
    image(pAx{2},B);
    colR = cat(3,A,B,A);
    image(pAx{3},uint8(colR));
    colT = colR * 0;
    colormap(pAx{1},gray(256))
    colormap(pAx{2},gray(256))
    colormap(pAx{3},gray(256))
    drawnow
end




Rxas = [0.5:.1:1]; % make X biggest contributor to A
Rxbs = [0:.1:1];
tracCC = zeros(length(Rxas),length(Rxbs));
for a = 1:length(Rxas);
    disp(sprintf('running ratio %d of %d',a,length(Rxas)))
    for b = 1:length(Rxbs);

        Rxa = Rxas(a);
        Rxb = Rxbs(b);
        Rya = 1-Rxa;
        Ryb = 1-Rxb;

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


        Ys = myNorms(Y,2);
        Xs = myNorms(X,2);

        cc = corrcoef(Ys,Xs);
        tracCC(a,b) = cc(1,2);

        if showProgress

            pcd = (length(Rxbs) * (a-1) + b ) / (length(Rxbs) * length(Rxas)) * 100;
            colT(:,:,1) = Ys;
            colT(:,:,2) = Xs;
            colT(:,:,3) = Ys;
            image(pAx{1},Ys);
            image(pAx{2},Xs);
            image(pAx{3},uint8(colT));
            colormap(pAx{1},gray(256));
            colormap(pAx{2},gray(256));
            globNM.app.textOut.Value = sprintf('Rxa = %0.2f, Rxb =%0.2f, cc = %0.2f, percent done = %0.1f',Rxa,Rxb,cc(1,2),pcd);
            drawnow
        end
    end
end

tracCC(isnan(tracCC)) = 1;
[bestA bestB] = find(abs(tracCC) == min(abs(tracCC(:))),1);

Rxa = Rxas(bestA);
Rxb = Rxbs(bestB);
Rya = 1 - Rxa;
Ryb = 1 - Rxb;


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
Ys = myNorms(Y,2);
Xs = myNorms(X,2);

cc = corrcoef(Ys,Xs);

if showProgress
    colT(:,:,1) = Ys;
    colT(:,:,2) = Xs;
    colT(:,:,3) = Ys;
    image(pAx{1},Ys);
    image(pAx{2},Xs);
    image(pAx{3},uint8(colT));
    colormap(pAx{1},gray(256))
    colormap(pAx{2},gray(256))
    drawnow
end

res.Y = Ys;
res.X = Xs;
res.Rxa = Rxa;
res.Rxb = Rxb; 
res.Rya = Rya;
res.Ryb = Ryb;


rRank = [Rxa Rxb Rya Ryb];
bigR = find(rRank == max(rRank),1);
if (bigR == 1) | (bigR == 4) %x == a
    res.flipChan = 0;
else
    res.flipChan = 1;
end

res.flipChan;

if (outID==1) 
    W = res.Y;
else
    W = res.X;
end


