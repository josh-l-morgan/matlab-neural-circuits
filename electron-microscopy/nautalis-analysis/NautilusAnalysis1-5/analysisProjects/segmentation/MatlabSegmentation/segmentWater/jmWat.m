colormap jet(256)

%% dummy variables
C = 4;
B = [ 20 20; 40 40; 20 40; 40 20];

field = rand(100,100);
gKern = gaus3d([15 15 1],3);
A = fastCon(field,gKern)*255;
A = 255-A;
subplot(1,2,1)
image(fitH(A))

%% Create crumb matrix
M = A * 0;
[ys xs] = size(M);

B(B(:)<1) = 1;
B(B(:,1)>ys,1)=ys;
B(B(:,2)>xs,2)=xs;
ind = B(:,1) + (B(:,2)-1)*xs;

ids = 1:size(B,1);
M(ind) = ids;

%% Iterative valleys
low = min(A(:));
high = max(A(:));
V = A;
t = 1;
while 1
    t = t + 1;
    if t>high, break,end
    [Lab numLab] = bwlabel(V<=t,conn);
    image(fitH(Lab)),pause
    for o = 1 : numLab
        vals = M(Lab == o);
        ids = vals(vals>0);
        if length(ids)>1
            V((Lab == o) & (V == t)) = high;
            t = t-1;
        end
        
    end
end









