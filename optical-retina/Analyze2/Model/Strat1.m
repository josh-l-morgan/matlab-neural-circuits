
Bsize=1000;
P=1; %Probablility of contact
G=fspecial('gaussian',Bsize,Bsize/4);
A=G(:,Bsize/2); %Axon
A=A/max(A);
D=A; %Dendrites
C=A.*D*P;  %Contacts



reps = 100
for r = 1:reps
    D=C;
    C=A.*D*P;

    plot(A,'b')
    hold on
    plot(D,'r')
    plot(C,'g')
    hold off
    pause
    

end