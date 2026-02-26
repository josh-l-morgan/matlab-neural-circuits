
clear D Tip siz reps Hyp Live Ext Ret As DCon

sTarget=10000; %Target Size
HomoS=100; %Strength of feedback (1 is maximum, 10 is on in ten switches

startB=4


VarSize=10000;
D=zeros(VarSize,2,2);
D(1:startB,:,:)=0;
Tip=zeros(VarSize,1)>0;
Tip(1:startB)=1;
siz=[200 200];
reps=10000;
Hyp=1;
Live=Tip;
Live(1:startB)=1;
DCon=zeros(VarSize,1)>0;
As=ones(VarSize,1);
As(1:startB)=(2*pi/startB):(2*pi/startB):(2*pi)
Ext=Tip;
Ret=Tip*0;


%make presynaptic field
[Prey Prex]=find(ones(50,50));
Pre=[Prey Prex];
Pre=Pre*4;
PreCon=zeros(size(Pre,1),1)>0;

for r = 1:reps;
    
   
   %Determine tip state
   Tips=find(Tip & ~DCon & Live);



   %% Retract
   for i = 1:size(Tips,1)
   if Ret(Tips(i))==0; %if not retracting, start
            Ret(Tips(i))=(sum(Live)-sTarget)/sum(Live)>(rand*HomoS);
   end    
       
       
       
         if Ret(Tips(i))
            Ret(Tips(i))=fix(rand*1.1)==0 ; %Stop retracting 
            Live(Tips(i))=0;
            RM=D(Tips(i),:,1); %Get Base
            aTip=(D(:,1,2)==RM(1)) & ((D(:,2,2) ==RM(2)) & Live);
            aBase=(D(:,1,1)==RM(1)) & ((D(:,2,1) ==RM(2)) & Live);
            if sum(aTip) && ~sum(aBase)
                Tip(aTip)=1; %enter all new tip positions as tips
                Ret(aTip)=Ret(Tips(i));
                Ext(aTip)=0;
            end
       end
    end

    %clean up info

%     D=D(Live,:,:);
%     Tip=Tip(Live);
%     Live=Live(Live);
%     id = size(D,1);
    
    %%Grow tips
    Tips=find(Tip & ~DCon & Live);
    
    for i = 1:size(Tips,1)
        
       if Ext(Tips(i))==0; %if not retracting, start
            Ext(Tips(i))=(sTarget-sum(Live))/sum(Live)>(rand*HomoS);
        end    
       
        
        
        
        if Ext(Tips(i))
            %make new node
            id = find(Live==0,1);
            if isempty(id),id=size(D,1)+1; end
            y=D(Tips(i),1,2); x = D(Tips(i),2,2);
            A=As(Tips(i))+(randn * 2*pi)/40; %pick random angle
            yd=sin(A) * Hyp;
            xd=cos(A) * Hyp;
            yn= y + yd;
            xn= x + xd;

            %enter new seg
            D(id,1,1)= y; D(id,2,1) = x;
            D(id,1,2)= yn; D(id,2,2) = xn;
            Tip(Tips(i))=0;
            Tip(id)=1;
            Live(id)=1;
            DCon(id)=0;
            As(id)=A;
            Ext(id)=fix(rand*1.1)==0;
            Ext(Tips(i))=0;
            Ret(id)=0;
        
        end %if in extended state

    end

    %Grow inner
    Bases=find(Live & ~DCon);
    Rnd=rand(size(Bases,1),1);
    Pik=Rnd==max(Rnd);
    Piks=Bases(find(Pik));
    for i = 1:size(Piks,1)

        Pn=(sum(Live)-100)/sum(Live)>(rand*1)+.02;
        if rand<Pn ;
            %make new node
            id = find(Live==0,1);
            if isempty(id),id=size(D,1)+1; end
            y=D(Piks(i),1,2); x = D(Piks(i),2,2);
            A=As(Piks(i))+pi/4*(fix(rand*2)*2-1); %pick random angle
            yd=sin(A) * Hyp;
            xd=cos(A) * Hyp;
            yn= y + yd;
            xn= x + xd;

            %enter new seg
            D(id,1,1)= y; D(id,2,1) = x;
            D(id,1,2)= yn; D(id,2,2) = xn;
            Tip(Piks(i))=0;
            Tip(id)=1;
            Live(id)=1;
            DCon(id)=0;
            As(id)=A;
            Ext(id)=1;
            Ret(id)=0;
        end
    end


    %%Check for connections
%     Tips=find(Tip & ~DCon & Live)
%     for i = 1: size(Tips,1)
%         Dists=sqrt((D(i,1,2)-Pre(:,1)).^2 + (D(i,2,2)-Pre(:,2)).^2);
%         Dist=Dists<10;
% 
%         if sum(~PreCon(Dist))
%             DCon(i)=1; %contact Post
%         end
%         PreCon(Dist)=1; %contact Pre
%     end
% 
%     sum(DCon)

    %% Show Result
    Mids=mean(D(Live,:,:),3);
    I=showvec(Mids,siz);
    image(I*100);pause(.01)



end


