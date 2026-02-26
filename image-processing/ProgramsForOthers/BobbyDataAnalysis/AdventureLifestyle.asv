
startPop = 2000;

age = rand(startPop,1)* 100;
male = rand(startPop,1)<.5;

for r = 1: 100000
    age = age + 1;
    old = (1./(101-age)).^1.5;
    plot(old)
    old = rand(length(old),1)<old;
    
    die = old;
    sum(die)/length(die);
    age = age(~die);
    male = male(~die);
    
    repF = sum(~male & (age>15 ) & (age <35)/3);
    repM = sum(male & (age>18) & (age < 60));
    rat = repF/repM;
    births = repF;
    births = round(births * .75); %succesful 
    age(length(age) + births) = -1;
    male = [ male ; rand(births,1)<.5];
    
    pop(r) = length(age);
    plot(pop)
    %hist(age,0:1:101),
    pause(.01)
    
end



%%
live = ones(males,1)>0;
age = ones(males,1);
for i = 1:20
   die = rand(length(live),1)<=.08;
   sum(die)/length(die);
   live(die)= 0 ;
   age(live) = age(live)+1;
end
sum(live)/length(live);
hist(age,min(age):1:max(age))

