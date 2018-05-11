%Sample x according the pdf p(x)=1/(1+x^2), 0<=x<=1
%x=[0:0.1:1]
%y=1./(1+x.^2)
%plot(x,y)--- This is the anlystical pdf 
%and what the sample histgram below should look like 


N=100000;
r1=rand(1,N);
r2=rand(1,N);
px=1./(1+r1.^2);
test=r2<px;

j=1;
for ii=1:length(test)
    if test(ii)
        x(j)=r1(ii);
        j=j+1;
    end
end
hist(x)



