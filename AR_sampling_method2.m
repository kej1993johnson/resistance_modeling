%Sample x according the pdf p(x)=1/(1+x^2), 0<=x<=1
%x=[0:0.1:1]
%y=1./(1+x.^2)
%plot(x,y)


N=10000;
r1=rand(1,N);
r2=rand(1,N);
px=1./(1+r1.^2);
test=r2<px;
x=r1.*test;
x(x==0)=[];
hist(x)



