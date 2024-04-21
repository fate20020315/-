%%%%%%%%%%%%%%%%%%%%%%%%%生成随机梯度牛顿法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gg = minigradient(x,i)  

%ns=1605;
%sns=ceil(0.5*ns);
lamed=0.01;
a1=load('data1a.mat');
b1=load('data1b.mat');
%a1=load('data2222a.mat');
%b1=load('data2222b.mat');
%a1=load('data3a.mat');
%b1=load('data3b.mat');
a=zscore(a1.A);
b=b1.b;
%f=@(x,i)(log(1+exp(-b(i)*(a(i,:)*x'))));
g=@(x,i)(-b(i)*exp(-b(i)*(a(i,:)*x'))*a(i,:)/(1+exp(-b(i)*(a(i,:)*x')))+lamed*x);
%h=@(x,i)(((b(i))^2*exp(-b(i)*(a(i,:)*x')))*(a(i,:)'*a(i,:))/((1+exp(-b(i)*(a(i,:)*x'))))^2+lamed*eye(95));
%hh=zeros(71);gg=zeros(1,71);
gg=zeros(1,71);
%hh=zeros(69);gg=zeros(1,69);
  %i=unidrnd(ns,sns,1,1);
for j=1:length(i)
    gg=gg+g(x,i(j));
end
gg=gg/length(i);

