%%%%%%%%%%%%%%%%%%%%%%%%%生成全梯度牛顿法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddd = hession(g,h,x,eta)  
%lamed=0.01;
ns=16033;
% a1=load('data1a.mat');
% b1=load('data1b.mat');
% %a1=load('data3a.mat');
% %b1=load('data3b.mat');
% %a1=load('data2222a.mat');
% %b1=load('data2222b.mat');
% a=zscore(a1.A);
% b=b1.b;
% f=@(x,i)(log(1+exp(-b(i)*(a(i,:)*x'))));
% g=@(x,i)(-b(i)*exp(-b(i)*(a(i,:)*x'))*a(i,:)/(1+exp(-b(i)*(a(i,:)*x')))+lamed*x);
% h=@(x,i)(((b(i))^2*exp(-b(i)*(a(i,:)*x')))*(a(i,:)'*a(i,:))/((1+exp(-b(i)*(a(i,:)*x'))))^2+lamed*eye(71));
hh=zeros(71);gg=zeros(1,71);
 %hh=zeros(69);gg=zeros(1,69);
 %hh=zeros(95);gg=zeros(1,95);
for i=1:ns
    hh=hh+h(x,i);
    gg=gg+g(x,i);
end
hh=hh/ns;
gg=(1+eta)*(gg/ns);
%dd=eig(hh);
%hhh=inv(hh);
ddd=hh\gg';

