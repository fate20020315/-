%%%%%%%%%%%%%%%%%%%%%%%%%生成全梯度牛顿法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddd = hession(g1,h1,g2,h2,x,eta,ns,n)  
% lamed=0.01;
 %ns=2353;
% %a1=load('data1a.mat');
% %b1=load('data1b.mat');
% %a1=load('data3a.mat');
% %b1=load('data3b.mat');
% a1=load('data2222a.mat');
% b1=load('data2222b.mat');
% a=zscore(a1.A);
% b=b1.b;
% %f=@(x,i)(1-tanh(b(i)*(a(i,:)*x'))+lamed*norm(x)^2/2);
% g=@(x,i)(-b(i)*(1-(tanh(b(i)*(a(i,:)*x')))^2)*a(i,:)+lamed*x);
% h=@(x,i)(2*(b(i))^2*(tanh(b(i)*(a(i,:)*x')))*(1-(tanh(b(i)*(a(i,:)*x')))^2)*(a(i,:)'*a(i,:))+lamed*eye(95));
% %hh=zeros(71);gg=zeros(1,71);
%  %hh=zeros(69);gg=zeros(1,69);
% hh1=zeros(1,n);gg1=0; hh2=zeros(1,n);gg2=0;
hh1=zeros(n,n);gg1=zeros(n,1); hh2=zeros(n,n);gg2=zeros(n,1);
for i=1:ns
    hh1=hh1+h1(x,i);
    gg1=gg1+g1(x,i)';
     hh2=hh2+h2(x,i);
    gg2=gg2+g2(x,i)';
end
hh1=hh1/ns;
gg1=gg1/ns;
hh2=hh2/ns;
gg2=gg2/ns;
 hh=[hh1;hh2];
 gg=[gg1;gg2];
% hh=hh1;
% gg=gg1;
ddd=(hh'*hh+eta*eye(n))\(hh'*gg);
%ddd=hh\gg';

