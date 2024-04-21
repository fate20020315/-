function SNM_line_search()
clear all;clc;
a1=load('data1a.mat');
b1=load('data1b.mat');
a=zscore(a1.A);

root=randn(1,71);
b=b1.b;
lamed=0.05; % 新函数 0.1算小
% f=@(x,i)(log(1+exp(-b(i)*(a(i,:)*x')))+lamed*norm(x)^2/2);
% g=@(x,i)(-b(i)*exp(-b(i)*(a(i,:)*x'))*a(i,:)/(1+exp(-b(i)*(a(i,:)*x')))+lamed*x);
% h=@(x,i)(((b(i))^2*exp(-b(i)*(a(i,:)*x')))*(a(i,:)'*a(i,:))/((1+exp(-b(i)*(a(i,:)*x'))))^2+lamed*eye(71));
f = @(x,i)((1-tanh(a(i,:)*x'+b(i)))^2+lamed*norm(x)^2/2);
g = @(x,i)(-2 * a(i,:) * tanh(a(i,:)*x' + b(i)) * (tanh(a(i,:)*x'+b(i))^2-1)+lamed*x);
h = @(x,i)(2 *(a(i,:)'*a(i,:)) * ((tanh(a(i,:)*x'+b(i))^2-1)^2) + 4 *(a(i,:)'*a(i,:)) *(tanh(a(i,:)*x'+b(i))-1) * (tanh(a(i,:)*x'+b(i))^2-1)+lamed*eye(71));
exec(1);
exec(2);
exec(3);
exec(4);

%%产生函数及其梯度和Hessian公式
function exec(sub)
rho=0.5; 
ns=16033;                              %%%%%%%%%%%%%%%%%%%%%%%%%%%样本数
sns=ceil(0.05*ns);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%小样本选取mini-batches
 eta=10^(-2);
 i=unidrnd(ns,sns,1,1); 
%eta=0;
iter=50;                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%迭代次数
init=randn(1,71);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始点定义   
k=0;q=1;

tic                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%时间计算  
while k<iter
    k=k+1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%算法开始运行
%eta=1/k;
  if sns>=15500
  d= hession(g,h,root(k,:),eta) ;
  if sub==1; s(k)=armijos(@(x)norm_g(g,x), root(k,:),k, rho, d') ;  end
  if sub==2; s(k) = armijo_gradient(@(x)f_func(f,x), @(x)f_func(g,x), root(k,:), k, d); end %是d不是d'
  if sub==3; s(k) = armijo_wolfe_gradient(@(x)f_func(f,x), @(x)f_func(g,x), root(k,:), k, d); end %是d不是d'
  if sub==4; s(k) = armijo_goldstein_gradient(@(x)f_func(f,x), @(x)f_func(g,x), root(k,:), k, d); end %是d不是d'
  else
 d= minihession(g,h,root(k,:),i,eta);
  sns=sns+ceil(0.05*ns);
  i=unidrnd(ns,sns,1,1); 
  if sub==1;    s(k)=armijos(@(x)norm_g1(g,x,i), root(k,:),k, rho, d') ;end
 if sub==2;  s(k) = armijo_gradient(@(x)minibatch_f_func(f,x,i), @(x)minibatch_f_func(g,x,i), root(k,:), k, d);end %是d不是d'
  if sub==3;  s(k) = armijo_wolfe_gradient(@(x)minibatch_f_func(f,x,i), @(x)minibatch_f_func(g,x,i), root(k,:), k, d);end %是d不是d'
 if sub==4;  s(k) = armijo_goldstein_gradient(@(x)minibatch_f_func(f,x,i), @(x)minibatch_f_func(g,x,i), root(k,:), k, d);end %是d不是d'
  end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%产生方向  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%步长准则
    root(k+1,:)=root(k,:)-s(k)*d';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成新的迭代点
 
 time(k)=toc;  
end
toc
for j=1:iter+1
  p=zeros(1,71);
  for t=1:ns
  p=p+g(root(j,:),t);
  end
  q(j)=norm(p)/ns;
end
plot(q);
title('梯度范数的变化趋势'); % 添加标题
xlabel('迭代次数'); % 横轴标签
ylabel('梯度范数'); % 纵轴标签
legend('armijios','traditional armijio ','armijios wolfe ','armijios goldstein '); % 为曲线添加图例
hold on;
end
end