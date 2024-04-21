function SNM_line_search()
clear all;clc;
a1=load('data1a.mat');
b1=load('data1b.mat');
a=zscore(a1.A);
b=b1.b;
lamed=0.05; % �º��� 0.1��С
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

%%�������������ݶȺ�Hessian��ʽ
function exec(sub)
rho=0.5; 
ns=16033;                              %%%%%%%%%%%%%%%%%%%%%%%%%%%������
sns=ceil(0.05*ns);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%С����ѡȡmini-batches
 eta=10^(-5);
 i=unidrnd(ns,sns,1,1); 
%eta=0;
iter=50;                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������
init=randn(1,71);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ�㶨��   
k=0;q=1;
root=[0.0775937009592298,0.227896537616627,-0.132091539044554,-0.933382269454150,1.24551563831399,-1.25604472596437,0.463100461831563,-0.317333713655938,0.441370453984454,-0.275181481082708,1.17513747488705,-0.0102439447897137,-1.60644717671675,0.669330276832457,0.827009372945023,0.880837440752042,-0.815142030235972,-0.0333057878453161,0.409469890683108,-2.09429163380619,-0.303395971086659,0.0465493743412715,0.157120755678474,0.764570382616161,0.225273046230969,-2.57357853587244,0.612102604515366,1.10835502945658,-1.22963206756836,-0.102745918484122,-0.688334825270670,-0.540160770307248,0.407483633033806,0.955525318072948,-1.45928965865122,0.599707898109061,1.83335581864339,0.698783461031708,-0.679221450722722,0.241488839671320,-0.242597862776445,-0.872821705487147,1.69205797836392,0.0154021340652591,1.85244804274710,2.11588333119277,-1.26513749482666,-0.359445871223749,-1.73706953188132,-0.143905618075390,-0.834064267006352,-1.25980030288843,-0.705869654729154,0.290019681861137,-1.08402763803933,0.597108041269361,0.0962645474917212,-0.365110260645678,-1.24516172693168,0.0735619798521243,-0.342310625683942,0.510266261458581,-0.242883909885955,-0.201737546936022,0.0361013902312149,-0.224905969277045,-0.307955055432012,1.44756339979869,-0.456183529819516,0.778887923440366,-0.602846526578952];
tic                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ʱ�����
while k<iter
    k=k+1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�㷨��ʼ����
%eta=1/k;
  if sns>=15500
  d= hession(g,h,root(k,:),eta) ;
  if sub==1; s(k)=armijos(@(x)norm_g(g,x), root(k,:),k, rho, d') ;  end
  if sub==2; s(k) = armijo_gradient(@(x)f_func(f,x), @(x)f_func(g,x), root(k,:), k, d); end %��d����d'
  if sub==3; s(k) = armijo_wolfe_gradient(@(x)f_func(f,x), @(x)f_func(g,x), root(k,:), k, d); end %��d����d'
  if sub==4; s(k) = armijo_goldstein_gradient(@(x)f_func(f,x), @(x)f_func(g,x), root(k,:), k, d); end %��d����d'
  else
 d= minihession(g,h,root(k,:),i,eta);
  sns=sns+ceil(0.05*ns);
  i=unidrnd(ns,sns,1,1); 
  if sub==1;    s(k)=armijos(@(x)norm_g1(g,x,i), root(k,:),k, rho, d') ;end
 if sub==2;  s(k) = armijo_gradient(@(x)minibatch_f_func(f,x,i), @(x)minibatch_f_func(g,x,i), root(k,:), k, d);end %��d����d'
  if sub==3;  s(k) = armijo_wolfe_gradient(@(x)minibatch_f_func(f,x,i), @(x)minibatch_f_func(g,x,i), root(k,:), k, d);end %��d����d'
 if sub==4;  s(k) = armijo_goldstein_gradient(@(x)minibatch_f_func(f,x,i), @(x)minibatch_f_func(g,x,i), root(k,:), k, d);end %��d����d'
  end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����׼��
    root(k+1,:)=root(k,:)-s(k)*d';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����µĵ�����
 
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
subplot(2,2, sub)
plot(s);
end
end