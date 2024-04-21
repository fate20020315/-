%%%%%%%%%%%%%%%%%%%%%%%%%生成随机梯度牛顿法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddd = minihession(g,h,x,i,eta)  
hh=zeros(71);gg=zeros(1,71);
for j=1:length(i)
    hh=hh+h(x,i(j));
    gg=gg+g(x,i(j));
end
hh=hh/length(i);
gg=(1+eta)*(gg/length(i));
ddd=hh\gg';
end
