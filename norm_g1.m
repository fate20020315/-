function y=norm_g1(g,x,i)
w=0;
for k=1:length(i)
    w=w+g(x,i(k));
end
y=norm(w)/length(i);