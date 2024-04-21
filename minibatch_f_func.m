function y=minibatch_f_func(g,x,i)
w=0;
for k=1:length(i)
    w=w+g(x,i(k));
end
y= (w)/length(i);