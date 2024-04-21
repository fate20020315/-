function y = f_func(f,x)
w=0;
ns=16033;
for k=1:ns
    w=w+f(x,k);
end
y= (w)/ns;
