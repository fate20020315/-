function y=norm_g(g,x)
w=0;
ns=16033;
for k=1:ns
    w=w+g(x,k);
end
y=norm(w)/ns;