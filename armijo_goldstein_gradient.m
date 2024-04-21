function alpha = armijo_goldstein_gradient(f, g, xk, k, d)
    gamma = 0.8;
    c = 0.1;
    alpha = 2;
    while f(xk - alpha * d') > f(xk) - c*alpha*g(xk)*d;
        alpha = gamma * alpha;
        if alpha <1e-1; break;end
    end
    while f(xk - alpha * d') < f(xk) - (1-c)*alpha*g(xk)*d;
        alpha = (alpha + alpha / gamma)/2;
         if alpha >1; break;end
    
    
end