function alpha = armijo_wolfe_gradient(f, g, xk, k, d)
    gamma = 0.8;
    c = 0.1;
    c2 = 0.7;
    alpha = 1;
    while f(xk - alpha * d') > f(xk) - c*alpha*g(xk)*d;
        alpha = gamma * alpha;
    end
    if  g(xk-alpha*d')*d < c2 * g(xk) *d
        alpha = (alpha + alpha / gamma)/2;
    end
    alpha = min(alpha,1);
end