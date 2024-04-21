function alpha = armijo_gradient(f, g, xk, k, d)
    gamma = 0.8;
    c = 0.3;
    alpha = 1;
    while f(xk - alpha * d') > f(xk) - c*alpha*g(xk)*d;
        alpha = gamma * alpha;
        if alpha <1e-1; break;end
    end
end
