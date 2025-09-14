function [t, X] = simulate_fixed_vec(dyn, tspan, x0, dt, method)

    t = (tspan(1):dt:tspan(end))';
    n = numel(t);
    m = numel(x0);
    X = zeros(n, m); 
    X(1,:) = x0(:).';

    switch lower(method)
        case 'euler'
            for k = 1:n-1
                X(k+1,:) = euler_step(dyn, t(k), X(k,:).', dt).';
            end

        case 'rk4'
            for k = 1:n-1
                X(k+1,:) = rk4_step(dyn, t(k), X(k,:).', dt).';
            end

        otherwise
            error('simulate_fixed_vec: method must be euler or rk4');
    end
end