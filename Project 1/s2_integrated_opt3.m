function [t_out, x_out] = s2_integrated_opt3(w0, th0, tspan, J1, b1, J2, b2, tau_fun, method, dt)

    Jtot = J1 + J2;
    btot = b1 + b2;
    dyn  = @(t,x) dyn_S12_combined(t, x, Jtot, btot, tau_fun);

    x0 = [w0; th0];

    switch lower(method)
        case 'ode45'
            [t_out, x_out] = ode45(dyn, tspan, x0);

        case 'euler'
            [t_out, x_out] = simulate_fixed_vec(dyn, tspan, x0, dt, 'euler');

        case 'rk4'
            [t_out, x_out] = simulate_fixed_vec(dyn, tspan, x0, dt, 'rk4');

        otherwise
            error('Unknown method: %s', method);
    end
end