function xdot = dyn_S12_combined(t, x, Jtot, btot, tau_fun)
    omega  = x(1);
    omegad = (tau_fun(t) - btot*omega) / Jtot;
    xdot   = [omegad; omega];
end