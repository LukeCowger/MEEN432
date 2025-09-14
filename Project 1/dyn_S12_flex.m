function xdot = dyn_S12_flex(t, x, J1, b1, J2, b2, k, tau_fun)

    omega1 = x(1);  theta1 = x(2);
    omega2 = x(3);  theta2 = x(4);

    tau_k   = k * (theta1 - theta2);
    omega1d = (tau_fun(t) - b1*omega1 - tau_k) / J1;
    omega2d = (tau_k - b2*omega2) / J2;

    xdot = [omega1d; omega1; omega2d; omega2];
end