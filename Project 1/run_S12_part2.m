clear; clc;

J1 = 100;  b1 = 1;
J2 = 1;    b2 = 1;
Tf = 25;                           
tspan = [0 Tf];
A_list = [1, 100];                  
K_list = [10, 100, 1000];           
dt_list = [0.1, 1.0];             
w01 = 0; th01 = 0;                
w02 = 0; th02 = 0;                 


solver_for_plots = 'euler_dt0.1';           % we can change solvers to any of these ode45, euler_dt0.1, euler_dt1, rk4_dt0.1, rk4_dt1
                                            % some solver
                               


rows = {};
CPU  = [];

fprintf('Part 2 Simulations (Options 1–3) \n');

for Ai = 1:numel(A_list)
    A = A_list(Ai);
    tau_fun = @(t) tau_const(t, A);

    fprintf('\n-- Step Input A = %g --\n', A);

    % Flexible shaft
    for kk = 1:numel(K_list)
        k = K_list(kk);

        % Variable-step 
        tic;
        dyn1 = @(t,x) dyn_S12_flex(t, x, J1, b1, J2, b2, k, tau_fun);
        [t_o45, X_o45] = ode45(dyn1, tspan, [w01; th01; w02; th02]);
        cpu_o45 = toc;
        rows{end+1,1} = sprintf('Opt1 (k=%g) | ode45', k);  CPU(end+1,1) = cpu_o45; 

        % Euler/RK4 at dt = 0.1 and 1.0
        for dt = dt_list
            % Euler
            tic; [t_eu, X_eu] = simulate_fixed_vec(dyn1, tspan, [w01; th01; w02; th02], dt, 'euler'); cpu_eu = toc;
            rows{end+1,1} = sprintf('Opt1 (k=%g) | euler dt=%.1f', k, dt); CPU(end+1,1) = cpu_eu;

            % RK4
            tic; [t_rk, X_rk] = simulate_fixed_vec(dyn1, tspan, [w01; th01; w02; th02], dt, 'rk4');   cpu_rk = toc;
            rows{end+1,1} = sprintf('Opt1 (k=%g) | rk4 dt=%.1f', k, dt);   CPU(end+1,1) = cpu_rk;

            % PLot
            key_eu = sprintf('euler_dt%.1f', dt);
            key_rk = sprintf('rk4_dt%.1f',   dt);
            if strcmpi(solver_for_plots,'ode45')
                S.Opt1(kk).t = t_o45;  S.Opt1(kk).w1 = X_o45(:,1); S.Opt1(kk).w2 = X_o45(:,3);
            elseif strcmpi(solver_for_plots, key_eu)
                S.Opt1(kk).t = t_eu;   S.Opt1(kk).w1 = X_eu(:,1);  S.Opt1(kk).w2 = X_eu(:,3);
            elseif strcmpi(solver_for_plots, key_rk)
                S.Opt1(kk).t = t_rk;   S.Opt1(kk).w1 = X_rk(:,1);  S.Opt1(kk).w2 = X_rk(:,3);
            end
        end
    end

    % Lumped inertia/damping 
    Jtot = J1 + J2;  btot = b1 + b2;
    dyn2 = @(t,x) dyn_S12_combined(t, x, Jtot, btot, tau_fun);

    tic; [t2_o45, X2_o45] = ode45(dyn2, tspan, [w01; th01]); cpu_o45 = toc;
    rows{end+1,1} = 'Opt2 (lumped) | ode45'; CPU(end+1,1) = cpu_o45;

    for dt = dt_list
        tic; [t2_eu, X2_eu] = simulate_fixed_vec(dyn2, tspan, [w01; th01], dt, 'euler'); cpu_eu = toc;
        rows{end+1,1} = sprintf('Opt2 (lumped) | euler dt=%.1f', dt); CPU(end+1,1) = cpu_eu;

        tic; [t2_rk, X2_rk] = simulate_fixed_vec(dyn2, tspan, [w01; th01], dt, 'rk4');   cpu_rk = toc;
        rows{end+1,1} = sprintf('Opt2 (lumped) | rk4 dt=%.1f', dt);   CPU(end+1,1) = cpu_rk;
    end

    % Integration inside S2 
    tic; [t3_o45, X3_o45] = s2_integrated_opt3(w01, th01, tspan, J1, b1, J2, b2, tau_fun, 'ode45', NaN); cpu_o45 = toc;
    rows{end+1,1} = 'Opt3 (S2 integrates) | ode45'; CPU(end+1,1) = cpu_o45;

    for dt = dt_list
        tic; [t3_eu, X3_eu] = s2_integrated_opt3(w01, th01, tspan, J1, b1, J2, b2, tau_fun, 'euler', dt); cpu_eu = toc;
        rows{end+1,1} = sprintf('Opt3 (S2 integrates) | euler dt=%.1f', dt); CPU(end+1,1) = cpu_eu;

        tic; [t3_rk, X3_rk] = s2_integrated_opt3(w01, th01, tspan, J1, b1, J2, b2, tau_fun, 'rk4', dt);   cpu_rk = toc;
        rows{end+1,1} = sprintf('Opt3 (S2 integrates) | rk4 dt=%.1f', dt);   CPU(end+1,1) = cpu_rk;
    end

    %  Plot
   
    for kk = 1:numel(K_list)
        k = K_list(kk);
        figure; hold on; grid on;
        plot(S.Opt1(kk).t, S.Opt1(kk).w1, 'r-',  'LineWidth',1.2); 
        plot(S.Opt1(kk).t, S.Opt1(kk).w2, 'm--', 'LineWidth',1.2); 

       
        plot(t2_o45, X2_o45(:,1), 'b-.', 'LineWidth',1.5);          
        plot(t3_o45, X3_o45(:,1), 'g:',  'LineWidth',1.5);          

        xlabel('Time [s]'); ylabel('\omega [rad/s]');
        title(sprintf('A = %g  |  Compare Options @ k=%g  (%s)', A, k, solver_for_plots));
        legend('Opt1 \omega_1','Opt1 \omega_2','Opt2 (lumped)','Opt3 (S2 integrates)','Location','best');
    end
end

fprintf('\n CPU Time Summary \n');
for i = 1:numel(rows)
    fprintf('%-35s : %.6fs\n', rows{i}, CPU(i));
end
