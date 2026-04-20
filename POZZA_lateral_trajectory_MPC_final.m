%%      Control of Cyber-physical systems:                            |
%                                                                     |
%       Lateral Vehicle Trajectory Optimization                       |
%       using Constrained Linear Time-Varying MPC                     |
%       -----------------------------------------                     |
%       Michelangelo Pozza - IN2300012                                |
%----------------------------------------------------------------------
%%
clear all
close all
clc

% Horizon
N = 20;

% Time step = 200ms
Ts = 0.2;

% Default simulation time
T = 35;

% State vector
n = 5;

% Input vector
m = 1;

% Disturbance vector
o = 1;

% Output vector
p = 4;

%% Vehicle parameters (BMW i3: wheelbase l = 2.57m)
l = 2.57;
r_circles = 1; % because r_circle < l/2 = 1.29
u_max = 0.25; % [1/ms]
u_min = -u_max;
k_max_delta = 0.25; % [1/m]
k_min_delta = -k_max_delta;
g = 9.81;
mu = 0.9;

%% User's menu

disp('------------------------------------------------------------------');
disp('      LATERAL VEHICLE TRAJECTORY OPTIMIZATION USING LMPC          ');
disp('------------------------------------------------------------------');
disp('Select the test you want to execute:');
disp('  1: Straight road with an obstacle');
disp('  2: Sine() road without obstacles');
disp('  3: Roundabout');
disp('  4: Parking scenario with more obstacles');
disp('  5: Lemniscate');
disp('  6: Simple curve');
disp('------------------------------------------------------------------');
number_test = input('Type the number of the test you want to execute: ');

%% Trajectory to follow
ds = 0.05;

switch number_test
    case 1
        T = 30;
        % waypoints of the track
        x_plot = 0:ds:40; 
        y_plot = zeros(1, length(x_plot));

        % Transformation between points (x,y) to (k,theta)
        % k = |x'y''-y'x''|/(x'^2+y'^2)^(3/2)
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);

        ddx = diff(dx); 
        ddy = diff(dy);
        
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(1.5);
        
        kappa_r_vec = num ./ den;
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        % computation of theta_r -> using [dy dy(end)] because diff() function
        % provide the intervals which are N_points-1 wrt original points. We need
        % to add the last piece as constat to have the same number of points
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        dtheta = diff(theta_r_vec);
        
        % kappa_r = ||dtheta|| ==> normalize using ds
        % other way: kappa_r_vec = dtheta./ds_step;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)]; % padding because of the differences
        
        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        % Casual obstacle
        obs_start = find(s_ref >= 7.5, 1);
        obs_end   = find(s_ref >= 9, 1);
        
        % Modification of the bounds to add the obstacle
        d_min_vec(obs_start:obs_end) =d_min_vec(obs_start:obs_end) + 1.3;
        
        % straight road but with an obstacle => not so fast
        v_max = 8;
        
        v_spatial = v_max*ones(1, length(kappa_r_vec));
    case 2
        T = 50;
        x_plot = 0:ds:160; 
        
        % amplitude and frequency of the sine
        amplitude = 20;
        lambda_sine = 80;
        y_plot = amplitude*sin((2*pi/lambda_sine)*x_plot);

        % Transformation between points (x,y) to (k,theta)
        % k = |x'y''-y'x''|/(x'^2+y'^2)^(3/2)
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);

        ddx = diff(dx); 
        ddy = diff(dy);
        
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(1.5);
        
        kappa_r_vec = num ./ den;
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        % computation of theta_r -> using [dy dy(end)] because diff() function
        % provide the intervals which are N_points-1 wrt original points. We need
        % to add the last piece as constat to have the same number of points
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        dtheta = diff(theta_r_vec);
        
        % kappa_r = ||dtheta|| ==> normalize using ds
        % other way: kappa_r_vec = dtheta./ds_step;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)]; % padding because of the differences
        
        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        % No obstacles => If possible go fast
        v_max = 10;
        
        v_spatial = v_max*ones(1, length(kappa_r_vec));
    case 3
        T = 30;

        % let's make the roundabout as if it was a gaussian
        round_mu = 30;
        round_alpha = 10;
        round_sigma = 10;
        
        x_plot = 0:ds:60;
        y_plot = -round_alpha*exp(-(x_plot-round_mu).^2/(2*round_sigma^2));
        
        % Transformation between points (x,y) to (k,theta)
        % k = |x'y''-y'x''|/(x'^2+y'^2)^(3/2)
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);

        ddx = diff(dx); 
        ddy = diff(dy);
        
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(1.5);
        
        kappa_r_vec = num ./ den;
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        % computation of theta_r -> using [dy dy(end)] because diff() function
        % provide the intervals which are N_points-1 wrt original points. We need
        % to add the last piece as constat to have the same number of points
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        dtheta = diff(theta_r_vec);
        
        % kappa_r = ||dtheta|| ==> normalize using ds
        % other way: kappa_r_vec = dtheta./ds_step;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)]; % padding because of the differences
        
        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        % No obstacles => If possible go fast
        v_max = 10;
        
        v_spatial = v_max*ones(1, length(kappa_r_vec));
    case 4
        % x1 and y1 represent the parking
        x1 = 0:ds:5; 
        y1 = zeros(1, length(x1));
        
        % 5m radius circle to describe the curve after the parking
        r_curve = 5.0; 
        
        % ds = r_circle*dl_curve => dl_curve is a little piece of the arc length
        dl_curve = ds/r_curve; 
        curve_angles_set = (-pi/2+dl_curve):dl_curve:0; 
        x2 = x1(end)+r_curve*cos(curve_angles_set);
        y2 = r_curve+r_curve*sin(curve_angles_set);
        
        % last straight part
        y3 = (r_curve+ds):ds:15;
        x3 = (x1(end)+r_curve)*ones(1, length(y3));
        
        x_plot = [x1, x2, x3];
        y_plot = [y1, y2, y3];
        
        % Transformation between points (x,y) to (k,theta)
        % k = |x'y''-y'x''|/(x'^2+y'^2)^(3/2)
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);

        ddx = diff(dx); 
        ddy = diff(dy);
        
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(1.5);
        
        kappa_r_vec = num ./ den;
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        % computation of theta_r -> using [dy dy(end)] because diff() function
        % provide the intervals which are N_points-1 wrt original points. We need
        % to add the last piece as constat to have the same number of points
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        dtheta = diff(theta_r_vec);
        
        % kappa_r = ||dtheta|| ==> normalize using ds
        % other way: kappa_r_vec = dtheta./ds_step;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)]; % padding because of the differences
        
        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        obs_start = find(s_ref >= 5.0, 1);
        obs_end   = find(s_ref >= 8.5, 1);
        d_max_vec(obs_start:obs_end) = d_max_vec(obs_start:obs_end) - 1.2;
        
        obs2_start = find(s_ref >= 10, 1);
        obs2_end   = find(s_ref >= 12, 1);
        d_min_vec(obs2_start:obs2_end) = d_min_vec(obs2_start:obs2_end) + 0.8;
        
        obs3_start = find(s_ref >= 14, 1);
        obs3_end   = find(s_ref >= 15.5, 1);
        d_min_vec(obs3_start:obs3_end) = d_min_vec(obs3_start:obs3_end) + 1.2;
        
        % parking situation, we may not going too fast
        v_max = 5;
        
        v_spatial = v_max*ones(1, length(kappa_r_vec));
    case 5
        T = 100;
        a = 60;
        t_param = -pi:ds:pi; 
        x_plot = (a*cos(t_param))./(1+sin(t_param).^2);
        y_plot = (a*sin(t_param).*cos(t_param))./(1+sin(t_param).^2);
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);
        ddx = diff(dx); 
        ddy = diff(dy);
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(3/2);
        kappa_r_vec = num ./ den;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)];

        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        v_max = 8;

        v_spatial = v_max*ones(1, length(kappa_r_vec));
    case 6
        T = 40;

        N_curve = 150;

        theta = linspace(0, pi/2, N_curve);
        r = 60;
        x_plot = r*cos(theta);
        y_plot = r*sin(theta);

        % Transformation between points (x,y) to (k,theta)
        % k = |x'y''-y'x''|/(x'^2+y'^2)^(3/2)
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);

        ddx = diff(dx); 
        ddy = diff(dy);
        
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(1.5);
        
        kappa_r_vec = num ./ den;
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        % computation of theta_r -> using [dy dy(end)] because diff() function
        % provide the intervals which are N_points-1 wrt original points. We need
        % to add the last piece as constat to have the same number of points
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        dtheta = diff(theta_r_vec);
        
        % kappa_r = ||dtheta|| ==> normalize using ds
        % other way: kappa_r_vec = dtheta./ds_step;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)]; % padding because of the differences
        
        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        v_max = 10;
        
        v_spatial = v_max*ones(1, length(kappa_r_vec));
    otherwise
        disp("The test you selected does not exist! Let's execute test 1");
        T = 30;
        % waypoints of the track
        x_plot = 0:ds:40; 
        y_plot = zeros(1, length(x_plot));

        % Transformation between points (x,y) to (k,theta)
        % k = |x'y''-y'x''|/(x'^2+y'^2)^(3/2)
        
        % dx=x' and dy=y' (differences at each step)
        dx = diff(x_plot);
        dy = diff(y_plot);

        ddx = diff(dx); 
        ddy = diff(dy);
        
        num = dx(1:end-1).*ddy - dy(1:end-1).*ddx;
        den = (dx(1:end-1).^2 + dy(1:end-1).^2).^(1.5);
        
        kappa_r_vec = num ./ den;
        
        % parametrization on the arc of the curves
        ds_step = sqrt(dx.^2 + dy.^2); % this is the derivative 
        s_ref = [0, cumsum(ds_step)];
        
        % computation of theta_r -> using [dy dy(end)] because diff() function
        % provide the intervals which are N_points-1 wrt original points. We need
        % to add the last piece as constat to have the same number of points
        theta_r_vec = unwrap(atan2([dy dy(end)], [dx dx(end)])); 
        dtheta = diff(theta_r_vec);
        
        % kappa_r = ||dtheta|| ==> normalize using ds
        % other way: kappa_r_vec = dtheta./ds;
        kappa_r_vec = [kappa_r_vec, kappa_r_vec(end), kappa_r_vec(end)]; % padding because of the differences
        
        N_curve = length(x_plot);
        
        % Lateral limits because of the road
        d_max_vec =  (2-r_circles)*ones(1, N_curve);
        d_min_vec = (-2+r_circles)*ones(1, N_curve);
        
        % Casual obstacle
        obs_start = find(s_ref >= 7.5, 1);
        obs_end   = find(s_ref >= 9, 1);
        
        % Modification of the bounds to add the obstacle
        d_min_vec(obs_start:obs_end) =d_min_vec(obs_start:obs_end) + 1.3;
        
        % straight road but with an obstacle => not so fast
        v_max = 8;
        
        v_spatial = v_max*ones(1, length(kappa_r_vec));
end

%% MPC

% vectors to save the path and the velocity profile of the vehicle
v_history = zeros(1, T+1);
s_history = zeros(1, T+1);
k_min_history = zeros(1, T+1);
k_max_history = zeros(1, T+1);

v_curr=v_spatial(1);

x0 = [0;
      theta_r_vec(1)-deg2rad(5);
      0;
      theta_r_vec(1);
      kappa_r_vec(1)];

x=x0;
u=[];

for time_instant = 1:T

    if time_instant == 1
        s_curr = 0;
    else
        s_curr = s_curr + v_curr * Ts; 
    end

    % If computing the next optimization I don't have any other data stop
    % the simulation
    if (s_curr >= s_ref(end)-N*ds(end))
        disp('Finished!');
        T = time_instant-1;
        break;
    end

    s_history(time_instant) = s_curr;
    
    % Prediction of position on the curve and velocity online
    s_pred = zeros(N+1, 1);
    v_pred = zeros(N+1, 1);
    
    s_pred(1) = s_curr;
    v_pred(1) = interp1(s_ref, v_spatial, s_pred(1), 'linear', 'extrap');
    
    for i = 1:N
        s_pred(i+1) = s_pred(i) + v_pred(i) * Ts;
        v_pred(i+1) = interp1(s_ref, v_spatial, s_pred(i+1), 'linear', 'extrap');
    end
    
    v_curr = v_pred(1);
    v_history(time_instant) = v_curr;
    
    % linear interpolation because of the sampling time of the vehicle wrt
    % the actual curve 
    kappa_fut = interp1(s_ref, kappa_r_vec, s_pred, 'linear', 'extrap');
    theta_fut = interp1(s_ref, theta_r_vec, s_pred, 'linear', 'extrap');
    % ZOH interpolation because the constraints must be exactly the same as
    % the curve. PROBLEM: If the constraint is "too small" wrt the sampling
    % time MPC will not see it (just put obstacles with length at least 1)
    d_max_fut = interp1(s_ref, d_max_vec, s_pred, 'previous', 'extrap');
    d_min_fut = interp1(s_ref, d_min_vec, s_pred, 'previous', 'extrap');
    
    % Definition of z(k)=(k_r(k+1)-k_r(k))/Ts
    z = zeros(N, 1);
    for k = 1:N
        z(k) = (kappa_fut(k+1) - kappa_fut(k)) / Ts;
    end

    % Continuous-time model    
    Ac = zeros(n, n, N+1);
    Bc = zeros(n, m, N+1);
    Ec = zeros(n, o, N+1);
    Cc = zeros(p, n, N+1);
    v = v_pred;
    for i = 1:N+1
        Ac(:,:,i) = [0      v(i)    0       -v(i)   0;
                     0      0       v(i)    0       0;
                     0      0       0       0       0;
                     0      0       0       0       v(i);
                     0      0       0       0       0];
    
        Bc(:,:,i) = [0;
                     0;
                     1;
                     0;
                     0];
    
        Ec(:,:,i) = [0;
                     0;
                     0;
                     0;
                     1];
    
        Cc(:,:,i) = [1      0       0       0       0;
                     1      l/2     0       -l/2    0;
                     1      l       0       -l      0;
                     0      0       1       0       0];
    end
    
    % Discrete-time model
    Ad = zeros(n, n, N+1);
    Bd = zeros(n, m, N+1);
    Ed = zeros(n, o, N+1);
    Cd = zeros(p, n, N+1);
    for i = 1:N+1
        A = Ac(:,:,i);
        B = Bc(:,:,i);
        C = Cc(:,:,i);
        E = Ec(:,:,i);
    
        sys_c = ss(A, B, C, 0);
        sys_d = c2d(sys_c, Ts);
    
        Ad(:,:,i) = sys_d.A;
        Bd(:,:,i) = sys_d.B;
        Cd(:,:,i) = sys_d.C;
    
        sys_c = ss(A, E, C, 0);
        sys_d = c2d(sys_c, Ts);
    
        Ed(:,:,i) = sys_d.B;
    end
    
    
    
    % Optimal control problem - cost function
    w_theta = 10;
    w_u = 10;
    
    Q = zeros(n, n, N);
    R = zeros(1, N);
    for i = 1:N+1
        w_d = 2/v_pred(i);
        w_k = v_pred(i)/10;
        Q(:,:,i) = [w_d          0              0           0               0;
                    0            w_theta        0           -w_theta        0;
                    0            0              w_k         0               0;
                    0            -w_theta       0           w_theta         0;
                    0            0              0           0               0];
        
        R(i) = w_u;
    end
    
    A_big = zeros(n*(N+1), n);
    A_big(1:n,:) = eye(n);
    A_prod = Ad(:,:,1);
    for i = 1:N
        idx_n = i*n;
        A_big(idx_n+1:idx_n+n,:) = A_prod;
        A_prod = Ad(:,:,i+1)*A_prod;
    end
    
    
    B_big = zeros(n*(N+1),m*N);
    for j = 1:N
        B_prod = Bd(:,:,j);
        k = (j-1)*m;
        for i = j:N
            idx_n = i*n;
            B_big(idx_n+1:idx_n+n,k+1:k+m) = B_prod;
            B_prod = Ad(:,:,i+1) * B_prod;
        end
    end
    
    C_big = zeros(p*(N+1), n*(N+1));
    for i = 1:N+1
        idx_n = (i-1)*n;
        idx_p = (i-1)*p;
        C_big(idx_p+1 : idx_p+p, idx_n+1 : idx_n+n) = Cd(:,:,i);
    end
    
    Eps_big = zeros(n*(N+1),o*N);
    for j = 1:N
        Eps_prod = Ed(:,:,j);
        k = (j-1)*o;
        for i = j:N
            idx_n = i*n;
            Eps_big(idx_n+1:idx_n+n,k+1:k+o) = Eps_prod;
            Eps_prod = Ad(:,:,i+1)*Eps_prod;
        end
    end
    
    % Now we insert x(N) on the Q-matrix as P
    w_dP = 1000;
    w_thetaP = 1000;
    w_kP = 1000;
    P =    [w_dP            0               0           0               0;
            0               w_thetaP        0           -w_thetaP       0;
            0               0               w_kP        0               0;
            0               -w_thetaP       0           w_thetaP        0;
            0               0               0           0               0];
    
    Q_big = zeros(n*(N+1), n*(N+1));
    for i = 1:N
        idx_n = (i-1)*n;
        Q_big(idx_n+1:idx_n+n,idx_n+1:idx_n+n) = Q(:,:,i);
    end
    idx_n = N*n;
    Q_big(idx_n+1:idx_n+n,idx_n+1:idx_n+n) = P;
    
    R_big = zeros(m*N, m*N);
    for i = 1:N
        idx_m = (i-1)*m;
        R_big(idx_m+1:idx_m+m,idx_m+1:idx_m+m) = R(i);
    end
    
    H_big = B_big'*Q_big*B_big + R_big;
    
    
    
    F_big = A_big'*Q_big*B_big;
    G_big = Eps_big'*Q_big*B_big;
    
    f = 2*(x0'*F_big+z'*G_big);
    f = f';
    
    % Softer problem
    Ns = 5;
    k_1 = [1e3;
           1e3];
    
    K_2 = diag([1e5; 1e5]);
    
    % now u_tilde = [u', [eps_u, eps_l]']'
    p_s = 3;
    p_h = p;
    
    H_tilde = blkdiag(H_big, K_2);
    F_tilde = [f; k_1];
    
    % Modify output sequence 
    % y_s are the soft outputs
    % y_h are the hard outputs
    Ch_big = [];
    C_h_Ks = zeros(1, n, Ns);
    for i = 1:Ns
        C_h_Ks(:,:,i) = [0, 0, 0, 1]*Cd(:,:,i);
        Ch_big = blkdiag(Ch_big, C_h_Ks(:,:,i));
    end
    C_h = zeros(p, n, N+1-Ns);
    for i = 1:N+1-Ns
        C_h(:,:,i) = Cd(:,:,i+Ns);
        Ch_big = blkdiag(Ch_big, C_h(:,:,i));
    end
    
    C_s = zeros(p_s, n, Ns);
    Cs_big = [];
    for i = 1:Ns
        C_s(:,:,i) = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0]*Cd(:,:,i);
        Cs_big = blkdiag(Cs_big, C_s(:,:,i));
    end
    
    
    % Constraints for the output
    k_max = zeros(N+1, 1);
    k_min = zeros(N+1, 1);
    
    for i = 1:N+1
        v = v_pred(i);

        % Kamm circle for the tire maximum grip that influences the max
        % steering action
        k_max_mu = (mu*g)/(v^2); 
        
        % close the bounds around the minimum and maximum steering actions
        k_max(i) = min(k_max_delta, k_max_mu);
        k_min(i) = max(k_min_delta, -k_max_mu);
    end

    k_min_history(time_instant) = k_min(1);
    k_max_history(time_instant) = k_max(1);
    
    s = 3;
    y_s_max = zeros(s*Ns, 1); y_s_min = zeros(s*Ns, 1);
    for i = 1:Ns
        is = (i-1)*s;
        y_s_max(is+1:is+s) = [d_max_fut(i); d_max_fut(i); d_max_fut(i)];
        y_s_min(is+1:is+s) = [d_min_fut(i); d_min_fut(i); d_min_fut(i)];
    end

    y_h_max_Ks = k_max(1:Ns);
    y_h_min_Ks = k_min(1:Ns);
    
    y_h_max_k = zeros(p*(N+1-Ns), 1);
    y_h_min_k = zeros(p*(N+1-Ns), 1);
    for i = 1:(N+1-Ns)
        ip = (i-1)*p;
        y_h_max_k(ip+1 : ip+p) = [d_max_fut(Ns+i); d_max_fut(Ns+i); d_max_fut(Ns+i); k_max(Ns+i)];
        y_h_min_k(ip+1 : ip+p) = [d_min_fut(Ns+i); d_min_fut(Ns+i); d_min_fut(Ns+i); k_min(Ns+i)];
    end
    
    y_h_max = [y_h_max_Ks; y_h_max_k];
    y_h_min = [y_h_min_Ks; y_h_min_k];
    
    % Final constraints
    ps = 3*Ns;
    ph = p*(N+1)-ps;
    r = 2;

    Cs_big = [Cs_big, zeros(ps, n*(N+1-Ns))];
    
    M_11 = [eye(N); -eye(N)];
    M_21 = [Ch_big*B_big; -Ch_big*B_big];
    M_12 = zeros(2*(N), r); 
    M_22 = zeros(2*ph, r); 
    M_31 = [Cs_big*B_big; -Cs_big*B_big];
    M_41 = zeros(r, N);
    M_32 = [-ones(ps, 1), zeros(ps, 1); zeros(ps, 1), -ones(ps, 1)];
    M_42 = -eye(r);
    
    Ac_tilde = [M_11, M_12;
                  M_21, M_22;
                  M_31, M_32;
                  M_41, M_42];
    
    V_1 = [u_max*ones(N, 1);
           -u_min*ones(N,1)];
    
    V_2 = [y_h_max-Ch_big*(A_big*x0+Eps_big*z);
           -y_h_min+Ch_big*(A_big*x0+Eps_big*z)];
    
    V_3 = [y_s_max-Cs_big*(A_big*x0+Eps_big*z);
           -y_s_min+Cs_big*(A_big*x0+Eps_big*z)];
    
    V_4 = zeros(r, 1);
    
    bc_tilde = [V_1; V_2; V_3; V_4];
    
    % I recall the function quadprog
    % UMPC contains the optimal input sequence
    % (H+H')/2 symmetrized a matrix in case of small numerical errors
    [UMPC,FVAL,EXITFLAG]=quadprog((H_tilde+H_tilde')/2,F_tilde,Ac_tilde,bc_tilde,[],[],[],[],[]);
    
    % exitflag tells what was the result of the optimization
    % if -2 the problem was not feasible
    if EXITFLAG==-2
        disp('Something wrong with the optimization')
        time_instant
        keyboard
    else
        u_curr=UMPC(1);
        u=[u u_curr];
        z_curr = z(1);
        
        [~, x_sim] = ode45(@(t, x) solution(t, x, u_curr, z_curr, v_curr), [time_instant*Ts, (time_instant+1)*Ts], x(:, end));
        
        xnext = x_sim(end, :)'; 
        x=[x xnext];
        x0 = x(:, end);
    end

end

% Because the simulation can stop before the original T, we may modified it
% due to the next horizon longer then the given data
% we shrink the vectors
s_history = s_history(1:T);
v_history = v_history(1:T);
k_min_history = k_min_history(1:T);
k_max_history = k_max_history(1:T);

% Let's add the last step to the space and velocity vectors
s_history(T+1) = s_history(T)+v_history(T)*Ts;
v_history(T+1) = v_history(T);
k_min_history(T+1) = k_min_history(T);
k_max_history(T+1) = k_max_history(T);
u = [u 0];

%% Plot of the results
d_r_plot = x(1, :); 
theta_plot = x(2, :);
kappa_plot = x(3, :);
theta_r_plot = x(4, :);

x_r_interp = interp1(s_ref, x_plot, s_history, 'linear', 'extrap');
y_r_interp = interp1(s_ref, y_plot, s_history, 'linear', 'extrap');

x_vehicle = x_r_interp - d_r_plot .* sin(theta_r_plot);
y_vehicle = y_r_interp + d_r_plot .* cos(theta_r_plot);

d_max_sense = interp1(s_ref, d_max_vec, s_history, 'linear', 'extrap');
d_min_sense = interp1(s_ref, d_min_vec, s_history, 'linear', 'extrap');

x_u = x_plot-d_max_vec.*sin(theta_r_vec); 
y_u = y_plot+d_max_vec.*cos(theta_r_vec);
x_l = x_plot-d_min_vec.*sin(theta_r_vec); 
y_l = y_plot+d_min_vec.*cos(theta_r_vec);

figure('Name', 'Bird''s eye view');
hold on; 
grid on; 
axis equal;

plot(x_plot, y_plot, 'k--', 'LineWidth', 1.5, 'DisplayName', '\Gamma');
plot(x_u, y_u, 'r-', 'LineWidth', 1.5, 'DisplayName', 'd_{max}');
plot(x_l, y_l, 'g-', 'LineWidth', 1.5, 'DisplayName', 'd_{min}');

plot(x_vehicle, y_vehicle, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Vehicle');

plot(x_vehicle(1), y_vehicle(1), 'bo', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(x_vehicle(end), y_vehicle(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'End');

title('Bird''s eye view', 'FontSize', 14);
xlabel('X [m]', 'FontSize', 12); ylabel('Y [m]', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11);


%% Plot of the time analysis
T_sim = 0:Ts:T*Ts;

d_r_des = zeros(1, length(T_sim));
kappa_des = interp1(s_ref, kappa_r_vec, s_history, 'linear', 'extrap');
u_des = zeros(1, length(T_sim));

figure('Name', 'Time analysis', 'Color', 'w');

subplot(3, 2, 1:2);
stairs(T_sim, u_des, 'k--', 'LineWidth', 1.5, 'DisplayName', 'u_{n}');
hold on; 
grid on;
stairs(T_sim, u, 'b-', 'LineWidth', 2, 'DisplayName', 'u'); 
plot(T_sim, u_max*ones(1, length(T_sim)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'u_{max}');
plot(T_sim, u_min*ones(1, length(T_sim)), 'g--', 'LineWidth', 1.5, 'DisplayName', 'u_{min}');

title('Applied control input', 'FontSize', 12); 
ylabel('u [1/(m\cdot s)]', 'FontSize', 11);
legend('Location', 'northeast'); 
axis tight
set(gca, 'FontSize', 10);

subplot(3, 2, 3:4);
plot(T_sim, d_r_des, 'k--', 'LineWidth', 1.5, 'DisplayName', 'd_{rn}');
hold on; 
grid on;
plot(T_sim, d_r_plot, 'b-', 'LineWidth', 2, 'DisplayName', 'd_r'); 
stairs(T_sim, d_max_sense, 'r--', 'LineWidth', 1.5, 'DisplayName', 'd_{max}');
stairs(T_sim, d_min_sense, 'g--', 'LineWidth', 1.5, 'DisplayName', 'd_{min}');

title('Distance from Reference Trajectory', 'FontSize', 12); 
ylabel('d_r [m]', 'FontSize', 11);
legend('Location', 'northeast');
axis tight
set(gca, 'FontSize', 10);

subplot(3, 2, 5);
plot(T_sim, kappa_des, 'k--', 'LineWidth', 1.5, 'DisplayName', '\kappa_{n}');
hold on; 
grid on;
plot(T_sim, kappa_plot, 'b-', 'LineWidth', 2, 'DisplayName', '\kappa'); 
plot(T_sim, k_max_history, 'r--', 'LineWidth', 1.5, 'DisplayName', '\kappa_{max}');
plot(T_sim, k_min_history, 'g--', 'LineWidth', 1.5, 'DisplayName', '\kappa_{min}');

title('Vehicle curvature and friction limits', 'FontSize', 12); 
ylabel('\kappa [1/m]', 'FontSize', 11);
legend('Location', 'northeast'); 
xlabel('Time (s)', 'FontSize', 11); 
axis tight
set(gca, 'FontSize', 10);

subplot(3, 2, 6);
plot(T_sim, zeros(1, length(T_sim)), 'k--', 'LineWidth', 1.5, 'DisplayName', '(\theta-\theta_r)_n');
hold on; 
grid on;
plot(T_sim, rad2deg(theta_plot-theta_r_plot), 'b-', 'LineWidth', 2, 'DisplayName', '(\theta-\theta_r)'); 
plot(T_sim, 20*ones(1,length(T_sim)), 'r--', 'LineWidth', 1.5, 'DisplayName', '(\theta-\theta_r)_{max}');
plot(T_sim, -20*ones(1,length(T_sim)), 'g--', 'LineWidth', 1.5, 'DisplayName', '(\theta-\theta_r)_{min}');

title('Heading Angle Error', 'FontSize', 12);
ylabel('\theta - \theta_r [deg]', 'FontSize', 11);
legend('Location', 'northeast'); 
xlabel('Time (s)', 'FontSize', 11); 
axis tight
set(gca, 'FontSize', 10);

%% Plot of lateral deviation and velocity
figure('Name', 'Lateral deviation and velocity');

subplot(2, 1, 1);
plot(T_sim, d_r_des, 'k--', 'LineWidth', 1.5, 'DisplayName', 'd_{rn}');
hold on; 
grid on;
plot(T_sim, d_r_plot, 'b-', 'LineWidth', 2, 'DisplayName', 'd_r'); 
stairs(T_sim, d_max_sense, 'r--', 'LineWidth', 1.5, 'DisplayName', 'd_{max}');
stairs(T_sim, d_min_sense, 'g--', 'LineWidth', 1.5, 'DisplayName', 'd_{min}');
title('Lateral Deviation', 'FontSize', 12); 
ylabel('d_r [m]', 'FontSize', 11);
legend('Location', 'northeast');
axis tight;
set(gca, 'FontSize', 10);

subplot(2, 1, 2);
plot(T_sim, v_history, 'b-', 'LineWidth', 2, 'DisplayName', 'v');
hold on; 
grid on;
plot(T_sim, 20 * ones(1, length(T_sim)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'v_{max}');
plot(T_sim, zeros(1, length(T_sim)), 'g--', 'LineWidth', 1.5, 'DisplayName', 'v_{min}');
title('Vehicle velocity', 'FontSize', 12); 
xlabel('Time (s)', 'FontSize', 11); 
ylabel('v [m/s]', 'FontSize', 11);
legend('Location', 'northeast'); 
axis tight;
set(gca, 'FontSize', 10);

%% Reference curve vs Waypoints and interpolation
figure('Name', 'Reference curve vs Waypoints and interpolation');
plot(x_plot, y_plot, 'k-', 'DisplayName', 'Curve'); 
hold on;
plot(x_r_interp, y_r_interp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Interpolated curve');
title('Reference curve vs Waypoints and interpolation');
xlabel('X [m]'); 
ylabel('Y [m]');
legend; 
axis equal; 
grid on;

%% Functions

function dxdt = solution(t, x, u, z, v)
    dxdt = zeros(5, 1);

    dxdt(1) = v*sin(x(2)-x(4));
    dxdt(2) = v*x(3);
    dxdt(3) = u;
    dxdt(4) = (v*cos(x(2)-x(4)))/(1-x(1)*x(5))*x(5);
    dxdt(5) = z;
end
