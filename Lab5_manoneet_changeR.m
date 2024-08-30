%% CE-677:Lab 5 (EKF)

syms x y x1 y1 x2 y2 x3 y3 x4 y4 
% Mathematical Model
d1 = sqrt((x-x1)^2 + (y-y1)^2);
d2 = sqrt((x-x2)^2 + (y-y2)^2);
d3 = sqrt((x-x3)^2 + (y-y3)^2);
d4 = sqrt((x-x4)^2 + (y-y4)^2);

disp(d1)
disp(d2)
disp(d3)
disp(d4)

Pki = diag([1, 1, 1, 1]); 

% Qi processed noise covariance
Qki = diag([0.01,0.01,0.01,0.01]);

% Ri Observed noise covariance
Ri = diag([2, 2, 2 , 2]); 

% Given values of CP are
x1 = -10;
y1 = 0;
x2 = 0;
y2 = -10;
x3 = 10;
y3 = 0;
x4 = 0;
y4 = 10;

P_store = [];

Zkk = [
    9.72, 9.72, 22.08, 22.08;
    8.77, 8.77, 20.74, 20.74;
    7.97, 7.97, 19.42, 19.42;
    7.34, 7.34, 18.10, 18.10;
    6.93, 6.93, 16.81, 16.81;
    6.79, 6.79, 15.53, 15.53;
    6.93, 6.93, 14.28, 14.28;
    7.34, 7.34, 13.06, 13.06;
    7.97, 7.97, 11.88, 11.88;
    8.77, 8.77, 10.76, 10.76;
    9.72, 9.72, 9.72, 9.72;
    10.76, 10.76, 8.77, 8.77;
    11.88, 11.88, 7.97, 7.97;
    13.06, 13.06, 7.34, 7.34;
    14.28, 14.28, 6.93, 6.93;
    15.53, 15.53, 6.79, 6.79;
    16.81, 16.81, 6.93, 6.93;
    18.10, 18.10, 7.34, 7.34;
    19.42, 19.42, 7.97, 7.97;
    20.74, 20.74, 8.77, 8.77;
    22.08, 22.08, 9.72, 9.72
];

% Measurement Model
X_kp_kp = [-0.28415;-0.28415;0;0]; % solving two simultaneous eqn we initialize

x=[];
y=[];
% for Loop of Prediction and Estimation in KF
for i =1:21
    
    syms x y vx vy
    % for equations
    d10 = sqrt((x + 10)^2 + y^2);
    d20 = sqrt(x^2 + (y + 10)^2);
    d30 = sqrt((x - 10)^2 + y^2);
    d40 = sqrt(x^2 + (y - 10)^2);
    Hk = jacobian([d10, d20, d30, d40], [x,y, vx, vy]);
    Hik = double(subs(Hk, [x,y,vx,vy], [X_kp_kp(1,1), X_kp_kp(2,1), X_kp_kp(3,1), X_kp_kp(4,1)]));

    Fk = [1,0,1,0;0,1,0,1;0,0,1,0;0,0,0,1];

    % Now we predict the values
    Pk = (Fk*Pki*Fk') + Qki;
    xk = X_kp_kp(1,1) + X_kp_kp(3,1) * 1;
    yk = X_kp_kp(2,1) + X_kp_kp(4,1) * 1;
    vxk = X_kp_kp(3,1) ;
    vyk = X_kp_kp(4,1) ;
    X_kp=[xk;yk;vxk;vyk];

     d1_0 = double(subs(d10, [x, y], [X_kp_kp(1,1),X_kp_kp(2,1)]));
     d2_0 = double(subs(d20, [x, y], [X_kp_kp(1,1),X_kp_kp(2,1)]));
     d3_0 = double(subs(d30, [x, y], [X_kp_kp(1,1),X_kp_kp(2,1)]));
     d4_0 = double(subs(d40, [x, y], [X_kp_kp(1,1),X_kp_kp(2,1)]));

  f0 = [d1_0;d2_0;d3_0;d4_0];
  gain = Pk * Hik' * inv(Hik * Pk * Hik' + Ri);
  Zk = Hik * X_kp + f0 - Hik * X_kp_kp;

   X_k_k = X_kp + gain * ([Zkk(i,1);Zkk(i,2);Zkk(i,3);Zkk(i,4)] - Zk);

   P_k_k = (diag([1;1;1;1]) - gain * Hik) * Pk;    

    % Updation of values
    X_kp_kp = X_k_k;
    Pki = P_k_k;

    trace_Pkk = trace(P_k_k);
    P_store(i)= trace_Pkk;

    L1 = store_P;
    L2 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21];

    disp("Values of position =")
    disp(X_kp_kp)

    x=[-0.28415;X_kp_kp(1)];
    y=[-0.28415;X_kp_kp(2)];

end


% for Trajectory graph plot
figure;
plot(x,y,'r-', 'LineWidth', 1, 'MarkerSize', 5)
xlabel(x)
ylabel(y)
title("Trajectory path estimated by EKF change")
grid on

% for Trajectory graph plot
figure;
plot(L2,P_store,'r-', 'LineWidth', 1, 'MarkerSize', 5)
xlabel('Time in sec');
ylabel('Trace of P')
title("Trace of P vs Time")
grid on
