%% Optimal Control of Pitch/Travel and Elevation with and without Feedback
% 10.4 .....
% Initalizing...

init;

% A and B matrix state space model. x = [lamda, r, p, p_dot,e,e_dot]
deltaT = 0.25;

A2 = [0 1 zeros(1,4);
    0 0 -K_2 zeros(1,3);
    0 0 0 1 0 0;
    0 0 -K_1*K_pp -K_1*K_pd 0 0;
    zeros(1,5) 1;
    zeros(1,4) -K_3*K_ep -K_3*K_ed];

A2 = eye(6) + deltaT*A2; % Discretized A matrix

B2 = [0 0 0 deltaT*K_1*K_pp 0 0;
    zeros(1,5) deltaT*K_3*K_ep]'; % Discretized B matrix

% Number of states and inputs
mx = size(A2,2); % Number of states (number of columns in A)
mu = size(B2,2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';          % Initial values

% Time horizon and initialization
global N;
N = 40;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Bounds
ul 	    = [-30*pi/180;-Inf];                   % Lower bound on control -- u1
uu 	    = [30*pi/180;Inf];                    % Upper bound on control -- u1

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = -30*pi/180;                           % Lower bound on state x3
xu(3)   = 30*pi/180;                           % Upper bound on state x3

%% Optional constraints on velocity and elevation rate
xl(2) = -5*pi/180;
xu(2) = 5*pi/180;

xl(4)   = -5*pi/180;  
xu(4)   = 5*pi/180; 
%
xl(6) = -1*pi/180;
xu(6) = 10*pi/180;


% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N,M,xl,xu,ul,uu); 
vlb(end-1:end)  = [0;0];                    % We want the last input to be zero
vub(end-1:end)  = [0;0];                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q2 = zeros(mx,mx);
Q2(1,1) = 1;                             % Weight on state x1 travel
Q2(2,2) = 0;                            % Weight on state x2 velocity
Q2(3,3) = 0;                             % Weight on state x3 pitch
Q2(4,4) = 0;                            % Weight on state x4 pitch rate
Q2(5,5) = 0;                            % Weight on state x5 elevation
Q2(6,6) = 0;                            % Weight on state x6 elevation rate

P = diag([1,1]);
Q2 = genq2(Q2,P,N,M,mu);              % Generate Q
c = zeros(N*mx+M*mu,1);                 % Generate c

%% Generate system matrixes for linear model
Aeq = gena2(A2,B2,N,mx,mu);           % Generate A, hint: gena2
beq = zeros(N*mx,1);        	  % Generate b
beq(1:mx) = A2*x0; % Initial value
z0 = [x0;zeros(N*(mx+mu)-mx,1)]; % Use initial from x0
%% Calculate Elevation

fmin = @(z) z'*Q2*z;
nonlcon = @nonlconst;
options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',100000);
[z,fval_el,exitflag,output,lambda,grad,hessian] = fmincon(fmin,z0,[],[],Aeq,beq,vlb,vub,nonlcon,options);

%% Extract control inputs and states
pc  = [z(N*mx+1:2:N*mx+M*mu);z(N*mx+M*mu-1)]; % Control input pitch from solution
ec  = [z(N*mx+2:2:N*mx+M*mu);z(N*mx+M*mu)]; % Control input elevation from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution    

num_variables = 5/deltaT;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

pc  = [zero_padding; pc; zero_padding];
ec  = [zero_padding; ec; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

t = 0:deltaT:deltaT*(length(pc)-1);

%% LQR
% A and B matrix from 10.2

Q_lqr = diag([1,1,1,1,1,1]); % Q matrix lqr
R_lqr = diag([1,1]); % R matrix lqr
x = [x1;x2;x3;x4;x5;x6];
[K,S,e] = dlqr(A2,B2,Q_lqr,R_lqr);
x_var.time = t;
x_var.signals.values = [x1,x2,x3,x4,x5,x6];
x_var.signals.dimensions = 6;
u_var.time = t;
u_var.signals.values = [pc,ec];
u_var.signals.dimensions = 2;