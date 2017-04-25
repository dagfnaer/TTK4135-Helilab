%% Linear Quadratic Control

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

num_variables = 5/deltaT;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

t = 0:deltaT:deltaT*(length(u)-1);


%% LQR
% A and B matrix from 10.2

Q_lqr = diag([1,0,0,0]); % Try using same weights as in 10.2
R_lqr = diag([5]); % Constant 1
x = [x1;x2;x3;x4];
[K,S,e] = dlqr(A1,B1,Q_lqr,R_lqr);
x_var.time = t;
x_var.signals.values = [x1,x2,x3,x4];
x_var.signals.dimensions = 4;
u_var.time = t;
u_var.signals.values = u;
u_var.signals.dimensions = 1;