% A script to produce approximate the minimum of a quadratic (zero-offset) function F(x) = x'g + 1/2 x'Jx given a set of constraints Ax = b (where n > m) and g_0(x)...g_k(x) <= 0

% for now, inequality constraint will be x_i <= 0 for all x, until the program is stable
% the specific values are from an example from Strang's textbook

% parts of function to minimize
g = [5;3;8];
J = zeros(3,3);

% constraints
A = [1 1 2];
b = [4];

% relevant dimensions
n = size(J, 1);
m = size(A, 1);

% slack parameter
theta = 4.0/3.0;
theta_lower_bound = 0.01;

% starting point
x = ones(3,1); % iterative x variable
y = 2; % iterative lagrange multiplier

% auxiliary vectors and matrices
B = zeros(n,m);
r = zeros(n+m,1);
dx = zeros(n,1); % displacement in x
dy = zeros(m,1); % displacement in y
s = zeros(n,1);
Z = zeros(m,m);
N = 10; % number of iterations
for i=1:20
	for k=1:n
		s(k) = -(g(k) + J(k,:)*x - theta/x(k)); % slack vector
		r(k) = theta - x(k)*s(k); % theta - x*s
		B(k,:) = -x(k)*(A(:,k))'; % stems from s dx + x ds = theta - x*s
	end
	K = [diag(s) B;A Z];
	d = K\r;
	dx = d(1:n);
	dy = d(n+1:n+m);
	x = x+dx;
	y = y+dy;
	theta = theta / 2;
end
