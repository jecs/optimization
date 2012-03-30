% A script to produce approximate the minimum of a quadratic (zero-offset) function F(x) = q'x + 1/2 x'Kx given a set of constraints Ax = b (where n > m) and x >= 0

% the specific values are from an example from Strang's textbook

% parts of function to minimize
q = [5;3;8];
K = zeros(3,3);

% constraints
A = [1 1 2];
b = [4];

% relevant dimensions
n = size(K, 1);
m = size(A, 1);

% barrier weight
theta = 1e-3;
theta_lower_bound = 0.01;

% starting point
x = ones(3,1); % iterative x variable
y = 2; % iterative lagrange multiplier

% auxiliary vectors and matrices
dx = x; % displacement in x
dy = y; % displacement in y
s = zeros(n,1);

Z = zeros(m,m);
r = zeros(n+m,1);
d = zeros(n+m,1); % displacement vector
O = ones(n,1); % ones vector
Th = zeros(n,1); % theta scalar vector
B = zeros(n+m,n+m); % big KKT matric

% diagonalized matrices
S = zeros(n, n);
X = zeros(n, n);

N = 100; % number of iterations
I = 0; % number of iterations
while max(dx)/max(x) > 1e-3
	Th = theta*O;
	s = q-Th./x;
	S = diag(s);
   	X = diag(x);
	r = [Th-(S*X*O); zeros(m, 1)];
	B = [(S+X*K) (-X*A');A Z];
	d = B\r;
	dx = d(1:n);
	dy = d(n+1:n+m);
	x = x+dx;
	y = y+dy;
	theta = theta;
	I = I+1;
end

x
y
I
