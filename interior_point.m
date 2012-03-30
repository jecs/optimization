% A script that approximates the minimum of a quadratic (zero-offset) function
% F(x) = q'x + 1/2 x'Kx given a set of constraints Dx-e=0 (Rp) and Cx-b >= 0 (Rm)

%%%% PARAMETER DECLARATIONS
% Strang LP Example (~30 iterations for convergence)
% function to minimize
q = [5;3;8];
K = zeros(3,3);

% constraints
C = [1 0 0;0 1 0;0 0 1]; b = [0 0 0]';
D = [1 1 2]; e = [4];

w = 4/3; % initial barrier weight
wm = 4/3*1/64; % (w minimum) w at which the program will terminate
wup = 1e-2; % (w update) ratio of dx/x, dy/y or dL/L at which w will be shrunk
ov = 0.2; % tolerated fractional overshoot
x = [1 1 1]'; y = [1 1 1]'; L = [0.5]; % initial point

% Mattingley QP Example
% function to minimize
% constraints

%%%% END PARAMETER DECLARATIONS

%%% DIMENSIONS
n = size(K, 1);
m = size(C, 1);
p = size(D, 1);

%%% DISPLACEMENT VECTORS
dx = x; dy = y; dL = L;

%%% AUXILIARY VECTORS & MATRICES
s = zeros(m,1); % Cx-b 
S = zeros(m,m); % diagonalized s
O = ones(m,1); % ones vector for the w
Y = zeros(m,m); % diagonalized y

Zpm = zeros(p,m); % static zero matrices
Zpp = zeros(p,p);
Zmp = zeros(m,p);

%%% LOOP
I = 0; % number of iterations
while w > wm
	s = C*x-b;
	S = diag(s);
	W = w*O;
	Y = diag(y);

	% assemble system
	B = [D Zpm Zpp;Y*C S Zmp;K -C' -D'];
	v = [e-D*x;W-Y*s;C'*y+D'*L-q-K*x];

	% solve
	d = B\v;
	dx = d(1:n);
	dy = d(n+1:n+m);
	dL = d(n+m+1:n+m+p);

	% check for overshoot / tentative
	px = max(dx./x);
	py = max(dy./y);
	pL = max(dL./L);
	pm = max([px, py, pL]);
	if pm > ov
		a = ov/pm;
		dx = a*dx;
		dy = a*dy;
		dL = a*dL;
	end

	% update
	x = x+dx;
	y = y+dy;
	L = L+dL;
	I = I+1;

	% occasionally shrink factor
	if pm < wup 
		w = w/2;
	end
end

x
y
L
I
