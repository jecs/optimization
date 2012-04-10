% A script that approximates the minimum of a quadratic function
% using a primal-dual interior point method
%
% f(x) = cf + lf'x + 1/2 x'Qfx (Rn -> R)
% given a set of constraints
% gk(x) = cgk + lgk'x + 1/2 x'Qgkx <= 0 (Rn -> R) for k=1...m
% g(x) <= 0 (Rn -> Rm)
% Ax = b (Rn -> Rp)
%
% (Convention: c_ - constant, l_ gradient, Q_ Hessian)

%%%% PARAMETER DECLARATIONS
% Strang LP Example (~30 iterations for convergence)
% function to minimize
cf = 0; lf = [5;3;8]; Qf = zeros(3,3);

% inequality constraints
cg = [0 0 0]; lg = -[[1 0 0]' [0 1 0]' [0 0 1]']; Qg = zeros(3,3,3);

% equality constraints
A = [1 1 2]; b = [4];

% Mattingley QP Example
% function to minimize
% constraints

%%% DIMENSIONS
n = size(lf, 1);
m = size(lg, 2);
p = size(A, 1);

%%% CHOOSE INITIAL POINT
% x = zeros(n, 1); % (all zeros for now)
x = [0.8 1 1]';
L = zeros(m,1);
for i=1:3
	L(i) = -100/(cg(i)+lg(:,i)'*x);
end
v = 1;

%%% PARAMETERS
mu = 10;      % barrier weight decay factor
alpha = 1e-2; % line search confidence
beta = 0.5;   % line search constant
epsf = 1e-5;  % primal-dual tolerance
eps = 1e-3;   % duality gap tolerance

%%% AUXILIARY VECTORS & MATRICES
rdua = zeros(n, 1);
rcen = zeros(m, 1);
rpri = zeros(p, 1);
rt = zeros(n+m+p, 1);

dgap = 0; % surrogate duality gap
t = 0;    % inverse barrier weight

dx = zeros(n, 1);
dL = zeros(m, 1);
dv = zeros(p, 1);
dy = zeros(n+m+p, 1);

g = zeros(m, 1);
LGS = zeros(m,m);

smax = zeros(m,1);

% calculate g(x), Dg, and sum of weighted g Hessians
for i=1:m
	xQg = x'*Qg(:,:,i);
	Dg(i,:) = lg(:,i)'+xQg;
	g(i) = cg(i)+lg(:,i)'*x+0.5*xQg*x;
	LGS = LGS + L(i)*Qg(:,:,i);
end

%%% LOOP
I = 0; % number of iterations
while (1)
	% calculate duality gap
	dgap = -g'*L;

	% determine t
	t = mu*m/dgap;
	invt = 1/t;

	% calculate rt
	rdua = lf+Qf*x+Dg'*L+A'*v;
	rcen = -diag(L)*g-invt*ones(m,1);
	rpri = A*x-b;

	% stop execution if optimum is found
	if (norm(rpri) <= epsf && norm(rdua) <= epsf && dgap <= eps)
		break;
	end

	% build block matrix and rt
	BX1 = Qf+LGS;
	BL1 = Dg';
	BV1 = A';

	BX2 = -diag(L)*Dg;
	BL2 = -diag(g);
	BV2 = zeros(m, p);
	
	BX3 = A;
	BL3 = zeros(p, m);
	BV3 = zeros(p, p);

	B = [BX1 BL1 BV1;BX2 BL2 BV2;BX3 BL3 BV3];
	rt = [rdua;rcen;rpri];

	% solve system
	dy = -(B\rt);
	dx = dy(1:n);
	dL = dy(n+1:n+m);
	dv = dy(n+m+1:n+m+p);

	% 5) line search
	% 5a) calculate smax
	for i=1:m
		if dL(i) < 0
	 		smax(i) = -L(i)/dL(i);
		else
			smax(i) = 1;
		end
	end
	s = 0.99*min([smax;1]);

	% 5b) decrease s until suitable step size is found
	oldnorm = norm(rt);
	while (1)
		xp = x+s*dx;
		Lp = L+s*dL;
		vp = v+s*dv;

		for i=1:m
			xQg = xp'*Qg(:,:,i);
			Dg(i,:) = lg(:,i)'+xQg;
			g(i) = cg(i)+lg(:,i)'*xp+0.5*xQg*xp;
		end
		rduap = lf+Qf*xp+Dg'*Lp+A'*vp;
		rcenp = -diag(Lp)*g-invt*ones(m,1);
		rprip = A*xp-b;
		rtp = [rduap;rcenp;rprip];
		
		if (norm(rtp) <= (1-alpha*s)*oldnorm)
			break;
		else
			s = beta*s;
		end
	end
	x = xp;
	L = Lp;
	v = vp;

	I = I+1;
	% DEBUG input('');
end
