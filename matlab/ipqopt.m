% [x, L, I, fo] = ipqopt(f, g, h, x, dgap, mu, eps, epsf, alpha, beta)
%
% Applies an interior-point primal-dual method to a function f
% subject to the inequalities g_1(x) <= 0, ..., g_m(x) <= 0
% and to the equality h(x) = Ax-b = 0.
%
% Returns x & L, the primal and dual optima, respectively; I,
% the number of major-loop iterations, and fo = f(x*).
%
% f, g, and h are structs that adhere to the following hierarchy
%
%		  f
%   c   l   q   n
%
%		  g
%   c   l   q   m
%
%         h
%     A   b   p
%
% The elements c, l, and q reflect the constant, linear, and quadratic
% terms of their respective function(s). l's must be column vectors.
% The q's of f and g can be left empty; ie, f.q = []; or g.q = [];
% Empty q's accelerate the execution of the function. h can also have
% an empty A and b, but its p must equal 0.
%
% The elements of each g_k are concatenated horizontally in the corres-
% ponding elements of g; ie,
% g.c = [c1 c2 c2 ... ];
% g.l = [l0 l1 l2 ... ];
% g.q = [q0 q1 q2 ... ];
%
% n, m, p reflect the length of x, the number of inequality constaints,
% and the number of equality constraints, respectively.
%
% Other arguments:
% x     - the starting point
% dgap  - the initial duality gap of the function (must not be too low)
% mu    - the duality gap reduction factor
% eps   - the duality gap tolerance
% epsf  - the tolerance on the norms of the residual primal and dual 
%         vectors
% alpha - the line search step threshold (usually > 0.01 and < 0.1)
% beta  - the factor by which the step is reduced during line search
%         (usually > 0.3 and < 0.8)
%
% Example of the instantiation of f, g, and h: 
% f = struct('c', [0], 'l', [5 3 8]', 'q', [], 'n', 3);
% g = struct('c', [0 0 0], 'l', -[[1 0 0]', [0 1 0]', [0 0 1]'], \
%			 'q', [], 'm', 5);
% h = struct('A', [1 1 2], 'b', [4], 'p', 1);

function [x, L, I, fo] = ipqopt(f, g, h, x, dgap, mu, eps, epsf, alpha, beta)
	%%%%%%%%%% INSTANTIATION
	isfq = ~isempty(f.q); % is f quadratic?
	isgq = ~isempty(g.q); % is g quadratic?
	iseq = ~isempty(h.A); % are there equality constraints?

	n = f.n;
	m = g.m;
	p = h.p;

	rdua = zeros(n, 1);
	rcen = zeros(m, 1);
	rpri = zeros(p, 1);
	rt = zeros(n+m+p, 1);
	B = sparse(zeros(n+m+p, n+m+p));

	t = 0;    % inverse barrier weight
	invt = 0;

	L = zeros(m, 1);
	v = zeros(p, 1);
	y = zeros(n+m+p, 1);

	dx = zeros(n, 1);
	dL = zeros(m, 1);
	dv = zeros(p, 1);
	dy = zeros(n+m+p, 1);

	ge = zeros(m, 1);
	Dg = zeros(m, m);
	lgq = zeros(n, n);
	gqx = zeros(m, 1);
	onesv = ones(m,1);
	
	smax = zeros(m, 1);

	I = 0; % number of iterations


	%%%%%%%%%% LOOP
	%%%%% calculate g(x), Dg, and, and sum of weigthed g Hessians
	[ge, gqx] = eval_g(g, x, isgq, ge, gqx);
	Dg = eval_Dg(g, x, gqx, isgq, Dg);
	L = -dgap./ge; % initialize L
	lgq = eval_lgq(g, L, isgq, lgq);
	while (1)
	%%%%% I. Determine t
		dgap = -ge'*L;
		t = mu*m/dgap;
		invt = 1/t;

	%%%%% II. Compute dy
		% calculate rt and build rt
		rdua = f.l+Dg'*L;
		if (iseq) 
			rdua = rdua+h.A'*v;
			rpri = h.A*x-h.b;
		end
		if (isfq) 
			rdua = rdua+f.q*x; 
		end
		rcen = -diag(L)*ge-invt*ones(m,1);
		rt = [rdua;rcen;rpri];

		% build block matrix and rt
		if (isfq) B(1:n,1:n) = (f.q+lgq); end
		if (isgq) B(1:n,1:n) = B(1:n,1:n)+lgq; end
		B(1:n,n+1:n+m) = Dg';
		
		B(n+1:n+m,1:n) = -diag(L)*Dg;
		B(n+1:n+m,n+1:n+m) = -diag(ge);
		
		if (iseq) 
			B(1:n,n+m+1:end) = h.A';
			B(n+m+1:end,1:n) = h.A;
		end

		% solve system
		dy = -(B\rt);
		dx = dy(1:n);
		dL = dy(n+1:n+m);
		dv = dy(n+m+1:n+m+p);

	%%%%% III. Line search
		% calculate smax
		for i=1:m
			if dL(i) < 0 && -dL(i) > L(i) % avoid dividing
				smax(i) = -L(i)/dL(i);
			else
				smax(i) = 1;
			end
		end
		s = 0.99*min(smax);

		% decrease s until g(xp) < 0
		while (1)
			xp = x+s*dx;
			[ge, gqx] = eval_g(g, xp, isgq, ge, gqx);
			if (all(ge < zeros(m,1))) break;
			else s = beta*s;
			end
		end

		% decrease s until suitable step size is found
		oldnorm = norm(rt);
		while (1)
			xp = x+s*dx;
			Lp = L+s*dL;
			vp = v+s*dv;

			[ge, gqx] = eval_g(g, xp, isgq, ge, gqx);
			Dg = eval_Dg(g, xp, gqx, isgq, Dg);

			rdua = f.l+Dg'*L;
			if (iseq) 
				rdua = rdua+h.A'*vp;
				rpri = h.A*xp-h.b;
			end
			if (isfq) 
				rdua = rdua+f.q*xp; 
			end
			rcen = -diag(Lp)*ge-invt*ones(m,1);
			rt = [rdua;rcen;rpri];
			
			if (norm(rt) <= (1-alpha*s)*oldnorm)
				break;
			else
				s = beta*s;
			end
		end
		x = xp;
		L = Lp;
		v = vp;

		I = I+1;

		%%%%% IV. Check if x is optimal enough
		if (norm(rpri) <= epsf && norm(rdua) <= epsf && dgap <= eps)
			break;
		end
	end

	fo = f.c+f.l'*x+0.5*x'*f.q*x;
end

% vectorized evaluations
function [ge, gqx] = eval_g(g, x, isgq, ge, gqx)
	bigx = repmat(x,g.m,1);
	if (isgq)
		gqx = g.q*bigx;
		ge = g.c' + (g.l' + 0.5*gqx')*x;
	else
		gqx = gqx;
		ge = g.c' + g.l'*x;
	end
end

function [Dg] = eval_Dg(g, x, gqx, isgq, Dg)
	if (isgq)
		Dg = g.l'+gqx;
	else
		Dg = g.l';
	end
end

function [lgq] = eval_lgq(g, L, isgq, lgq)
	if (isgq)
		Lr = kron(L,speye(g.m,g.m));
		lgq = g.q*Lr;
	else
		lgq = lgq;
	end
end	
