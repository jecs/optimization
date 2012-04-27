% script dedicated to testing ipqopt

test = 2;

switch test
case 0
	% from Strang's textbook
	f = struct('c', 0, 'l', [5 3 8]', 'q', [], 'n', 3);
	g = struct('c', [0 0 0], 'l', -[[1 0 0]' [0 1 0]' [0 0 1]'], 'q', [], 'm', 3);
	h = struct('A', [1 1 2], 'b', [4], 'p', 1);

	x_init = [0.8 1 1]';
	[x, L, I, fo] = ipqopt(f, g, h, x_init, 10, 10, 1e-3, 1e-3, 0.1, 0.8);
case 1
	% from http://www.akiti.ca/QuadProgEx0Constr.html
	f = struct('c', 0, 'l', [-2 -6]', 'q', [1 -1;-1 2], 'n', 2);
	g = struct('c', [-2 -2 -3 0 0], 'l', [[1 1]', [-1 2]', [2 1]', [-1 0]', [0 -1]'], 'q', [], 'm', 5);
	h = struct('A', [], 'b', [], 'p', 0);
	x_init = [0.8 0.8]';
	[x, L, I, fo] = ipqopt(f, g, h, x_init, 1, 10, 1e-5, 1e-5, 0.01, 0.8);
case 2
	% from http://trac.openopt.org/openopt/browser/PythonPackages/OpenOpt/openopt/examples/qp_1.py
	f = struct('c', 0, 'l', [15 8 80]', 'q', [1 0 0;0 2 0;0 0 3], 'n', 3);
	g = struct('c', [-150 -800 -15], 'l', [[1 2 3]' [8 15 80]' [1 0 0]'], 'q', [], 'm', 3);
	h = struct('A', [0 1 -1], 'b', 25.5, 'p', 1);
	x_init = [0 25.5 0]';
	[x, L, I, fo] = ipqopt(f, g, h, x_init, 1, 10, 1e-5, 1e-5, 0.01, 0.8);
end

x
L
I
fo
