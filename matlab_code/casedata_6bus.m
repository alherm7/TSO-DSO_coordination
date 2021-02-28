function mpc = casedata_6bus
%CASE5  Power flow data for modified 5 bus, 5 gen case based on PJM 5-bus system
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from ...
%     F.Li and R.Bo, "Small Test Systems for Power System Economic Studies",
%     Proceedings of the 2010 IEEE Power & Energy Society General Meeting

%   Created by Rui Bo in 2006, modified in 2010, 2014.
%   Distributed with permission.

%   MATPOWER

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	1		2		3	4	5	6	7		8	9	10		11		12		13
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1		3		0	0	0	0	0		1	0	1		1		1.1		0.9;
	2		1		1	0	0	0	0		1	0	1		1		1.1		0.9;
	3		1		0	0	0	0	0		1	0	1		1		1.1		0.9;
	4		1		0	0	0	0	1		1	0	1		1		1.1		0.9;
	5		1		1	0	0	0	1		1	0	1		1		1.1		0.9;
	6		1		1	1	0	0	1		1	0	1		1		1.1		0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
% 	1	1	1	4		-4		1.1	1		1		200		0;
	3	1	1	3		-3		1	1		1		15		0;
	4	1	1	3		-3		1	1		1		2		0;
];

%% demand data
%	bus	Pd	Qd	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	
mpc.dem = [
% 	2	1	1	4		0		1.1	1		1		200		0.8;
% 	3	1	1	3		0		1	1		1		3		0.2;
% 	4	1	1	3		0		1	1		1		5		1.1;
	6	1	1	3		0		1	1		1		7		3;
];

%% branch data
%	fbus	tbus	r		x		b	Sl_max	ratio	angle	status
mpc.branch = [
	1		2		0.5		0.04	0	1		0		0		1;
	2		3		0.2		0.04	0	1		0		0		1;
	1		3		0.2		0.04 	0	1		0		0		1;
	1		4		0.006	0.005	0	5.5		0		0		1;
	4		5		0.006	0.005	0	0		0		0		1;
	5		6		0.006	0.005	0	0		0		0		1;
];

%%-----  OPF Data  -----%%
%% generator cost data in day-ahead
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
% 	2	0	0	3	0	20	0;
	2	0	0	3	0	10	0;
	2	0	0	3	0	40	0;
];

%% generator cost data in Real-time
%	1	upreg_price	downreg_price
%	2	upreg_price	downreg_price
mpc.gencost_RT = [
% 	2	5			4;
	2	5.5			5;
	2	10			3;
];

%% demand bid data in day-ahead
mpc.demcost = [
% 	2	0	0	3	0	19	0;
% 	2	0	0	3	0	22	0;
% 	2	0	0	3	0	21	0;
	2	0	0	3	0	25	0;
];

%% demand bid data in Real-time
%	1	upreg_price	downreg_price
%	2	upreg_price	downreg_price
mpc.demcost_RT = [
% 	2	13			2;
% 	2	12			3;
% 	2	11			4;
	2	2			5;
];

%% wind generation scenarios in real-time
%	bus	n_scen	Wmax1	Wmin1	prob1	Wmax2	Wmin2	prob2	Wmax33	Wmin3	prob3 ....
mpc.RTscen = [
	2	3	0.1	0	0.4		0.3	0	0.4		1	0	0.2;
% 	5	3	0.1	0	0.4		0.4	0	0.4		0.5	0	0.2;
% 	6	3	0.1	0	0.4		0.5	0	0.4		0.6	0	0.2;
	];

