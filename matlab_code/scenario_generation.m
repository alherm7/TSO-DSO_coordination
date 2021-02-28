clc
clear all
close all


%%

test_case = 'case_24TSO_3DSO';
mpc = eval(test_case);

nwgen = length(mpc.wind_loc(:,1));
x_loc = mpc.wind_loc(:,2);
y_loc = mpc.wind_loc(:,3);

for w1 = 1:nwgen
	for w2 = 1:nwgen
		dist_ft_wgen(w1,w2) = norm([x_loc(w1) - x_loc(w2), y_loc(w1) - y_loc(w2)]);
	end
end

variance = [24	27	28	32	45	23 33];

for w1 = 1:nwgen
	for w2 = 1:nwgen
		sigma(w1,w2) = (variance(w1) + variance(w2))/2*1/(exp(dist_ft_wgen(w1,w2)));
	end
end

mean = [170	130	140	100	90	90	22];

gm = gmdistribution(mean,sigma);

rng('default');
for k = 1:100
	y(k,:) = random(gm);
end

plot(y)

