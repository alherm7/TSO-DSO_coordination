function [scenarios,varargout] = scenario_generator(data,var_forcast,mean_forcast,nscen,rand_status,if_plot)

nout = nargout;
nwgen = length(data.wind_loc(:,1));
x_loc = data.wind_loc(:,2);
y_loc = data.wind_loc(:,3);
scenarios = zeros(nscen,nwgen);
corr_length = 1;

for w1 = 1:nwgen
	for w2 = 1:nwgen
		dist_ft_wgen(w1,w2) = norm([x_loc(w1) - x_loc(w2), y_loc(w1) - y_loc(w2)]);
	end
end


for w1 = 1:nwgen
	for w2 = 1:nwgen
		sigma(w1,w2) = sqrt(var_forcast(w1) * var_forcast(w2))*(1.1/(exp(dist_ft_wgen(w1,w2)/corr_length))-0.1);
	end
end

for w1 = 1:nwgen
	for w2 = 1:nwgen
		normalized_sigma(w1,w2) = sigma(w1,w2)/(sqrt(var_forcast(w1))*sqrt(var_forcast(w2)));
	end
end
if nout > 1
    varargout{1} = normalized_sigma;
end

gm = gmdistribution(mean_forcast,sigma);

rng(rand_status);
for k = 1:nscen
	scenarios(k,:) = random(gm);
end
scenarios(scenarios<0) = 0;

if if_plot
	figure
	plot(scenarios)
end


end