function [results_PCCoptim,results_PCCoptim_avg] ...
	= run_insample_optim_EVPI(mpc,base_PCC_capacity,PCC_lim_factor_grided,...
	wind_factor_grided,line_factor_gridded,line_var_cap,variance_base,mean_base,PCC_branch_id,nscen,kk)

mean_wind = mean_base;
variance_wind = variance_base;

line_cap_base = mpc.branch(line_var_cap,6);

disp([' '])
disp([' '])
disp([' '])
disp([' IN Sample optimization '])
disp(['Wind-penetration point ' num2str(kk) ])
disp([' '])

% variance_wind(1:3) = variance_base(1:3)*wind_factor_grided(kk);
% mean_wind(1:3) = mean_base(1:3)*wind_factor_grided(kk);
variance_wind = variance_base*wind_factor_grided(kk);
mean_wind = mean_base*wind_factor_grided(kk);

% 	variance = [0.01]*wind_factor(kk);
% 	mean = [0.8]*wind_factor(kk);

scenarios = scenario_generator(mpc,variance_wind,mean_wind,300,'default',false);

for k = 1:nscen
	mpc.RTscen(:,3+3*(k-1)) = scenarios(k,:);
	mpc.RTscen(:,5+3*(k-1)) = 1;
	mpc.RTscen(:,4+3*(k-1)) = 0;
	mpc.RTscen(:,2) = nscen;
end

mpc.branch(PCC_branch_id,6) = base_PCC_capacity*PCC_lim_factor_grided(kk);

mpc.branch(line_var_cap,6) = line_factor_gridded(kk)*line_cap_base + 0.01;

tsim = tic;
parfor s = 1:nscen
    mpc2 = mpc;
    mpc2.RTscen(:,3) = scenarios(s,:);
	mpc2.RTscen(:,5) = 1;
	mpc2.RTscen(:,4) = 0;
	mpc2.RTscen(:,2) = 1;
    results_PCCoptim(s) = PCC_optim(mpc2,kk);
end

results_PCCoptim_avg = []
for s = 1:nscen
    results_PCCoptim_avg = [results_PCCoptim_avg results_PCCoptim(s).DA_outcome.cost];
end
results_PCCoptim_avg = mean(results_PCCoptim_avg);

time_PCC = toc(tsim);
tsim = tic;
% results_full_coordination = TSO_DSO_cooptim(mpc);
time_FC = toc(tsim);
tsim = tic;
% results_conv_market = conventional_market_clearing(mpc,kk);
time_CM = toc(tsim);
results_full_coordination = 1;
results_conv_market = 1;
disp(['WPP-IS: ' num2str(kk) ', Times: PCC-optim ' num2str(time_PCC) ' , FC ' num2str(time_FC) ', CM ' num2str(time_CM)])
	
end
