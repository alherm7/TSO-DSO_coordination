function [results_PCCoptim, results_full_coordination, results_conv_market,results_DA_PF_market,results_DA_scen_market, varargout] ...
	= run_insample_optim(mpc,base_PCC_capacity,PCC_lim_factor_grided,...
	wind_factor_grided,line_factor_gridded,line_var_cap,variance_base,mean_base,PCC_branch_id,nscen,kk)
nout=nargout;
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
variance_wind = variance_base;
mean_wind = mean_base*wind_factor_grided(kk);

% 	variance = [0.01]*wind_factor(kk);
% 	mean = [0.8]*wind_factor(kk);

[scenarios, normalized_sigma] = scenario_generator(mpc,variance_wind,mean_wind,300,'default',false);
% figure(1)
% plot(scenarios(1:30,1:3))
% figure(2)
% plot(scenarios(1:30,4:7))
if nargout>5
    varargout{1} = normalized_sigma;
end

for k = 1:nscen
	mpc.RTscen(:,3+3*(k-1)) = scenarios(k,:);
	mpc.RTscen(:,5+3*(k-1)) = 1/nscen;
	mpc.RTscen(:,4+3*(k-1)) = 0;
	mpc.RTscen(:,2) = nscen;
end

mpc.branch(PCC_branch_id,6) = base_PCC_capacity*PCC_lim_factor_grided(kk);

mpc.branch(line_var_cap,6) = line_factor_gridded(kk)*line_cap_base + 0.01;

tsim = tic;
results_PCCoptim = PCC_optim(mpc,kk);
% results_PCCoptim.DA = [1];
time_PCC = toc(tsim);
tsim = tic;
results_full_coordination = TSO_DSO_cooptim(mpc);
% results_full_coordination.dummy = [1];
time_FC = toc(tsim);
tsim = tic;
results_conv_market = conventional_market_clearing(mpc,kk);
% results_conv_market.dummy = [1];
time_CM = toc(tsim);
results_DA_PF_market = market_clearing_DA_PF(mpc,kk);
% results_DA_PF_market.dummy = 1;
results_DA_scen_market = market_clearing_DA_SCEN(mpc,kk);
% results_DA_scen_market.dummy = [1];
disp('-----------------------------')
disp('-----------------------------')
disp('-----------------------------')
disp('-----------------------------')
disp('')
disp(['WPP-IS  FINISHED: ' num2str(kk) ', Times: PCC-optim ' num2str(time_PCC) ' , FC ' num2str(time_FC) ', CM ' num2str(time_CM)])
disp('')
disp('-----------------------------')
disp('-----------------------------')
disp('-----------------------------')
disp('-----------------------------')
disp('-----------------------------')

end
