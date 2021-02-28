function [DA_market_outcome] = PCC_optimizer_DA_KKT_simple_EVPI(data_tso,dual_DA_gen,dual_DA_dem,...
	dual_day_ahead_wind,cost_RT,iter,p_gen_DA_hat,p_dem_DA_hat,wind_DA_hat,kk,s)


%% load parameters
VOLL_DA = 10000;
alpha_min = -100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('results_full_coord.mat');
% fcoo_pg_DA = results_full_coordination.p_gen_DA;
% fcoo_pd_DA = results_full_coordination.p_dem_DA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc{1} = double.empty(1,0);
[nb,busL,busN,busType,Pd,Qd,Vmax,Vmin,statusG,activeG,activeD,...
	ng,nd,genL,genN,incidentG,demL,incidentD,incidentW,Qmax_gen,Qmin_gen,...
    Pmax_gen,Pmin_gen,Qmax_dem,Qmin_dem,Pmax_dem,Pmin_dem,statusL,...
	activeL,nl,fbusL,tbusL,SlmMax,fbusN,tbusN,incidentF,incidentT,...
    Yf,Yt,YfP,YtP,Ybus,edges,offer_gen_DA,offer_gen_upreg,...
	offer_gen_downreg,bid_dem_DA,bid_dem_upreg,bid_dem_downreg,...
	nr_vals_cost,p,cost_type,ys,areas,bus_areacodes,....
	incident_singlephase,Y_ft,rBranch,xBranch,narea,prob_wscen,Wmax,Wmin,...
	n_wgen,nscen,Wmax_mean_DA,windgL,...
	offer_wind_DA,offer_wind_up,offer_wind_dn] = Data_Reader(data_tso,cc);


count = 0;
err_count = 0;
while count == err_count
    try
        cvx_solver Mosek
    catch
        pause(20)
        err_count = err_count + 1;
    end
    count = count + 1;
end
disp(['Error count: ', num2str(err_count)])

nw = length(windgL);

[area_codes, n_areas, DSO_codes, fbusL_area, tbusL_area, branch_num_area, ...
	fbusL_ext_area, tbusL_ext_area, branch_num_ext_area, overlap,...
	overlap_t, neigh_area, ext_area, ext_area_singlephase,fbus_local,tbus_local,incidence_area,PCC_branch_id] ...
	= find_overlap(bus_areacodes, fbusL, tbusL, areas);

ne = n_areas-1;
gens_DSO = [];
dems_DSO = [];
windg_DSO = [];

for DSO_num = 1:length(DSO_codes)
	gens_area{DSO_num} = find(ismember(genL,areas{DSO_num+1}));
	dems_area{DSO_num} =  find(ismember(demL,areas{DSO_num+1}));
	windg_area{DSO_num} = find(ismember(windgL,areas{DSO_num+1}));
	
	gens_DSO = [gens_DSO; gens_area{DSO_num}];
	dems_DSO = [dems_DSO; dems_area{DSO_num}];
	windg_DSO = [windg_DSO; windg_area{DSO_num}];
end



windg_area_T = find(ismember(windgL,areas{1}));
gens_T = find(ismember(genL,areas{1}));
dems_T = find(ismember(demL,areas{1}));


ng_DSO = length(gens_DSO);
nd_DSO = length(dems_DSO);

% p_gen_DA_tilde = [1; 1; 1.5];
% p_dem_DA_tilde = [2; 2; 3; 3];

%% define variables
cvx_begin quiet

% cvx_precision high
variable p_gen_DA_tilde(ng)
variable p_dem_DA_tilde(nd)

variable varsigma_g_plus(ng)
variable varsigma_g_minus(ng)

variable sigma_g_plus(ng)
variable sigma_g_minus(ng)

variable varsigma_d_plus(nd)
variable varsigma_d_minus(nd)

variable sigma_d_plus(nd)
variable sigma_d_minus(nd)

variable lambda_DA

variable nu_minus(nw)
variable nu_plus(nw)

variable rho_minus
variable rho_plus

variable p_gen_DA(ng)
variable p_dem_DA(nd)
variable wind_DA(nw)
variable shed_p

variable alpha_cut(nscen)

variable wind_comp(nw,2) binary
variable gen_comp(ng,4) binary
variable dem_comp(nd,4) binary
variable shed_comp(2) binary
%% define cost function
expression cost_gen(ng)
expression cost_dem(nd)
for m = 1:ng
	if cost_type(m) == 2
		for k = 1:nr_vals_cost
			cost_gen(m) = cost_gen(m) + sum(offer_gen_DA(m,k).*p_gen_DA(m).^(nr_vals_cost-k));
		end
	elseif cost_type(m) == 1
		error('Wrong cost function, piecewise linear is not convex');
	end
end

for m = 1:nd
	cost_dem(m) = bid_dem_DA(m)*p_dem_DA(m);
end

% cost_shed_T = s_shed_T_DA*VOLL_T;
% cost_shed_E = sum(s_shed_E_DA*VOLL_E);
cost_shed = VOLL_DA*shed_p;

% cost_wind_e = sum(p_wind_E_DA*offer_wind);
% cost_wind_T = p_wind_T_DA*offer_wind;
cost_wind = wind_DA.*offer_wind_DA;

cost_DA = sum(cost_gen)-sum(cost_dem) + cost_shed + sum(cost_wind);


cost_cuts = prob_wscen(1,:) * alpha_cut;
cost = cost_DA + cost_cuts;
%% objective statement
minimize(cost)

subject to
%% stationary constraints

% for DSO_num = 1:length(DSO_codes)
for g = 1:ng
	if any(g == gens_DSO)
		offer_gen_DA(g,2) - (varsigma_g_minus(g) - varsigma_g_plus(g)) - lambda_DA == 0;
	elseif any(g == gens_T)
		offer_gen_DA(g,2) - (sigma_g_minus(g) - sigma_g_plus(g)) - lambda_DA == 0;
	end
end
% end
% for DSO_num = 1:length(DSO_codes)
for d = 1:nd
	if any(d == dems_DSO)
		-bid_dem_DA(d) - (varsigma_d_minus(d) - varsigma_d_plus(d)) + lambda_DA - rho_plus == 0;
	elseif any(d == dems_T)
		 -bid_dem_DA(d) - (sigma_d_minus(d) - sigma_d_plus(d)) + lambda_DA - rho_plus == 0;
	end
end
% end

VOLL_DA - lambda_DA - rho_minus + rho_plus == 0;

for w = 1:nw
% 	if any(w == windg_DSO)
		offer_wind_DA - lambda_DA - nu_minus(w) + nu_plus(w) == 0;
% 	end
end
%% complimentarity constraints
for g = 1:ng
	if any(g == gens_DSO)
% 		gc = find(gens_DSO == g);
		
		varsigma_g_minus(g) <= VOLL_DA * gen_comp(g,1);
		(p_gen_DA(g) - Pmin_gen(g)) <= VOLL_DA *  (1 - gen_comp(g,1));

		varsigma_g_plus(g)  <= VOLL_DA * gen_comp(g,2);
		p_gen_DA_tilde(g) - p_gen_DA(g) <= VOLL_DA * (1 - gen_comp(g,2));
	elseif any(g == gens_T)
		sigma_g_minus(g) <= VOLL_DA * gen_comp(g,3);
		(p_gen_DA(g) - Pmin_gen(g)) <= VOLL_DA *  (1 - gen_comp(g,3));

		sigma_g_plus(g)  <= VOLL_DA * gen_comp(g,4);
		Pmax_gen(g) - p_gen_DA(g) <= VOLL_DA * (1 - gen_comp(g,4));
	end
end

for d = 1:nd
	if any(d == dems_DSO)
% 		dc = find(dems_DSO == d);
		
		varsigma_d_minus(d) <= VOLL_DA * dem_comp(d,1);
		(p_dem_DA(d) - Pmin_dem(d)) <= VOLL_DA *  (1 - dem_comp(d,1));

		varsigma_d_plus(d)  <= VOLL_DA * dem_comp(d,2);
		p_dem_DA_tilde(d) - p_dem_DA(d) <= VOLL_DA * (1 - dem_comp(d,2));
	elseif any(d == dems_T)
		sigma_d_minus(d) <= VOLL_DA * dem_comp(d,3);
		(p_dem_DA(d) - Pmin_dem(d)) <= VOLL_DA *  (1 - dem_comp(d,3));

		sigma_d_plus(d)  <= VOLL_DA * dem_comp(d,4);
		Pmax_dem(d) - p_dem_DA(d) <= VOLL_DA * (1 - dem_comp(d,4));
	end
end

nu_minus <= VOLL_DA *wind_comp(:,1);
wind_DA <= VOLL_DA * (1 - wind_comp(:,1));

nu_plus <= VOLL_DA * wind_comp(:,2);
Wmax_mean_DA - wind_DA <= VOLL_DA * (1 - wind_comp(:,2));

rho_minus <= 1.5*VOLL_DA * shed_comp(1);
shed_p <= VOLL_DA * (1 - shed_comp(1));

rho_plus <= VOLL_DA * shed_comp(2);
sum(p_dem_DA) - shed_p <= VOLL_DA * (1 - shed_comp(2));

%% primal constrains

sum(p_gen_DA) + sum(wind_DA) + shed_p == sum(p_dem_DA);


% for DSO_num = 1:length(DSO_codes)
	for k = 1:ng
		if any(k == gens_DSO)
% 			gc = find(gens_DSO == k);
			
			Pmin_gen(k) <= p_gen_DA(k) <= p_gen_DA_tilde(k);
		else
			Pmin_gen(k) <= p_gen_DA(k) <= Pmax_gen(k);
			p_gen_DA_tilde(k) == 0;
		end
	end
	
	for k = 1:nd
% 		for m = 1:length(dems_area{DSO_num})
		if any(k == dems_DSO)
% 			dc = find(dems_DSO == k);

			Pmin_dem(k) <= p_dem_DA(k) <= p_dem_DA_tilde(k);
		else
			Pmin_dem(k) <= p_dem_DA(k) <= Pmax_dem(k);
			p_dem_DA_tilde(k) == 0;
		end
% 		end
% 	end
	

% 	0 <= p_wind_E_DA(DSO_num) <= sum(Wmax_mean_DA(windg_area{DSO_num}));
end

0 <= wind_DA <= Wmax_mean_DA;
0 <= shed_p <= sum(p_dem_DA);
% 0 <= p_wind_T_DA <= sum(Wmax_mean_DA(windg_area_T));

% p_gen_DA == fcoo_pg_DA;
% p_dem_DA == fcoo_pd_DA;

p_gen_DA_tilde <= Pmax_gen;
p_dem_DA_tilde <= Pmax_dem;

%%%%%%%%%%%%%%%%%

%% dual feasibility
varsigma_g_minus >= 0;
varsigma_g_plus >= 0;
varsigma_d_minus >= 0;
varsigma_d_plus >= 0;

sigma_g_minus >= 0;
sigma_g_plus >= 0;
sigma_d_minus >= 0;
sigma_d_plus >= 0;

nu_minus >= 0;
nu_plus >= 0;

rho_minus >= 0;
rho_plus >= 0;

%% Benders cuts
alpha_min <= alpha_cut;

for k = 1:iter
	for s = 1:nscen
		alpha_cut(s) >= cost_RT{k}(s) + sum(dual_DA_gen{k}(:,s) .* ( p_gen_DA - p_gen_DA_hat{k}  ) )...
			+ sum(dual_DA_dem{k}(:,s) .* ( p_dem_DA - p_dem_DA_hat{k} ) ) ...
			+ sum(dual_day_ahead_wind{k}(:,s) .* (wind_DA - wind_DA_hat{k} ) );
	end
end
time_master_solve = tic;
cvx_end
fin_time_master_solve = toc(time_master_solve);
%% output

DA_market_outcome.p_gen_DA = p_gen_DA;
DA_market_outcome.p_dem_DA = p_dem_DA;
DA_market_outcome.cost_DA = cost_DA;
DA_market_outcome.lambda = lambda_DA;
DA_market_outcome.cost_cuts = cost_cuts;
DA_market_outcome.cost = cost;
DA_market_outcome.alpha_cut = alpha_cut;
DA_market_outcome.wind_DA = wind_DA;
% DA_market_outcome.p_wind_T_DA = p_wind_T_DA;
% DA_market_outcome.p_wind = p_wind_E_DA + p_wind_T_DA;
DA_market_outcome.shed_p = shed_p;
% DA_market_outcome.p_gen_DA_tilde = zeros(ng,1);
% DA_market_outcome.p_dem_DA_tilde = zeros(nd,1);
DA_market_outcome.p_gen_DA_tilde = p_gen_DA_tilde;
DA_market_outcome.p_dem_DA_tilde = p_dem_DA_tilde;
DA_market_outcome.fin_time_master_solve = fin_time_master_solve;


if DA_market_outcome.shed_p > 10e-5
	warning(['WPP: ' num2str(kk) ', PCC optimizer, Benders Master Iteration: ' num2str(iter) ', There is load shedding in the Day-Ahead market']);
	pause(3)
end


end