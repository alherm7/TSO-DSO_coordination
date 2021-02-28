function [DA_market_outcome] = DA_market_clearing(data_tso,dual_DA_gen,dual_DA_dem,...
	dual_day_ahead_wind, cost_RT,iter,p_gen_DA_hat,p_dem_DA_hat,wind_DA_hat,p_gen_DA_tilde,p_dem_DA_tilde,kk)


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

cvx_solver Gurobi


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

nw = length(windgL);

%% define variables
cvx_begin quiet

% cvx_precision high

variable p_gen_DA(ng)
variable p_dem_DA(nd) nonnegative
variable wind_DA(nw) nonnegative
variable shed_p nonnegative

variable alpha_cut(nscen)

dual variable lambda

expression cost_gen(ng)
expression cost_dem(nd)
%% define cost function
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
cost = cost_DA;
%% objective statement
minimize(cost)

subject to
%% constraints
% Day ahead market
lambda : sum(p_gen_DA) + sum(wind_DA) + shed_p == sum(p_dem_DA);

for DSO_num = 1:length(DSO_codes)
	for k = 1:ng
		if any(k == gens_area{DSO_num})
			Pmin_gen(k) <= p_gen_DA(k) <= p_gen_DA_tilde(k);
		elseif any(k == gens_T)
			Pmin_gen(k) <= p_gen_DA(k) <= Pmax_gen(k);
		end
	end
	
	for k = 1:nd
		for m = 1:length(dems_area{DSO_num})
			if any(k == dems_area{DSO_num}(m))
				Pmin_dem(k) <= p_dem_DA(k) <= p_dem_DA_tilde(k);
			elseif any(k == dems_T)
				Pmin_dem(k) <= p_dem_DA(k) <= Pmax_dem(k);
			end
		end
	end
	

end


0 <= wind_DA <= Wmax_mean_DA;

% p_gen_DA == fcoo_pg_DA;
% p_dem_DA == fcoo_pd_DA;

%%%%%%%%%%%%%%%%%
%% Benders cuts
alpha_min <= alpha_cut;

for k = 1:iter
	for s = 1:nscen
		alpha_cut(s) >= cost_RT{k}(s) + sum(dual_DA_gen{k}(:,s) .* ( p_gen_DA - p_gen_DA_hat{k}  ) )...
			+ sum(dual_DA_dem{k}(:,s) .* ( p_dem_DA - p_dem_DA_hat{k} ) ) ...
			+ sum(dual_day_ahead_wind{k}(:,s) .* (wind_DA - wind_DA_hat{k} ) );
		
	end
end

cvx_end

%% output

DA_market_outcome.p_gen_DA = p_gen_DA;
DA_market_outcome.p_dem_DA = p_dem_DA;
DA_market_outcome.cost_DA = cost_DA;
DA_market_outcome.lambda = lambda;
DA_market_outcome.cost_cuts = cost_cuts;
DA_market_outcome.cost = cost;
DA_market_outcome.alpha_cut = alpha_cut;
DA_market_outcome.wind_DA = wind_DA;
DA_market_outcome.shed_p = shed_p;
DA_market_outcome.cost_dem = cost_dem;
DA_market_outcome.cost_gen = cost_gen;


if DA_market_outcome.shed_p > 10e-5
	warning(['WPP: ' num2str(kk) ', Conventional Market, Benders Master Iteration: ' num2str(iter) ', There is load shedding in the Day-Ahead market']);
	pause(3)
end

end

