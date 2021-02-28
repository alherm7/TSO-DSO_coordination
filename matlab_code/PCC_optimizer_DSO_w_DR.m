function [result_DSO_lookahead,DA_market_outcome] = PCC_optimizer_DSO_w_DR(data_tso,dual_DA_gen,dual_DA_dem,...
	dual_DA_DR,cost_RT,iter,p_gen_DA_hat,p_dem_DA_hat,p_DR_DA_hat,npart_in)


%% load parameters
VOLL = 1000;
offer_wind = 0;
alpha_min = -300;

cc{1} = double.empty(1,0);
[nb,~,busN,~,~,~,Vmax,Vmin,~,~,~,...
	ng,nd,genL,~,incidentG,demL,incidentD,incidentW,Qmax_gen,Qmin_gen,...
    Pmax_gen,Pmin_gen,Qmax_dem,Qmin_dem,Pmax_dem,Pmin_dem,~,...
	~,~,fbusL,tbusL,SlmMax,~,~,~,~,...
    ~,~,~,~,~,~,offer_gen_DA,offer_gen_upreg,...
	offer_gen_downreg,bid_dem_DA,bid_dem_upreg,bid_dem_downreg,...
	nr_vals_cost,~,cost_type,~,areas,bus_areacodes,....
	incident_singlephase,~,rBranch,xBranch,~,prob_wscen,Wmax,Wmin,...
	~,nscen,Wmax_mean_DA,bus_wgen] = Data_Reader(data_tso,cc);

if ~isempty(data_tso.DR_DSO)
	bus_DR = data_tso.DR_DSO(:,1);
	Pmax_DR = data_tso.DR_DSO(:,2);
	Pmin_DR = data_tso.DR_DSO(:,3);
	bid_DA_DR = data_tso.DR_DSO(:,4);
	bid_up_DR = data_tso.DR_DSO(:,5);
	bid_dn_DR = data_tso.DR_DSO(:,6);

	DR_N = busN(bus_DR);

	nDR = length(bus_DR);
	incidentDR =  sparse(1:nDR, DR_N, 1 , nDR, nb);
	incidentDR = full(incidentDR);
else
	bus_DR = [];
	Pmax_DR = [];
	Pmin_DR = [];
	bid_DA_DR = [];
	bid_up_DR = [];
	bid_dn_DR = [];

	DR_N = [];

	nDR = length(bus_DR);
	incidentDR =  sparse(1:nDR, DR_N, 1 , nDR, nb);
	incidentDR = full(incidentDR);
end

cvx_solver gurobi
cvx_solver_settings('TimeLimit',300)


[~, n_areas, DSO_codes, fbusL_area, tbusL_area, branch_num_area, ...
	~, ~, ~, overlap,...
	~, ~, ~, ~,~,~,~] ...
	= find_overlap(bus_areacodes, fbusL, tbusL, areas);

ne = n_areas-1;

% [~,~,windg_area_T] = intersect(areas{1},bus_wgen);
[~,~,gens_T] = intersect(areas{1},genL);
[~,~,dems_T] = intersect(areas{1},demL);


nw = length(bus_wgen);
gens_DSO = [];
dems_DSO = [];
for DSO_num = 1:length(DSO_codes)
	[~,~,gens_area{DSO_num}] = intersect(areas{DSO_num+1},genL);
	[~,~,dems_area{DSO_num}] = intersect(areas{DSO_num+1},demL);
	[~,~,windg_area{DSO_num}] = intersect(areas{DSO_num+1},bus_wgen);
	
	gens_DSO = [gens_DSO; gens_area{DSO_num}];
	dems_DSO = [dems_DSO; dems_area{DSO_num}];
end
npart = npart_in;
npart2 = npart_in;
%% define variables
cvx_begin %quiet

% cvx_precision high

	%% PCC optimizer variables
	variable pi_pcc_DA(ne)
	variable fe_up(ne)
	variable fe_dn(ne)
% 	variable pi_pcc_up(ne)
% 	variable pi_pcc_dn(ne)
	
%% DSO market variables
	variable p_gen_DA_tilde(ng)
	variable p_dem_DA_tilde(nd)
	variable p_DR_DA(nDR)

	variable wind_DA_E(ne) nonnegative
	variable shed_p_DA_E(ne) nonnegative
	variable p_pcc_DA_E(ne)

	expression cost_DR_DA(nDR)
	
	%% dual variables
	
	variable lambda_DA_e(ne)
	variable varsigma_gen_DA_minus_e(ng) nonnegative
	variable varsigma_gen_DA_plus_e(ng) nonnegative
	variable varsigma_dem_DA_minus_e(nd) nonnegative
	variable varsigma_dem_DA_plus_e(nd) nonnegative
	variable varsigma_DR_DA_minus_e(nDR) nonnegative
	variable varsigma_DR_DA_plus_e(nDR) nonnegative
	variable iota_minus_DA(ne) nonnegative
	variable iota_plus_DA(ne) nonnegative
	variable rho_DA_minus(ne) nonnegative
	variable rho_DA_plus(ne) nonnegative
	variable Upsilon_DA_minus(ne) nonnegative
	variable Upsilon_DA_plus(ne) nonnegative
	

	variable bin_pcclim_DA(ne,2) binary
	variable bin_windlim_DA(ne,2) binary
	variable bin_dem(nd,2) binary 
	variable bin_gen(ng,2) binary
	
	expression dual_obj_linecap(ne,nscen)
	%% McCormick Variables for the linear envelopes
	variable w_piDA_pDA(ne)
	variable w_fdn_rhominus_DA(ne)
	variable w_fup_rhoplus_DA(ne)

	variable fdn_part(ne,npart2) binary
    variable rho_DA_plus_part(ne,npart2) binary
    
    variable rho_DA_minus_hat(ne,npart2)
    variable fup_hat(ne,npart2)

    variable p_pcc_DA_hat(ne,npart)
    
    variable pi_pcc_DA_part(ne,npart) binary
%% DA market variables on TSO level


variable varsigma_g_plus(ng)
variable varsigma_g_minus(ng)

variable sigma_g_plus(ng)
variable sigma_g_minus(ng)

variable varsigma_d_plus(nd)
variable varsigma_d_minus(nd)

variable sigma_d_plus(nd)
variable sigma_d_minus(nd)

variable lambda_DA

variable nu_minus
variable nu_plus

variable rho_minus
variable rho_plus

variable p_gen_DA(ng)
variable p_dem_DA(nd)
variable wind_DA
variable shed_p

variable alpha_cut(nscen)

variable wind_comp(2) binary
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

cost_shed = VOLL*shed_p;
cost_wind = wind_DA*offer_wind;

for m = 1:nDR
	cost_DR_DA(m) = bid_DA_DR(m)*p_DR_DA(m);
end

cost_DA = sum(cost_gen)-sum(cost_dem) + cost_shed + cost_wind - sum(cost_DR_DA);

cost_cuts = prob_wscen(1,:) * alpha_cut;
cost = cost_DA + cost_cuts;
%% objective statement
minimize(cost)

subject to
%% stationary constraints

for g = 1:ng
	if any(g == gens_DSO)
		offer_gen_DA(g,2) - (varsigma_g_minus(g) - varsigma_g_plus(g)) - lambda_DA == 0;
	elseif any(g == gens_T)
		offer_gen_DA(g,2) - (sigma_g_minus(g) - sigma_g_plus(g)) - lambda_DA == 0;
	end
end

for d = 1:nd
	if any(d == dems_DSO)
		-bid_dem_DA(d) - (varsigma_d_minus(d) - varsigma_d_plus(d)) + lambda_DA - rho_plus == 0;
	elseif any(d == dems_T)
		 -bid_dem_DA(d) - (sigma_d_minus(d) - sigma_d_plus(d)) + lambda_DA - rho_plus == 0;
	end
end

VOLL - lambda_DA - rho_minus + rho_plus == 0;

offer_wind - lambda_DA - nu_minus + nu_plus == 0;

%% complimentarity constraints
for g = 1:ng
	if any(g == gens_DSO)		
		varsigma_g_minus(g) <= VOLL * gen_comp(g,1);
		(p_gen_DA(g) - Pmin_gen(g)) <= VOLL *  (1 - gen_comp(g,1));

		varsigma_g_plus(g)  <= VOLL * gen_comp(g,2);
		p_gen_DA_tilde(g) - p_gen_DA(g) <= VOLL * (1 - gen_comp(g,2));
	elseif any(g == gens_T)
		sigma_g_minus(g) <= VOLL * gen_comp(g,3);
		(p_gen_DA(g) - Pmin_gen(g)) <= VOLL *  (1 - gen_comp(g,3));

		sigma_g_plus(g)  <= VOLL * gen_comp(g,4);
		Pmax_gen(g) - p_gen_DA(g) <= VOLL * (1 - gen_comp(g,4));
	end
end

for d = 1:nd
	if any(d == dems_DSO)		
		varsigma_d_minus(d) <= VOLL * dem_comp(d,1);
		(p_dem_DA(d) - Pmin_dem(d)) <= VOLL *  (1 - dem_comp(d,1));

		varsigma_d_plus(d)  <= VOLL * dem_comp(d,2);
		p_dem_DA_tilde(d) - p_dem_DA(d) <= VOLL * (1 - dem_comp(d,2));
	elseif any(d == dems_T)
		sigma_d_minus(d) <= VOLL * dem_comp(d,3);
		(p_dem_DA(d) - Pmin_dem(d)) <= VOLL *  (1 - dem_comp(d,3));

		sigma_d_plus(d)  <= VOLL * dem_comp(d,4);
		Pmax_dem(d) - p_dem_DA(d) <= VOLL * (1 - dem_comp(d,4));
	end
end

nu_minus <= VOLL *wind_comp(1);
wind_DA <= VOLL * (1 - wind_comp(1));

nu_plus <= VOLL * wind_comp(2);
sum(Wmax_mean_DA) - wind_DA <= VOLL * (1 - wind_comp(2));

rho_minus <= 1.5*VOLL * shed_comp(1);
shed_p <= VOLL * (1 - shed_comp(1));

rho_plus <= VOLL * shed_comp(2);
sum(p_dem_DA) - shed_p <= VOLL * (1 - shed_comp(2));

%% primal constrains

sum(p_gen_DA) + wind_DA + shed_p == sum(p_dem_DA);


for DSO_num = 1:length(DSO_codes)
	for k = 1:ng
		if any(k == gens_area{DSO_num})			
			Pmin_gen(k) <= p_gen_DA(k) <= p_gen_DA_tilde(k);
		else
			Pmin_gen(k) <= p_gen_DA(k) <= Pmax_gen(k);
		end
	end
	
	for k = 1:nd
		if any(k == dems_area{DSO_num})
			Pmin_dem(k) <= p_dem_DA(k) <= p_dem_DA_tilde(k);
		else
			Pmin_dem(k) <= p_dem_DA(k) <= Pmax_dem(k);
		end
	end
	
end

0 <= wind_DA <= sum(Wmax_mean_DA);
0 <= shed_p <= sum(p_dem_DA);

p_gen_DA_tilde <= 20;
p_dem_DA_tilde <= 20;

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
			+ sum(dual_DA_DR{k}(:,s) .* (p_DR_DA - p_DR_DA_hat{k} ) );
	end
end


%% all DSO lookahead constraints
	for e = 1:ne
		[~,~,gens_area] = intersect(areas{e+1},genL);
		[~,~,dems_area] = intersect(areas{e+1},demL);
		[~,~,windg_area] = intersect(areas{e+1},bus_wgen);
		[~,~,DR_area] = intersect(areas{e+1},bus_DR);
		
		ng_d = length(gens_area);
		nd_d = length(dems_area);
		nDR_d = length(DR_area);

		
		
		%% day ahead cost 
% 		expression cost_gen_DA(ng)
% 		expression cost_dem_DA(nd)
% 		for m = 1:ng_d
% 			if cost_type(m) == 2
% 				for k = 1:nr_vals_cost
% 					cost_gen_DA(gens_area(m)) = cost_gen_DA(gens_area(m)) + sum(offer_gen_DA(gens_area(m),k).*p_gen_DA_tilde(gens_area(m)).^(nr_vals_cost-k));
% 				end
% 			elseif cost_type(m) == 1
% 				error('Wrong cost function, piecewise linear is not convex');
% 			end
% 		end
% 
% 		for m = 1:nd_d
% 			cost_dem_DA(dems_area(m)) = bid_dem_DA(dems_area(m))*p_dem_DA_tilde(dems_area(m));
% 		end
% 
% 		cost_shed_DA = VOLL*shed_p_DA_E(e);
% 		cost_wind_DA = wind_DA_E(e)*offer_wind;
% % 		cost_PCC_DA = pi_pcc_DA(e) * p_pcc_DA_E(e);
% 		cost_PCC_DA = w_piDA_pDA(e);
% 
% 
% 		DSO_cost_DA(e) = sum(cost_gen_DA(gens_area) )-sum(cost_dem_DA(dems_area)) - sum(cost_DR_DA(DR_area))...
% 			+ cost_shed_DA + cost_wind_DA + cost_PCC_DA;
% 		
%         %% cost function statement
% % 		DSO_cost(e) = cost_RT_expec(e) + DSO_cost_DA(e);
% 		DSO_cost(e) = DSO_cost_DA(e);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		%% McCormick expansion of bi-linear terms


% 		sum(pi_pcc_DA_part(e,:)) == 1;
% 		
% 		
% 		pi_pcc_DA_lo = 5;
% 		pi_pcc_DA_up = 51;
% 		
% 		
% 		p_pcc_DA_E_lo = -10;
% 		p_pcc_DA_E_up = 10;
% 	
% 		pi_pcc_DA_lo_part = zeros(npart,1);
% 		pi_pcc_DA_up_part = zeros(npart,1);
% 
% 		
% 		for pl = 1:npart			
% 			pi_pcc_DA_lo_part(pl) = pi_pcc_DA_lo - (pi_pcc_DA_lo - pi_pcc_DA_up)*(pl - 1)/npart;
% 			pi_pcc_DA_up_part(pl) = pi_pcc_DA_lo + (pi_pcc_DA_up - pi_pcc_DA_lo)*(pl)/npart;
%         end
% 		
% 		
% 		p_pcc_DA_E_lo * pi_pcc_DA_part(e,:) <= p_pcc_DA_hat(e,:) <= p_pcc_DA_E_up * pi_pcc_DA_part(e,:);
% 		sum( pi_pcc_DA_lo_part(:) .* pi_pcc_DA_part(e,:)') <= pi_pcc_DA(e) <= sum( pi_pcc_DA_up_part(:) .* pi_pcc_DA_part(e,:)');
% 		p_pcc_DA_E(e) == sum(p_pcc_DA_hat(e,:));
% 			
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% 		w_piDA_pDA(e) >= sum( pi_pcc_DA_lo_part(:) .* p_pcc_DA_hat(e,:)' ...
% 			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_lo_part(:) * p_pcc_DA_E_lo ) + pi_pcc_DA(e) * p_pcc_DA_E_lo;
% 		w_piDA_pDA(e) >= sum( pi_pcc_DA_up_part(:) .* p_pcc_DA_hat(e,:)' ...
% 			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_up_part(:) * p_pcc_DA_E_up ) + pi_pcc_DA(e) * p_pcc_DA_E_up;
% 		w_piDA_pDA(e) <= sum( pi_pcc_DA_up_part(:) .* p_pcc_DA_hat(e,:)'...
% 			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_up_part(:) * p_pcc_DA_E_lo ) + pi_pcc_DA(e) * p_pcc_DA_E_lo;
% 		w_piDA_pDA(e) <= sum(pi_pcc_DA_lo_part(:) .* p_pcc_DA_hat(e,:)' ...
% 			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_lo_part(:) * p_pcc_DA_E_up ) + pi_pcc_DA(e) * p_pcc_DA_E_up;
% 		
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 		sum(fdn_part(e,:)) == 1;
% 		sum(rho_DA_plus_part(e,:)) == 1;
% 				
% 		fdn_lo = -10;
% 		fdn_up = 5;
% 		
% 		rhominus_lo = 0;
% 		rhominus_up = 25;
% 		
% 		fup_lo = -5;
% 		fup_up = 10;
% 		
% 		rhoplus_lo = 0;
% 		rhoplus_up = 25;
% 		
% 		fdn_lo_part = zeros(npart2,1);
% 		fdn_up_part = zeros(npart2,1);
% 		fup_lo_part = zeros(npart2,1);
% 		fup_up_part = zeros(npart2,1);
% 		rhominus_lo_part = zeros(npart2,1);
% 		rhominus_up_part = zeros(npart2,1);
% 		rhoplus_lo_part = zeros(npart2,1);
% 		rhoplus_up_part = zeros(npart2,1);
% 		
% 		for pl = 1:npart2
% 			fdn_lo_part(pl) = fdn_lo - (fdn_lo - fdn_up)*(pl - 1)/npart2;
% 			fdn_up_part(pl) = fdn_lo + (fdn_up - fdn_lo)*(pl)/npart2;
% 			
% 			fup_lo_part(pl) = fup_lo - (fup_lo - fup_up)*(pl - 1)/npart2;
% 			fup_up_part(pl) = fup_lo + (fup_up - fup_lo)*(pl)/npart2;
% 			
% 			rhominus_lo_part(pl) = rhominus_lo - (rhominus_lo - rhominus_up)*(pl - 1)/npart2;
% 			rhominus_up_part(pl) = rhominus_lo + (rhominus_up - rhominus_lo)*(pl)/npart2;
% 			
% 			rhoplus_lo_part(pl) = rhoplus_lo - (rhoplus_lo - rhoplus_up)*(pl - 1)/npart2;
% 			rhoplus_up_part(pl) = rhoplus_lo + (rhoplus_up - rhoplus_lo)*(pl)/npart2;
% 		end
% 		
% 		rhominus_lo * fdn_part(e,:) <= rho_DA_minus_hat(e,:) <= rhominus_up * fdn_part(e,:);
% 		sum( fdn_lo_part(:) .* fdn_part(e,:)') <= fe_dn(e) <= sum( fdn_up_part(:) .* fdn_part(e,:)');
% 		rho_DA_minus(e) == sum(rho_DA_minus_hat(e,:));
% 		
% 		fup_lo * rho_DA_plus_part(e,:) <= fup_hat(e,:) <= fup_up * rho_DA_plus_part(e,:);
% 		sum( rhoplus_lo_part(:) .* rho_DA_plus_part(e,:)') <= rho_DA_plus(e) <= sum( rhoplus_up_part(:) .* rho_DA_plus_part(e,:)');
% 		fe_up(e) == sum(fup_hat(e,:));				
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		w_fdn_rhominus_DA(e) >= sum( fdn_lo_part(:) .* rho_DA_minus_hat(e,:)' ...
% 			- fdn_part(e,:)' .* fdn_lo_part(:) * rhominus_lo ) + fe_dn(e) * rhominus_lo;
% 		w_fdn_rhominus_DA(e) >= sum( fdn_up_part(:) .* rho_DA_minus_hat(e,:)' ...
% 			- fdn_part(e,:)' .* fdn_up_part(:) * rhominus_up ) + fe_dn(e) * rhominus_up ;
% 		w_fdn_rhominus_DA(e) <= sum( fdn_up_part(:) .* rho_DA_minus_hat(e,:)'  ...
% 			- fdn_part(e,:)' .* fdn_up_part(:) * rhominus_lo ) + fe_dn(e) * rhominus_lo;
% 		w_fdn_rhominus_DA(e) <=  fe_dn(e) * rhominus_up + sum( fdn_lo_part(:) .* rho_DA_minus_hat(e,:)' ...
% 			- fdn_part(e,:)' .* fdn_lo_part(:) * rhominus_up );
% 		
% 		w_fup_rhoplus_DA(e) >= fup_lo * rho_DA_plus(e) + sum( fup_hat(e,:)' .* rhoplus_lo_part(:) ...
% 			- rho_DA_plus_part(e,:)' .* fup_lo .* rhoplus_lo_part(:) );
% 		w_fup_rhoplus_DA(e) >= fup_up * rho_DA_plus(e) + sum( fup_hat(e,:)' .* rhoplus_up_part(:) ...
% 			- rho_DA_plus_part(e,:)' .* fup_up .* rhoplus_up_part(:) );
% 		w_fup_rhoplus_DA(e) <= fup_up * rho_DA_plus(e) + sum( fup_hat(e,:)' .* rhoplus_lo_part(:) ...
% 			- rho_DA_plus_part(e,:)' .* fup_up .* rhoplus_lo_part(:) );
% 		w_fup_rhoplus_DA(e) <= sum( fup_hat(e,:)' .* rhoplus_up_part(:) ...
% 			- rho_DA_plus_part(e,:)' .* fup_lo .* rhoplus_up_part(:) ) + fup_lo * rho_DA_plus(e);		
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%% primal constraints	
% 		subject to
		%% day-ahead primal constraints
		sum(p_gen_DA_tilde(gens_area)) + wind_DA_E(e) + shed_p_DA_E(e) + p_pcc_DA_E(e) == sum(p_dem_DA_tilde(dems_area)) + sum(p_DR_DA(DR_area));


		Pmin_gen(gens_area) <= p_gen_DA_tilde(gens_area) <= Pmax_gen(gens_area);
		Pmin_dem(dems_area) <= p_dem_DA_tilde(dems_area) <= Pmax_dem(dems_area);
		0 <= wind_DA_E(e) <= sum(Wmax_mean_DA(windg_area));
		0 <= shed_p_DA_E(e) <= sum(p_dem_DA_tilde(dems_area));

		fe_dn(e) <= p_pcc_DA_E(e) <= fe_up(e);
		
		Pmin_DR(DR_area) <= p_DR_DA(DR_area) <= Pmax_DR(DR_area);

		
		%% KKT stationary conditions
		for g = 1:ng_d
			offer_gen_DA(gens_area(g),2) - lambda_DA_e(e) ...
				- varsigma_gen_DA_minus_e(gens_area(g)) + varsigma_gen_DA_plus_e(gens_area(g))  == 0; % p_gen_DA
		end

		for d = 1:nd_d
			-bid_dem_DA(dems_area(d)) + lambda_DA_e(e) - varsigma_dem_DA_minus_e(dems_area(d))  ...
				+ varsigma_dem_DA_plus_e(dems_area(d)) - Upsilon_DA_plus(e) == 0; % p_dem_DA
		end
		
		for d = 1:nDR_d
			-bid_DA_DR(DR_area(d)) + lambda_DA_e(e)...
				- varsigma_DR_DA_minus_e(DR_area(d)) + varsigma_DR_DA_plus_e(DR_area(d)) == 0; % p_DR_DA
		end


		pi_pcc_DA(e)  - lambda_DA_e(e) - rho_DA_minus(e)...
			+ rho_DA_plus(e) == 0; % p_pcc_DA

		VOLL - lambda_DA_e(e) - Upsilon_DA_minus(e) + Upsilon_DA_plus(e) == 0; % shed_p_DA
		offer_wind - lambda_DA_e(e) - iota_minus_DA(e) + iota_plus_DA(e) == 0; % wind_DA
		
		
		%% DSO compimentarity constraints
	for g = 1:ng_d
		varsigma_gen_DA_minus_e(gens_area(g)) <= VOLL * bin_gen(gens_area(g),1);
		p_gen_DA_tilde(gens_area(g)) - Pmin_gen(gens_area(g)) <= VOLL *  (1 - bin_gen(gens_area(g),1));

		varsigma_gen_DA_plus_e(gens_area(g)) <= VOLL * bin_gen(gens_area(g),2);
		Pmax_gen(gens_area(g)) - p_gen_DA_tilde(gens_area(g)) <= VOLL *  (1 - bin_gen(gens_area(g),2));
		
	end
	
	
	for d = 1:nd_d
		varsigma_dem_DA_minus_e(dems_area(d)) <= VOLL * bin_dem(dems_area(d),1);
		p_dem_DA_tilde(dems_area(d)) - Pmin_dem(dems_area(d)) <= VOLL *  (1 - bin_dem(dems_area(d),1));

		varsigma_dem_DA_plus_e(dems_area(d)) <= VOLL * bin_dem(dems_area(d),2);
		Pmax_dem(dems_area(d)) - p_dem_DA_tilde(dems_area(d)) <= VOLL *  (1 - bin_dem(d,2));
	end
 	
	iota_minus_DA(e) <= VOLL * bin_windlim_DA(e,1);
	wind_DA_E(e) <= VOLL * (1 - bin_windlim_DA(e,1));
	
	iota_plus_DA(e) <= VOLL * bin_windlim_DA(e,2);
	sum(Wmax_mean_DA(windg_area)) - wind_DA_E(e) <= VOLL * (1 - bin_windlim_DA(e,2));
	
	rho_DA_minus(e) <= VOLL * bin_pcclim_DA(e,1);
	p_pcc_DA_E(e) - fe_dn(e) <= VOLL * (1 - bin_pcclim_DA(e,1));
	rho_DA_plus(e) <= VOLL * bin_pcclim_DA(e,2);
	fe_up(e) - p_pcc_DA_E(e) <= VOLL * (1 - bin_pcclim_DA(e,2));


% 		%% strong duality objective
% 
% 		dual_objective_DA(e) = sum(Pmin_gen(gens_area)' * varsigma_gen_DA_minus_e(gens_area)...
% 			- Pmax_gen(gens_area)' * varsigma_gen_DA_plus_e(gens_area) ) ...
% 			+ sum(Pmin_dem(dems_area)' * varsigma_dem_DA_minus_e(dems_area) ...
% 			- Pmax_dem(dems_area)' * varsigma_dem_DA_plus_e(dems_area) )...
% 			+ sum(Pmin_DR(DR_area)' * varsigma_DR_DA_minus_e(DR_area) ...
% 			- Pmax_DR(DR_area)' * varsigma_DR_DA_plus_e(DR_area) )...
% 			- iota_plus_DA(e) * sum(Wmax_mean_DA(windg_area)) ...
% 			+ w_fdn_rhominus_DA(e) - w_fup_rhoplus_DA(e);
% % 			+ fe_dn(e) * rho_DA_minus(e) - fe_up(e) * rho_DA_plus(e);
% 
% 		dual_objective(e) = dual_objective_DA(e);% + sum(dual_objective_RT(e,:)) + sum(dual_obj_linecap(e,:));
% 
% 		DSO_cost(e) == dual_objective(e);			
	end
cvx_end

%% output

DA_market_outcome.p_gen_DA = p_gen_DA;
DA_market_outcome.p_dem_DA = p_dem_DA;
DA_market_outcome.cost_DA = cost_DA;
DA_market_outcome.lambda = lambda_DA;
DA_market_outcome.cost_cuts = cost_cuts;
DA_market_outcome.cost = cost;
DA_market_outcome.alpha_cut = alpha_cut;
DA_market_outcome.wind_DA = wind_DA;
DA_market_outcome.shed_p = shed_p;
DA_market_outcome.p_gen_DA_tilde = p_gen_DA_tilde;
DA_market_outcome.p_dem_DA_tilde = p_dem_DA_tilde;


if DA_market_outcome.shed_p > 10e-5
	warning('There is load shedding in the Day-Ahead market');
	pause(3)
end

% result_DSO_lookahead.DSO_cost = DSO_cost;
result_DSO_lookahead.cost_DR_DA = cost_DR_DA;
result_DSO_lookahead.p_gen_DA_tilde = p_gen_DA_tilde;
result_DSO_lookahead.p_dem_DA_tilde = p_dem_DA_tilde;
result_DSO_lookahead.p_DR_DA = p_DR_DA;
result_DSO_lookahead.wind_DA_E = wind_DA_E;
result_DSO_lookahead.p_pcc_DA_E = p_pcc_DA_E;
result_DSO_lookahead.fe_up = fe_up;
result_DSO_lookahead.fe_dn = fe_dn;
result_DSO_lookahead.pi_pcc_DA = pi_pcc_DA;
end