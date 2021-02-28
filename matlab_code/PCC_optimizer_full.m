function [result_DSO_lookahead,DA_market_outcome] = PCC_optimizer_full(data_tso,dual_DA_gen,dual_DA_dem,...
	dual_DA_DR,cost_RT,iter,p_gen_DA_hat,p_dem_DA_hat,p_DR_DA_hat,npart_in)


%% load parameters
VOLL = 1000;
offer_wind = 0;
alpha_min = -200;

cc{1} = double.empty(1,0);
[nb,~,busN,~,~,~,Vmax,Vmin,~,~,~,...
	ng,nd,genL,~,incidentG,demL,incidentD,incidentW,Qmax_gen,Qmin_gen,...
    Pmax_gen,Pmin_gen,Qmax_dem,Qmin_dem,Pmax_dem,Pmin_dem,~,...
	~,~,fbusL,tbusL,SlmMax,~,~,~,~,...
    ~,~,~,~,~,~,offer_gen_DA,offer_gen_upreg,...
	offer_gen_downreg,bid_dem_DA,bid_dem_upreg,bid_dem_downreg,...
	nr_vals_cost,~,cost_type,~,areas,bus_areacodes,....
	incident_singlephase,~,rBranch,xBranch,~,prob_wscen,Wmax,Wmin,...
	~,nscen,Wmax_mean_DA,windgL,...
	offer_wind_DA,offer_wind_up,offer_wind_dn] = Data_Reader(data_tso,cc);

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

windg_area_T = find(ismember(windgL,areas{1}));
gens_T = find(ismember(genL,areas{1}));
dems_T = find(ismember(demL,areas{1}));

nw = length(windgL);
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
npart = npart_in;
npart2 = npart_in;
%% define variables
cvx_begin %quiet

% cvx_precision high

	%% PCC optimizer variables
	variable pi_pcc_DA(ne)
	variable fe_up(ne)
	variable fe_dn(ne)
	variable pi_pcc_up(ne)
	variable pi_pcc_dn(ne)
	
%% DSO lookahead variables
	variable p_gen_DA_tilde(ng)
	variable p_dem_DA_tilde(nd)
	variable p_DR_DA(nDR)

% 	variable p_gen_DA(ng)
% 	variable p_dem_DA(nd) nonnegative
	variable wind_DA_E(ne) nonnegative
	variable shed_p_DA_E(ne) nonnegative
	variable p_pcc_DA_E(ne)

	
	
	variable p_DR_RT(nDR,nscen)
	variable pup_DR(nDR,nscen) nonnegative
	variable pdn_DR(nDR,nscen) nonnegative
	variable pup_g(ng,nscen) nonnegative
	variable pdn_g(ng,nscen) nonnegative
	variable pup_d(nd,nscen) nonnegative
	variable pdn_d(nd,nscen) nonnegative
	variable p_pcc_rt(ne,nscen)
	variable q_pcc_rt(ne,nscen)
	variable pup_pcc(ne,nscen) nonnegative
	variable pdn_pcc(ne,nscen) nonnegative
	variable shed_p_RT(nb,nscen) nonnegative
	variable shed_q_RT(nb,nscen) nonnegative
	variable wind_RT(nw,nscen) nonnegative

	variable p_inj_rt(nb,nscen)
	variable q_inj_rt(nb,nscen)
	variable p_flow(nb,nb,nscen)
	variable q_flow(nb,nb,nscen)
	variable p_gen_RT(ng,nscen)
	variable p_dem_RT(nd,nscen)
	variable q_gen_RT(ng,nscen)
	variable q_dem_RT(nd,nscen)

	variable i_sq(nb,nb,nscen) nonnegative
	variable v_sq(nb,nscen) nonnegative
    
    variable p_pcc_comp(ne,nscen) binary
	variable p_gen_comp(ng,nscen) binary
	variable p_dem_comp(nd,nscen) binary
	variable p_DR_comp(nDR,nscen) binary

	expression cost_g_RT(ng,nscen)
	expression cost_d_RT(nd,nscen)
	expression cost_DR_RT(nDR,nscen)
	expression cost_DR_DA(nDR)
% 	expression cost_RT
	
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
	
	
	variable zeta_p_gen_RT(ng,nscen)
	variable zeta_p_dem_RT(nd,nscen)
	variable zeta_DR_RT(nDR,nscen)
	variable zeta_pcc_RT(ne,nscen)
	variable lambda_p_RT(nb,nscen)
	variable lambda_q_RT(nb,nscen)
	
	variable gamma_RT_p(nb,nb,nscen)
	variable gamma_RT_q(nb,nb,nscen)
	variable gamma_RT_l(nb,nb,nscen)
	variable gamma_RT_w(nb,nb,nscen) nonnegative
	
	variable mu_p_RT(nb,nb,nscen)
	variable mu_q_RT(nb,nb,nscen)
	
	variable eta_RT_p(nb,nb,nscen) 
	variable eta_RT_q(nb,nb,nscen) 
	variable eta_RT_s(nb,nb,nscen) nonnegative
	
	variable beta_RT(nb,nb,nscen)
	variable sigma_RT_minus(nb,nscen) nonnegative
	variable sigma_RT_plus(nb,nscen) nonnegative
	variable nu_RT_minus(nw,nscen) nonnegative
	variable nu_RT_plus(nw,nscen) nonnegative
	variable varsigma_gen_RT_minus(ng,nscen) nonnegative
	variable varsigma_gen_RT_plus(ng,nscen) nonnegative
	variable varsigma_dem_RT_minus(nd,nscen) nonnegative
	variable varsigma_dem_RT_plus(nd,nscen) nonnegative
	variable varsigma_DR_RT_minus(nDR,nscen) nonnegative
	variable varsigma_DR_RT_plus(nDR,nscen) nonnegative
	variable kappa_gen_RT_minus(ng,nscen) nonnegative
	variable kappa_gen_RT_plus(ng,nscen) nonnegative
	variable kappa_dem_RT_minus(nd,nscen) nonnegative
	variable kappa_dem_RT_plus(nd,nscen) nonnegative
	variable rho_RT_minus(ne,nscen) nonnegative
	variable rho_RT_plus(ne,nscen) nonnegative
	
	variable epsilon_up_RT(ng,nscen) nonnegative
	variable epsilon_dn_RT(ng,nscen) nonnegative
	variable varepsilon_up_RT(nd,nscen) nonnegative
	variable varepsilon_dn_RT(nd,nscen) nonnegative
	variable epsilon_DR_up_RT(nDR,nscen) nonnegative
	variable epsilon_DR_dn_RT(nDR,nscen) nonnegative
	variable epsilon_pcc_up_RT(ne,nscen) nonnegative
	variable epsilon_pcc_dn_RT(ne,nscen) nonnegative
	
	variable Upsilon_RT_minus(nb,nscen) nonnegative
	variable Upsilon_RT_plus(nb,nscen) nonnegative
	
	variable Upsilon_Q_RT_minus(nb,nscen) nonnegative
	variable Upsilon_Q_RT_plus(nb,nscen) nonnegative
	
	variable rho_DA_comp(ne) binary
	variable rho_RT_comp(ne,nscen) binary
	
	variable varsigma_gen_DA_comp(ng) binary
	variable varsigma_dem_DA_comp(nd) binary
	
	variable varsigma_gen_RT_comp(ng,nscen) binary
	variable varsigma_dem_RT_comp(nd,nscen) binary
	variable varsigma_DR_RT_comp(nDR,nscen) binary
	variable sigma_RT_comp(nb,nscen) binary
	variable Upsilon_RT_comp(nb,nscen) binary
	variable Upsilon_DA_comp(ne) binary
	
	expression dual_obj_linecap(ne,nscen)
	%% McCormick Variables for the linear envelopes
	variable w_piDA_pRT(ne,nscen)
	variable w_piDA_pDA(ne)
	variable w_piup_pup(ne,nscen)
	variable w_pidn_pdn(ne,nscen)

	variable w_fdn_rhominus_DA(ne)
	variable w_fup_rhoplus_DA(ne)
	variable w_fdn_rhominus_RT(ne,nscen)
	variable w_fup_rhoplus_RT(ne,nscen)
    
    variable fdn_part(ne,npart2) binary
    variable rho_DA_plus_part(ne,npart2) binary
    variable rho_RT_minus_part(ne,npart2,nscen) binary
    variable fup_part(ne,npart2,nscen) binary
    
    variable rho_DA_minus_hat(ne,npart2)
    variable fup_hat(ne,npart2)
    variable fdn_hat(ne,npart2,nscen)
    variable rho_RT_plus_hat(ne,npart2,nscen)

    variable pi_pcc_DA_hat(ne,npart,nscen)
    variable p_pcc_DA_hat(ne,npart)
    variable pup_pcc_hat(ne,npart,nscen)
    variable pi_pcc_dn_hat(ne,npart,nscen)
    
    variable p_pcc_rt_part(ne,npart,nscen) binary
    variable pi_pcc_DA_part(ne,npart) binary
    variable pi_pcc_up_part(ne,npart,nscen) binary
    variable pdn_pcc_part(ne,npart,nscen) binary
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
cost = cost_DA + cost_cuts ;
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

VOLL - lambda_DA - rho_minus + rho_plus == 0;

offer_wind - lambda_DA - nu_minus + nu_plus == 0;

%% complimentarity constraints
for g = 1:ng
	if any(g == gens_DSO)
% 		gc = find(gens_DSO == g);
		
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
% 		dc = find(dems_DSO == d);
		
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
% 			gc = find(gens_DSO == k);
			
			Pmin_gen(k) <= p_gen_DA(k) <= p_gen_DA_tilde(k);
		else
			Pmin_gen(k) <= p_gen_DA(k) <= Pmax_gen(k);
		end
	end
	
	for k = 1:nd
		if any(k == dems_area{DSO_num})
% 			dc = find(dems_DSO == k);

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

	incident_contra = ones(nb,nb) - incident_singlephase;
	incident_contra = logical(incident_contra);

	[n1,n2] = ind2sub([nb nb],find(incident_contra));
	for k = 1:length(n1)
		p_flow(n1(k),n2(k),:) == 0;
		q_flow(n1(k),n2(k),:) == 0;
	end
	
	buses_DSO = cat(1,areas{2:end});
	nb_DSO = length(buses_DSO);

	incident_contra = ones(nb,nb);
	incident_contra(buses_DSO,buses_DSO) = ones(nb_DSO,nb_DSO) - incident_singlephase(buses_DSO,buses_DSO);
	incident_contra = logical(incident_contra);
	[n3,n4] = ind2sub([nb nb],find(incident_contra));
	
% 	for s = 1:nscen
		for k = 1:length(n3)
			gamma_RT_l(n3(k),n4(k),:) == 0;
			gamma_RT_w(n3(k),n4(k),:) == 0;
			beta_RT(n3(k),n4(k),:) == 0;
			eta_RT_s(n3(k),n4(k),:) == 0;
		end
% 	end
	

	
	for e = 1:ne
		gens_area = find(ismember(genL,areas{e+1}));
		dems_area =  find(ismember(demL,areas{e+1}));
		windg_area = find(ismember(windgL,areas{e+1}));
		DR_area = find(ismember(bus_DR,areas{e+1}));

		buses = areas{e+1};
		branch_num = branch_num_area{e+1};
		fbus = fbusL_area{e+1};
		tbus = tbusL_area{e+1};

		
		ng_d = length(gens_area);
		nd_d = length(dems_area);
		nw_d = length(windg_area);
		nDR_d = length(DR_area);
		nb_d = length(buses);
		nl_d = length(fbus);
		
		pcc_global = intersect(overlap{1,e+1},areas{e+1});
		pcc_local = areas{e+1} == pcc_global;
		incident_pcc = zeros(nb_d,1);
		incident_pcc(pcc_local) = 1;
		
		%% Real time objective function
		for m = 1:ng_d
			cost_g_RT(gens_area(m),:) = offer_gen_DA(gens_area(m),2).*(p_gen_RT(gens_area(m),:) - p_gen_DA_tilde(gens_area(m))) ...
				+ offer_gen_upreg(gens_area(m))*pup_g(gens_area(m),:)...
				+ offer_gen_downreg(gens_area(m)) *pdn_g(gens_area(m),:);
		end

		for m = 1:nd_d
			cost_d_RT(dems_area(m),:) = bid_dem_DA(dems_area(m)).*(p_dem_DA_tilde(dems_area(m)) - p_dem_RT(dems_area(m),:)) ...
				+ bid_dem_upreg(dems_area(m)).*pup_d(dems_area(m),:) ...
				+ bid_dem_downreg(dems_area(m)) .* pdn_d(dems_area(m),:);
		end
		cost_shed_RT = VOLL*(shed_p_RT(buses,:)+shed_q_RT(buses,:));

		cost_wind_RT = wind_RT(windg_area,:)*offer_wind;

% 		shed_p_RT == 0;
% 		shed_q_RT == 0;
% 		cost_PCC_RT = pi_pcc_DA(e) * (p_pcc_rt(e,:) - p_pcc_DA_E(e)) ...
% 				+ pi_pcc_up(e) * pup_pcc(e,:)...
% 				+ pi_pcc_dn(e) * pdn_pcc(e,:);
			
		cost_PCC_RT(e,:) = w_piDA_pRT(e,:) - w_piDA_pDA(e) ...
			+ w_piup_pup(e,:) + w_pidn_pdn(e,:);
		
		for m = 1:nDR_d
			cost_DR_RT(DR_area(m),:) = bid_DA_DR(DR_area(m)).*(p_DR_DA(DR_area(m)) - p_DR_RT(DR_area(m),:)) ...
				+ bid_up_DR(DR_area(m)).*pup_DR(DR_area(m),:) ...
				+ bid_dn_DR(DR_area(m)) .* pdn_DR(DR_area(m),:);
		end

		DSO_cost_RT(:,e) = sum(cost_g_RT(gens_area,:),1) + sum(cost_d_RT(dems_area,:),1)...
			+ sum(cost_shed_RT,1) + sum(cost_wind_RT,1)+cost_PCC_RT(e,:)...
			+ sum(cost_DR_RT(DR_area,:),1);
		cost_RT_expec(e) = prob_wscen(1,:)*DSO_cost_RT(:,e);
		%% day ahead cost 
		expression cost_gen_DA(ng)
		expression cost_dem_DA(nd)
		for m = 1:ng_d
			if cost_type(m) == 2
				for k = 1:nr_vals_cost
					cost_gen_DA(gens_area(m)) = cost_gen_DA(gens_area(m)) + sum(offer_gen_DA(gens_area(m),k).*p_gen_DA_tilde(gens_area(m)).^(nr_vals_cost-k));
				end
			elseif cost_type(m) == 1
				error('Wrong cost function, piecewise linear is not convex');
			end
		end

		for m = 1:nd_d
			cost_dem_DA(dems_area(m)) = bid_dem_DA(dems_area(m))*p_dem_DA_tilde(dems_area(m));
		end

		cost_shed_DA = VOLL*shed_p_DA_E(e);
		cost_wind_DA = wind_DA_E(e)*offer_wind;
% 		cost_PCC_DA = pi_pcc_DA(e) * p_pcc_DA_E(e);
		cost_PCC_DA = w_piDA_pDA(e);
% 		shed_p_DA_E == 0;
		
		DSO_cost_DA(e) = sum(cost_gen_DA(gens_area) )-sum(cost_dem_DA(dems_area)) - sum(cost_DR_DA(DR_area))...
			+ cost_shed_DA + cost_wind_DA + cost_PCC_DA;
		
        %% cost function statement
		DSO_cost(e) = cost_RT_expec(e) + DSO_cost_DA(e);
% 		minimize(1)

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		%% McCormick expansion of bi-linear terms

		for s = 1:nscen
			sum(p_pcc_rt_part(e,:,s)) == 1;
			sum(pi_pcc_up_part(e,:,s)) == 1;
			sum(pdn_pcc_part(e,:,s)) == 1;
		end
		sum(pi_pcc_DA_part(e,:)) == 1;
		
		
		pi_pcc_DA_lo = 15;
		pi_pcc_DA_up = 42;
		
		p_pcc_rt_lo = -4;
		p_pcc_rt_up = 7;
		
		p_pcc_DA_E_lo = -4;
		p_pcc_DA_E_up = 7;
		
		pi_pcc_up_lo = 5;
		pi_pcc_up_up = 20;
		
		pup_pcc_lo = 0;
		pup_pcc_up = 6;
		
		pi_pcc_dn_lo = 5;
		pi_pcc_dn_up = 20;
		
		pdn_pcc_lo = 0;
		pdn_pcc_up = 6;
		
		p_pcc_rt_lo_part = zeros(npart,1);
		p_pcc_rt_up_part = zeros(npart,1);
		pi_pcc_DA_lo_part = zeros(npart,1);
		pi_pcc_DA_up_part = zeros(npart,1);
		pi_pcc_up_lo_part = zeros(npart,1);
		pi_pcc_up_up_part = zeros(npart,1);
		pdn_pcc_lo_part = zeros(npart,1);
		pdn_pcc_up_part = zeros(npart,1);
		
		for pl = 1:npart
			p_pcc_rt_lo_part(pl) = p_pcc_rt_lo - (p_pcc_rt_lo - p_pcc_rt_up)*(pl - 1)/npart;
			p_pcc_rt_up_part(pl) = p_pcc_rt_lo + (p_pcc_rt_up - p_pcc_rt_lo)*(pl)/npart;
			
			pi_pcc_DA_lo_part(pl) = pi_pcc_DA_lo - (pi_pcc_DA_lo - pi_pcc_DA_up)*(pl - 1)/npart;
			pi_pcc_DA_up_part(pl) = pi_pcc_DA_lo + (pi_pcc_DA_up - pi_pcc_DA_lo)*(pl)/npart;
			
			pi_pcc_up_lo_part(pl) = pi_pcc_up_lo - (pi_pcc_up_lo - pi_pcc_up_up)*(pl - 1)/npart;
			pi_pcc_up_up_part(pl) = pi_pcc_up_lo + (pi_pcc_up_up - pi_pcc_up_lo)*(pl)/npart;
			
			pdn_pcc_lo_part(pl) = pdn_pcc_lo - (pdn_pcc_lo - pdn_pcc_up)*(pl - 1)/npart;
			pdn_pcc_up_part(pl) = pdn_pcc_lo + (pdn_pcc_up - pdn_pcc_lo)*(pl)/npart;
        end
		
		for s = 1:nscen
			pi_pcc_DA_lo * p_pcc_rt_part(e,:,s) <= pi_pcc_DA_hat(e,:,s) <= pi_pcc_DA_up * p_pcc_rt_part(e,:,s);
			sum( p_pcc_rt_lo_part(:) .* p_pcc_rt_part(e,:,s)') <= p_pcc_rt(e,s) <= sum( p_pcc_rt_up_part(:) .* p_pcc_rt_part(e,:,s)');
			pi_pcc_DA(e) == sum(pi_pcc_DA_hat(e,:,s));
			
			pup_pcc_lo * pi_pcc_up_part(e,:,s) <= pup_pcc_hat(e,:,s) <= pup_pcc_up * pi_pcc_up_part(e,:,s);
			sum( pi_pcc_up_lo_part(:) .* pi_pcc_up_part(e,:,s)') <= pi_pcc_up(e) <= sum( pi_pcc_up_up_part(:) .* pi_pcc_up_part(e,:,s)');
			pup_pcc(e,s) == sum(pup_pcc_hat(e,:,s));
			
			pi_pcc_dn_lo * pdn_pcc_part(e,:,s) <= pi_pcc_dn_hat(e,:,s) <= pi_pcc_dn_up * pdn_pcc_part(e,:,s);
			sum( pdn_pcc_lo_part(:) .* pdn_pcc_part(e,:,s)') <= pdn_pcc(e,s) <= sum( pdn_pcc_up_part(:) .* pdn_pcc_part(e,:,s)');
			pi_pcc_dn(e) == sum(pi_pcc_dn_hat(e,:,s));
		end
		
		p_pcc_DA_E_lo * pi_pcc_DA_part(e,:) <= p_pcc_DA_hat(e,:) <= p_pcc_DA_E_up * pi_pcc_DA_part(e,:);
		sum( pi_pcc_DA_lo_part(:) .* pi_pcc_DA_part(e,:)') <= pi_pcc_DA(e) <= sum( pi_pcc_DA_up_part(:) .* pi_pcc_DA_part(e,:)');
		p_pcc_DA_E(e) == sum(p_pcc_DA_hat(e,:));
			
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for s = 1:nscen
			w_piDA_pRT(e,s) >= sum( pi_pcc_DA_hat(e,:,s)' .* p_pcc_rt_lo_part(:)...
				- pi_pcc_DA_lo .* p_pcc_rt_lo_part(:).*p_pcc_rt_part(e,:,s)' ) + pi_pcc_DA_lo * p_pcc_rt(e,s) ;
			
			w_piDA_pRT(e,s) >= sum( pi_pcc_DA_hat(e,:,s)' .* p_pcc_rt_up_part(:)...
				- pi_pcc_DA_up * p_pcc_rt_up_part(:).*p_pcc_rt_part(e,:,s)' ) +  pi_pcc_DA_up * p_pcc_rt(e,s);
			
			w_piDA_pRT(e,s) <= sum( pi_pcc_DA_hat(e,:,s)' .* p_pcc_rt_lo_part(:)...
				- pi_pcc_DA_up * p_pcc_rt_lo_part(:).*p_pcc_rt_part(e,:,s)' ) + pi_pcc_DA_up * p_pcc_rt(e,s);
			
			w_piDA_pRT(e,s) <= sum( pi_pcc_DA_hat(e,:,s)' .* p_pcc_rt_up_part(:)...
				- pi_pcc_DA_lo * p_pcc_rt_up_part(:) .*p_pcc_rt_part(e,:,s)' ) + pi_pcc_DA_lo * p_pcc_rt(e,s);
		end
		
		w_piDA_pDA(e) >= sum( pi_pcc_DA_lo_part(:) .* p_pcc_DA_hat(e,:)' ...
			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_lo_part(:) * p_pcc_DA_E_lo ) + pi_pcc_DA(e) * p_pcc_DA_E_lo;
		w_piDA_pDA(e) >= sum( pi_pcc_DA_up_part(:) .* p_pcc_DA_hat(e,:)' ...
			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_up_part(:) * p_pcc_DA_E_up ) + pi_pcc_DA(e) * p_pcc_DA_E_up;
		w_piDA_pDA(e) <= sum( pi_pcc_DA_up_part(:) .* p_pcc_DA_hat(e,:)'...
			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_up_part(:) * p_pcc_DA_E_lo ) + pi_pcc_DA(e) * p_pcc_DA_E_lo;
		w_piDA_pDA(e) <= sum(pi_pcc_DA_lo_part(:) .* p_pcc_DA_hat(e,:)' ...
			- pi_pcc_DA_part(e,:)' .* pi_pcc_DA_lo_part(:) * p_pcc_DA_E_up ) + pi_pcc_DA(e) * p_pcc_DA_E_up;
		
		for s = 1:nscen
			w_piup_pup(e,s) >= sum( pi_pcc_up_lo_part(:) .* pup_pcc_hat(e,:,s)'  ...
				- pi_pcc_up_part(e,:,s)' .* pi_pcc_up_lo_part(:) * pup_pcc_lo ) + pi_pcc_up(e) * pup_pcc_lo;
			w_piup_pup(e,s) >= sum( pi_pcc_up_up_part(:) .* pup_pcc_hat(e,:,s)' ...
				- pi_pcc_up_part(e,:,s)' .* pi_pcc_up_up_part(:) * pup_pcc_up ) + pi_pcc_up(e) * pup_pcc_up;
			w_piup_pup(e,s) <= sum( pi_pcc_up_up_part(:) .* pup_pcc_hat(e,:,s)' ...
				-  pi_pcc_up_part(e,:,s)' .* pi_pcc_up_up_part(:) * pup_pcc_lo ) + pi_pcc_up(e) * pup_pcc_lo;
			w_piup_pup(e,s) <= sum( pi_pcc_up_lo_part(:) .* pup_pcc_hat(e,:,s)' ...
				-  pi_pcc_up_part(e,:,s)' .* pi_pcc_up_lo_part(:) * pup_pcc_up ) + pi_pcc_up(e) * pup_pcc_up;
		end
		
		
		for s = 1:nscen
			w_pidn_pdn(e,s) >= pi_pcc_dn_lo * pdn_pcc(e,s) ...
				+ sum(pi_pcc_dn_hat(e,:,s)' .* pdn_pcc_lo_part(:) - pdn_pcc_part(e,:,s)' * pi_pcc_dn_lo .* pdn_pcc_lo_part(:));
			w_pidn_pdn(e,s) >= pi_pcc_dn_up * pdn_pcc(e,s) ...
				+ sum( pi_pcc_dn_hat(e,:,s)' .* pdn_pcc_up_part(:) - pdn_pcc_part(e,:,s)' * pi_pcc_dn_up .* pdn_pcc_up_part(:) );
			w_pidn_pdn(e,s) <= pi_pcc_dn_up * pdn_pcc(e,s) ...
				+ sum( pi_pcc_dn_hat(e,:,s)' .* pdn_pcc_lo_part(:) - pdn_pcc_part(e,:,s)' * pi_pcc_dn_up .* pdn_pcc_lo_part(:));
			w_pidn_pdn(e,s) <= pi_pcc_dn_lo * pdn_pcc(e,s) ...
				+ sum( pi_pcc_dn_hat(e,:,s)' .* pdn_pcc_up_part(:) - pdn_pcc_part(e,:,s)' * pi_pcc_dn_lo .* pdn_pcc_up_part(:));
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		sum(fdn_part(e,:)) == 1;
		sum(rho_DA_plus_part(e,:)) == 1;
		for s = 1:nscen
			sum(rho_RT_minus_part(e,:,s)) == 1;
			sum(fup_part(e,:,s)) == 1;
		end
		
		
		fdn_lo = -4;
		fdn_up = 2;
		
		rhominus_lo = 0;
		rhominus_up = 15;
		
		fup_lo = -2;
		fup_up = 7;
		
		rhoplus_lo = 0;
		rhoplus_up = 15;
		
		fdn_lo_part = zeros(npart2,1);
		fdn_up_part = zeros(npart2,1);
		fup_lo_part = zeros(npart2,1);
		fup_up_part = zeros(npart2,1);
		rhominus_lo_part = zeros(npart2,1);
		rhominus_up_part = zeros(npart2,1);
		rhoplus_lo_part = zeros(npart2,1);
		rhoplus_up_part = zeros(npart2,1);
		
		for pl = 1:npart2
			fdn_lo_part(pl) = fdn_lo - (fdn_lo - fdn_up)*(pl - 1)/npart2;
			fdn_up_part(pl) = fdn_lo + (fdn_up - fdn_lo)*(pl)/npart2;
			
			fup_lo_part(pl) = fup_lo - (fup_lo - fup_up)*(pl - 1)/npart2;
			fup_up_part(pl) = fup_lo + (fup_up - fup_lo)*(pl)/npart2;
			
			rhominus_lo_part(pl) = rhominus_lo - (rhominus_lo - rhominus_up)*(pl - 1)/npart2;
			rhominus_up_part(pl) = rhominus_lo + (rhominus_up - rhominus_lo)*(pl)/npart2;
			
			rhoplus_lo_part(pl) = rhoplus_lo - (rhoplus_lo - rhoplus_up)*(pl - 1)/npart2;
			rhoplus_up_part(pl) = rhoplus_lo + (rhoplus_up - rhoplus_lo)*(pl)/npart2;
        end
		
		rhominus_lo * fdn_part(e,:) <= rho_DA_minus_hat(e,:) <= rhominus_up * fdn_part(e,:);
		sum( fdn_lo_part(:) .* fdn_part(e,:)') <= fe_dn(e) <= sum( fdn_up_part(:) .* fdn_part(e,:)');
		rho_DA_minus(e) == sum(rho_DA_minus_hat(e,:));
		
		fup_lo * rho_DA_plus_part(e,:) <= fup_hat(e,:) <= fup_up * rho_DA_plus_part(e,:);
		sum( rhoplus_lo_part(:) .* rho_DA_plus_part(e,:)') <= rho_DA_plus(e) <= sum( rhoplus_up_part(:) .* rho_DA_plus_part(e,:)');
		fe_up(e) == sum(fup_hat(e,:));
		
		for s = 1:nscen	
			fdn_lo * rho_RT_minus_part(e,:,s) <= fdn_hat(e,:,s) <= fdn_up * rho_RT_minus_part(e,:,s);
			sum( rhominus_lo_part(:) .* rho_RT_minus_part(e,:,s)') <= rho_RT_minus(e,s) <= sum( rhominus_up_part(:) .* rho_RT_minus_part(e,:,s)');
			fe_dn(e) == sum(fdn_hat(e,:,s));
			
			rhoplus_lo * fup_part(e,:,s) <= rho_RT_plus_hat(e,:,s) <= rhoplus_up * fup_part(e,:,s);
			sum( fup_lo_part(:) .* fup_part(e,:,s)') <= fe_up(e) <= sum( fup_up_part(:) .* fup_part(e,:,s)');
			rho_RT_plus(e,s) == sum(rho_RT_plus_hat(e,:,s));
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		w_fdn_rhominus_DA(e) >= sum( fdn_lo_part(:) .* rho_DA_minus_hat(e,:)' ...
			- fdn_part(e,:)' .* fdn_lo_part(:) * rhominus_lo ) + fe_dn(e) * rhominus_lo;
		w_fdn_rhominus_DA(e) >= sum( fdn_up_part(:) .* rho_DA_minus_hat(e,:)' ...
			- fdn_part(e,:)' .* fdn_up_part(:) * rhominus_up ) + fe_dn(e) * rhominus_up ;
		w_fdn_rhominus_DA(e) <= sum( fdn_up_part(:) .* rho_DA_minus_hat(e,:)'  ...
			- fdn_part(e,:)' .* fdn_up_part(:) * rhominus_lo ) + fe_dn(e) * rhominus_lo;
		w_fdn_rhominus_DA(e) <=  fe_dn(e) * rhominus_up + sum( fdn_lo_part(:) .* rho_DA_minus_hat(e,:)' ...
			- fdn_part(e,:)' .* fdn_lo_part(:) * rhominus_up );
		
		w_fup_rhoplus_DA(e) >= fup_lo * rho_DA_plus(e) + sum( fup_hat(e,:)' .* rhoplus_lo_part(:) ...
			- rho_DA_plus_part(e,:)' .* fup_lo .* rhoplus_lo_part(:) );
		w_fup_rhoplus_DA(e) >= fup_up * rho_DA_plus(e) + sum( fup_hat(e,:)' .* rhoplus_up_part(:) ...
			- rho_DA_plus_part(e,:)' .* fup_up .* rhoplus_up_part(:) );
		w_fup_rhoplus_DA(e) <= fup_up * rho_DA_plus(e) + sum( fup_hat(e,:)' .* rhoplus_lo_part(:) ...
			- rho_DA_plus_part(e,:)' .* fup_up .* rhoplus_lo_part(:) );
		w_fup_rhoplus_DA(e) <= sum( fup_hat(e,:)' .* rhoplus_up_part(:) ...
			- rho_DA_plus_part(e,:)' .* fup_lo .* rhoplus_up_part(:) ) + fup_lo * rho_DA_plus(e);
		
		for s = 1:nscen
			w_fdn_rhominus_RT(e,s) >= fdn_lo * rho_RT_minus(e,s) + fdn_hat(e,:,s) * rhominus_lo_part(:)...
				- rho_RT_minus_part(e,:,s) * fdn_lo * rhominus_lo_part(:);
			w_fdn_rhominus_RT(e,s) >= fdn_up * rho_RT_minus(e,s) + fdn_hat(e,:,s) * rhominus_up_part(:)...
				- rho_RT_minus_part(e,:,s) * fdn_up * rhominus_up_part(:);
			w_fdn_rhominus_RT(e,s) <= fdn_up * rho_RT_minus(e,s) + fdn_hat(e,:,s) * rhominus_lo_part(:)...
				- rho_RT_minus_part(e,:,s) * fdn_up * rhominus_lo_part(:);
			w_fdn_rhominus_RT(e,s) <= fdn_hat(e,:,s) * rhominus_up_part(:) + fdn_lo * rho_RT_minus(e,s)...
				- rho_RT_minus_part(e,:,s) * fdn_lo * rhominus_up_part(:);
		
			w_fup_rhoplus_RT(e,s) >= rho_RT_plus_hat(e,:,s) * fup_lo_part(:) + fe_up(e) * rhoplus_lo ...
				- fup_part(e,:,s) * fup_lo_part(:) * rhoplus_lo;
			w_fup_rhoplus_RT(e,s) >=  rho_RT_plus_hat(e,:,s) * fup_up_part(:) + fe_up(e) * rhoplus_up ...
				- fup_part(e,:,s) * fup_up_part(:) * rhoplus_up;
			w_fup_rhoplus_RT(e,s) <= rho_RT_plus_hat(e,:,s) * fup_up_part(:) + fe_up(e) * rhoplus_lo...
				- fup_part(e,:,s) * fup_up_part(:) * rhoplus_lo;
			w_fup_rhoplus_RT(e,s) <= fe_up(e) * rhoplus_up +  rho_RT_plus_hat(e,:,s) * fup_lo_part(:) ...
				- fup_part(e,:,s) * fup_lo_part(:) * rhoplus_up;
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
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

		%% real time primal constraints



% 		branch_num = branch_num_area{e+1};

% 		fbus = fbus_local{e+1};
% 		tbus = tbus_local{e+1};

		A = [2	0	0	0;
			0	2	0	0;
			0	0	1	-1];
		b = [0	0	1	1];

		for k = 1:nl_d
			for s = 1:nscen
				test_vec{s}(:,k) = [p_flow(fbus(k),tbus(k),s); q_flow(fbus(k),tbus(k),s); ...
					i_sq(fbus(k),tbus(k),s); v_sq(fbus(k),s)];
				norm(A*test_vec{s}(:,k)) <= b*test_vec{s}(:,k);

				test_vec2{s}(:,k) = [p_flow(tbus(k),fbus(k),s); q_flow(tbus(k),fbus(k),s);...
					i_sq(tbus(k),fbus(k),s); v_sq(tbus(k),s)];
				norm(A*test_vec2{s}(:,k)) <= b*test_vec2{s}(:,k);
			end
		end

		for s = 1:nscen
			for k = 1:length(fbus)
				p_flow(fbus(k),tbus(k),s) + p_flow(tbus(k),fbus(k),s) == rBranch(branch_num(k)).*i_sq(fbus(k),tbus(k),s);
				p_flow(tbus(k),fbus(k),s) + p_flow(fbus(k),tbus(k),s) == rBranch(branch_num(k)).*i_sq(tbus(k),fbus(k),s);

				q_flow(fbus(k),tbus(k),s) + q_flow(tbus(k),fbus(k),s) == xBranch(branch_num(k)).*i_sq(fbus(k),tbus(k),s);
				q_flow(tbus(k),fbus(k),s) + q_flow(fbus(k),tbus(k),s) == xBranch(branch_num(k)).*i_sq(tbus(k),fbus(k),s);

				v_sq(tbus(k),s) == v_sq(fbus(k),s) - 2*(rBranch(branch_num(k)).*p_flow(fbus(k),tbus(k),s)...
					+ xBranch(branch_num(k)).*q_flow(fbus(k),tbus(k),s) ) ...
					+ (rBranch(branch_num(k))^2+xBranch(branch_num(k))^2).*i_sq(fbus(k),tbus(k),s);
			end
		end

		for s = 1:nscen
			Vmin(buses).^2 <= v_sq(buses,s) <= Vmax(buses).^2;

			sum(q_flow(buses,buses,s),2) == q_inj_rt(buses,s);

			q_inj_rt(buses,s) == incidentG(gens_area,buses)' * q_gen_RT(gens_area,s) ...
						- incidentD(dems_area,buses)' * q_dem_RT(dems_area,s) + shed_q_RT(buses,s) + incident_pcc * q_pcc_rt(e,s);


			for k = 1:nl_d
				if SlmMax(branch_num(k)) ~= 0
					p_flow(fbus(k),tbus(k),s).^2 + q_flow(fbus(k),tbus(k),s).^2 <= SlmMax(branch_num(k)).^2;
					p_flow(tbus(k),fbus(k),s).^2 + q_flow(tbus(k),fbus(k),s).^2 <= SlmMax(branch_num(k)).^2;
				end
			end
        end
        
		p_flow_matrix = reshape(sum(p_flow(buses,buses,:),2),size(p_inj_rt(buses,:)));
		p_flow_matrix ==  p_inj_rt(buses,:);
		for s = 1:nscen
			p_DR_RT(DR_area,s) == p_DR_DA(DR_area) - pup_DR(DR_area,s) + pdn_DR(DR_area,s);
			p_pcc_rt(e,s) == p_pcc_DA_E(e) + pup_pcc(e,s) - pdn_pcc(e,s);
			p_gen_RT(gens_area,s) == p_gen_DA_tilde(gens_area) + pup_g(gens_area,s) - pdn_g(gens_area,s);
			p_dem_RT(dems_area,s) == p_dem_DA_tilde(dems_area) - pup_d(dems_area,s) + pdn_d(dems_area,s);

			p_inj_rt(buses,s) == incidentG(gens_area,buses)' * p_gen_RT(gens_area,s) - incidentD(dems_area,buses)' * p_dem_RT(dems_area,s)...
					+ incidentW(windg_area,buses)' * wind_RT(windg_area,s) + shed_p_RT(buses,s) + incident_pcc * p_pcc_rt(e,s)...
					- incidentDR(DR_area,buses)' * p_DR_RT(DR_area,s);

			0 <= shed_p_RT(buses,s) <= incidentD(dems_area,buses)' * p_dem_RT(dems_area,s);
		end


		for s = 1:nscen
			Pmin_gen(gens_area) <= p_gen_RT(gens_area,s) <= Pmax_gen(gens_area);
			Qmin_gen(gens_area) <= q_gen_RT(gens_area,s) <= Qmax_gen(gens_area);

			Pmin_dem(dems_area) <= p_dem_RT(dems_area,s) <= Pmax_dem(dems_area);
			Qmin_dem(dems_area) <= q_dem_RT(dems_area,s) <= Qmax_dem(dems_area);

			Wmin(windg_area,s) <= wind_RT(windg_area,s) <= Wmax(windg_area,s);
			
			Pmin_DR(DR_area) <= p_DR_RT(DR_area,s) <= Pmax_DR(DR_area);
		end

		fe_dn(e) <= p_pcc_rt(e,:) <= fe_up(e);
	% 	0 <= shed_p_RT <= incidentD' * p_dem_RT;
	
% 	for s = 1:nscen
% 		pup_pcc(e,s) <= VOLL * p_pcc_comp(e,s);
% 		pdn_pcc(e,s) <= VOLL * (1 - p_pcc_comp(e,s));
% 		
% 		pup_g(gens_area,s) <= VOLL * p_gen_comp(gens_area,s);
% 		pdn_g(gens_area,s) <= VOLL * (1 - p_gen_comp(gens_area,s));
% 		
% 		pup_d(dems_area,s) <= VOLL * p_dem_comp(dems_area,s);
% 		pdn_d(dems_area,s) <= VOLL * (1 - p_dem_comp(dems_area,s));
% 		
% 		pup_DR(DR_area,s) <= VOLL * p_DR_comp(DR_area,s);
% 		pdn_DR(DR_area,s) <= VOLL * (1 - p_DR_comp(DR_area,s));
% 	end
		%% KKT stationary conditions
		for g = 1:ng_d
			offer_gen_DA(gens_area(g),2) - lambda_DA_e(e) - varsigma_gen_DA_minus_e(gens_area(g)) + varsigma_gen_DA_plus_e(gens_area(g)) ...
					- sum(prob_wscen(1,:)*(offer_gen_DA(gens_area(g),2))) + sum(zeta_p_gen_RT(gens_area(g),:))  == 0; % p_gen_DA
			for s = 1:nscen
				prob_wscen(1,s) * offer_gen_upreg(gens_area(g)) + zeta_p_gen_RT(gens_area(g),s) - epsilon_up_RT(gens_area(g),s) == 0; % pup_g
				prob_wscen(1,s) * offer_gen_downreg(gens_area(g)) - zeta_p_gen_RT(gens_area(g),s) - epsilon_dn_RT(gens_area(g),s) == 0; % pdn_g

				prob_wscen(1,s)*offer_gen_DA(gens_area(g),2) - zeta_p_gen_RT(gens_area(g),s) - varsigma_gen_RT_minus(gens_area(g),s) ...
					+ varsigma_gen_RT_plus(gens_area(g),s) - incidentG(gens_area(g),buses) *lambda_p_RT(buses,s) == 0; % p_gen_RT

				- kappa_gen_RT_minus(gens_area(g),s) ...
					+ kappa_gen_RT_plus(gens_area(g),s) - incidentG(gens_area(g),buses) *lambda_q_RT(buses,s) == 0; % q_gen_RT
			end
		end

		for d = 1:nd_d
			-bid_dem_DA(dems_area(d)) + lambda_DA_e(e) - varsigma_dem_DA_minus_e(dems_area(d)) + varsigma_dem_DA_plus_e(dems_area(d)) ...
					+  sum(prob_wscen(1,:)*bid_dem_DA(dems_area(d)) ) + sum(zeta_p_dem_RT(dems_area(d),:)) - Upsilon_DA_plus(e) == 0; % p_dem_DA
			for s = 1:nscen
				prob_wscen(1,s) * bid_dem_upreg(dems_area(d)) - zeta_p_dem_RT(dems_area(d),s) - varepsilon_up_RT(dems_area(d),s) == 0; % pup_d
				prob_wscen(1,s) * bid_dem_downreg(dems_area(d)) + zeta_p_dem_RT(dems_area(d),s) - varepsilon_dn_RT(dems_area(d),s) == 0; % pdn_d

				-prob_wscen(1,s)*bid_dem_DA(dems_area(d)) - zeta_p_dem_RT(dems_area(d),s) - varsigma_dem_RT_minus(dems_area(d),s) ...
					+ varsigma_dem_RT_plus(dems_area(d),s) + incidentD(dems_area(d),buses) * (lambda_p_RT(buses,s) - Upsilon_RT_plus(buses,s)) == 0; % p_dem_RT

				- kappa_dem_RT_minus(dems_area(d),s) + kappa_dem_RT_plus(dems_area(d),s) ...
					+ incidentD(dems_area(d),buses) *lambda_q_RT(buses,s) - incidentD(dems_area(d),buses)*Upsilon_Q_RT_plus(buses,s) == 0; % q_dem_RT
			end
		end
		
		for d = 1:nDR_d
			-bid_DA_DR(DR_area(d)) + lambda_DA_e(e) - varsigma_DR_DA_minus_e(DR_area(d)) + varsigma_DR_DA_plus_e(DR_area(d)) ...
				+ sum(prob_wscen(1,:)*bid_DA_DR(DR_area(d)) ) + sum(zeta_DR_RT(DR_area(d),:)) == 0; % p_DR_DA
			for s = 1:nscen
				prob_wscen(1,s) * bid_up_DR(DR_area(d)) - zeta_DR_RT(DR_area(d),s) - epsilon_DR_up_RT(DR_area(d),s) == 0; % pup_DR
				prob_wscen(1,s) * bid_dn_DR(DR_area(d)) + zeta_DR_RT(DR_area(d),s) - epsilon_DR_dn_RT(DR_area(d),s) == 0; % pdn_DR

				-prob_wscen(1,s)*bid_DA_DR(DR_area(d)) - zeta_DR_RT(DR_area(d),s) - varsigma_DR_RT_minus(DR_area(d),s) ...
					+ varsigma_DR_RT_plus(DR_area(d),s) + incidentDR(DR_area(d),buses) * lambda_p_RT(buses,s) == 0; % p_DR_RT
			end
		end

		for n = 1:nb_d
			for s = 1:nscen
				VOLL - lambda_p_RT(buses(n),s) - Upsilon_RT_minus(buses(n),s) + Upsilon_RT_plus(buses(n),s) == 0; % shed_p_RT
				VOLL - lambda_q_RT(buses(n),s) - Upsilon_Q_RT_minus(buses(n),s) + Upsilon_Q_RT_plus(buses(n),s) == 0;
			end
		end

		for w = 1:nw_d
			for s = 1:nscen
				offer_wind - incidentW(windg_area(w),buses) * lambda_p_RT(buses,s) - nu_RT_minus(windg_area(w),s) + nu_RT_plus(windg_area(w),s) == 0; % wind_RT
			end
		end


		for l = 1:nl_d
			for s = 1:nscen

				lambda_p_RT(fbus(l),s) + eta_RT_p(fbus(l),tbus(l),s) + 2*gamma_RT_p(fbus(l),tbus(l),s) ...
					- mu_p_RT(fbus(l),tbus(l),s) - mu_p_RT(tbus(l),fbus(l),s) ...
					- 2*beta_RT(fbus(l),tbus(l),s)*rBranch(branch_num(l)) == 0; % p_flow -->


				lambda_p_RT(tbus(l),s) + eta_RT_p(tbus(l),fbus(l),s) + 2*gamma_RT_p(tbus(l),fbus(l),s) ...
					- mu_p_RT(tbus(l),fbus(l),s) - mu_p_RT(fbus(l),tbus(l),s) ...
					- 2*beta_RT(tbus(l),fbus(l),s)*rBranch(branch_num(l)) == 0; % p_flow <--


				lambda_q_RT(fbus(l),s) + eta_RT_q(fbus(l),tbus(l),s) + 2*gamma_RT_q(fbus(l),tbus(l),s) ...
					- mu_q_RT(fbus(l),tbus(l),s) - mu_q_RT(tbus(l),fbus(l),s) ...
					- 2*beta_RT(fbus(l),tbus(l),s)*xBranch(branch_num(l)) == 0; % q_flow -->
				lambda_q_RT(tbus(l),s) + eta_RT_q(tbus(l),fbus(l),s) + 2*gamma_RT_q(tbus(l),fbus(l),s) ...
					- mu_q_RT(tbus(l),fbus(l),s) - mu_q_RT(fbus(l),tbus(l),s) ...
					- 2*beta_RT(tbus(l),fbus(l),s)*xBranch(branch_num(l)) == 0; % q_flow <--


				gamma_RT_l(fbus(l),tbus(l),s) - gamma_RT_w(fbus(l),tbus(l),s) + mu_p_RT(fbus(l),tbus(l),s)...
					*rBranch(branch_num(l)) + mu_q_RT(fbus(l),tbus(l),s)*xBranch(branch_num(l)) ...
					+ beta_RT(fbus(l),tbus(l),s)*(rBranch(branch_num(l))^2 + xBranch(branch_num(l))^2) == 0; % i_sq -->
				gamma_RT_l(tbus(l),fbus(l),s) - gamma_RT_w(tbus(l),fbus(l),s) + mu_p_RT(tbus(l),fbus(l),s)...
					*rBranch(branch_num(l)) + mu_q_RT(tbus(l),fbus(l),s)*xBranch(branch_num(l)) ...
					+ beta_RT(tbus(l),fbus(l),s)*(rBranch(branch_num(l))^2 + xBranch(branch_num(l))^2) == 0; % i_sq -->



				-sigma_RT_minus(fbus(l),s) + sigma_RT_plus(fbus(l),s)...
					- sum(gamma_RT_l(fbus(l),:,s) + gamma_RT_w(fbus(l),:,s)...
					- beta_RT(fbus(l),:,s) + beta_RT(:,fbus(l),s)' )== 0; % v_sq
				-sigma_RT_minus(tbus(l),s) + sigma_RT_plus(tbus(l),s)...
					- sum(gamma_RT_l(tbus(l),:,s) + gamma_RT_w(tbus(l),:,s)...
					- beta_RT(tbus(l),:,s) + beta_RT(:,tbus(l),s)' )== 0; % v_sq

				norm([eta_RT_p(fbus(l),tbus(l),s); eta_RT_q(fbus(l),tbus(l),s)]) <= eta_RT_s(fbus(l),tbus(l),s); % SOCP constraint for eta
				norm([eta_RT_p(tbus(l),fbus(l),s); eta_RT_q(tbus(l),fbus(l),s)]) <= eta_RT_s(tbus(l),fbus(l),s); % SOCP constraint for eta

				norm([gamma_RT_p(fbus(l),tbus(l),s); gamma_RT_q(fbus(l),tbus(l),s); gamma_RT_l(fbus(l),tbus(l),s)])...
					<= gamma_RT_w(fbus(l),tbus(l),s); % SOCP constraint for gamma
				norm([gamma_RT_p(tbus(l),fbus(l),s); gamma_RT_q(tbus(l),fbus(l),s); gamma_RT_l(tbus(l),fbus(l),s)])...
					<= gamma_RT_w(tbus(l),fbus(l),s); % SOCP constraint for gamma

			end
		end

		for s = 1:nscen

			prob_wscen(1,s)*pi_pcc_DA(e) - incident_pcc' * lambda_p_RT(buses,s) - zeta_pcc_RT(e,s) - rho_RT_minus(e,s)...
				+ rho_RT_plus(e,s) == 0; % p_pcc_RT

			incident_pcc' * lambda_q_RT(buses,s) == 0; % q_pcc_RT

			prob_wscen(1,s)*pi_pcc_up(e) + zeta_pcc_RT(e,s) - epsilon_pcc_up_RT(e,s) == 0; % pup_pcc
			prob_wscen(1,s)*pi_pcc_dn(e) - zeta_pcc_RT(e,s) - epsilon_pcc_dn_RT(e,s) == 0; % pdn_pcc
		end
		pi_pcc_DA(e) - sum(prob_wscen(1,:)*pi_pcc_DA(e)) - lambda_DA_e(e) - rho_DA_minus(e)...
			+ rho_DA_plus(e) + sum(zeta_pcc_RT(e,:)) == 0; % p_pcc_DA

		VOLL - lambda_DA_e(e) - Upsilon_DA_minus(e) + Upsilon_DA_plus(e) == 0; % shed_p_DA
		offer_wind - lambda_DA_e(e) - iota_minus_DA(e) + iota_plus_DA(e) == 0; % wind_DA
		
		
		
        
        %%
		rho_DA_minus(e) <= VOLL* rho_DA_comp(e);
		rho_DA_plus(e) <= VOLL * (1 - rho_DA_comp(e));
		
		rho_RT_minus(e,:) <= VOLL* rho_RT_comp(e,:);
		rho_RT_plus(e,:) <= VOLL * (1 - rho_RT_comp(e,:));
		


		%% dual constraints
			
% 		for s = 1:nscen
% 			for l = 1:nl_d
% 				
% 				ui{s}(:,l) = A'*[gamma_RT_p(fbus(l),tbus(l),s); gamma_RT_q(fbus(l),tbus(l),s); gamma_RT_l(fbus(l),tbus(l),s)];
% 				bi{s}(:,l) = b' * gamma_RT_w(fbus(l),tbus(l),s);
% 				
% % 				sum( ui{s}(:,l) - bi{s}(:,l) ) == 0;
% 			end
% 		end

		%% strong duality objective

		dual_objective_DA(e) = sum(Pmin_gen(gens_area)' * varsigma_gen_DA_minus_e(gens_area)...
			- Pmax_gen(gens_area)' * varsigma_gen_DA_plus_e(gens_area) ) ...
			+ sum(Pmin_dem(dems_area)' * varsigma_dem_DA_minus_e(dems_area) ...
			- Pmax_dem(dems_area)' * varsigma_dem_DA_plus_e(dems_area) )...
			+ sum(Pmin_DR(DR_area)' * varsigma_DR_DA_minus_e(DR_area) ...
			- Pmax_DR(DR_area)' * varsigma_DR_DA_plus_e(DR_area) )...
			- iota_plus_DA(e) * sum(Wmax_mean_DA(windg_area)) ...
			+ w_fdn_rhominus_DA(e) - w_fup_rhoplus_DA(e);
% 			+ fe_dn(e) * rho_DA_minus(e) - fe_up(e) * rho_DA_plus(e);

		for s = 1:nscen
			dual_objective_RT(e,s) =  sum(Vmin(buses)'.^2 * sigma_RT_minus(buses,s))  - sum(Vmax(buses)'.^2 * sigma_RT_plus(buses,s)) ...
				- sum(nu_RT_plus(windg_area,s)' *  Wmax(windg_area,s) )...
				+ sum(Pmin_gen(gens_area)' * varsigma_gen_RT_minus(gens_area,s) - Pmax_gen(gens_area)' * varsigma_gen_RT_plus(gens_area,s)) ...
				+ sum(Pmin_dem(dems_area)' * varsigma_dem_RT_minus(dems_area,s) - Pmax_dem(dems_area)' * varsigma_dem_RT_plus(dems_area,s)) ...
				+ sum(Pmin_DR(DR_area)' * varsigma_DR_RT_minus(DR_area,s) - Pmax_DR(DR_area)' * varsigma_DR_RT_plus(DR_area,s)) ...
				+ sum(Qmin_gen(gens_area)' * kappa_gen_RT_minus(gens_area,s) - Qmax_gen(gens_area)' * kappa_gen_RT_plus(gens_area,s)) ...
				+ sum(Qmin_dem(dems_area)' * kappa_dem_RT_minus(dems_area,s) - Qmax_dem(dems_area)' * kappa_dem_RT_plus(dems_area,s)) ...
				+ w_fdn_rhominus_RT(e,s) - w_fup_rhoplus_RT(e,s);
% 				+ fe_dn(e) * rho_RT_minus(e,s) - fe_up(e) *
% 				rho_RT_plus(e,s);
				
			for l = 1:nl_d
				if SlmMax(branch_num(l)) ~= 0
					dual_obj_linecap(e,s) = dual_obj_linecap(e,s) - sum(  (eta_RT_s(fbus(l),tbus(l),s)+eta_RT_s(tbus(l),fbus(l),s)) * SlmMax(branch_num(l)) ) ;
				else
					eta_RT_s(fbus(l),tbus(l),s) + eta_RT_s(tbus(l),fbus(l),s) == 0;
				end
			end
		end
		dual_objective(e) = dual_objective_DA(e) + sum(dual_objective_RT(e,:)) + sum(dual_obj_linecap(e,:));

		DSO_cost(e) == dual_objective(e);
		
		sigma_RT_minus(buses,:) <= VOLL* sigma_RT_comp(buses,:);
		sigma_RT_plus(buses,:) <= VOLL * (1 - sigma_RT_comp(buses,:));
% 		Upsilon_RT_minus(buses,:) <= VOLL* Upsilon_RT_comp(buses,:);
% 		Upsilon_RT_plus(buses,:) <= VOLL * (1 - Upsilon_RT_comp(buses,:));
		
		
	end

	varsigma_gen_DA_minus_e <= VOLL* varsigma_gen_DA_comp;
	varsigma_gen_DA_plus_e <= VOLL * (1 - varsigma_gen_DA_comp);

	varsigma_dem_DA_minus_e <= VOLL* varsigma_dem_DA_comp;
	varsigma_dem_DA_plus_e <= VOLL * (1 - varsigma_dem_DA_comp);
	
	varsigma_gen_RT_minus <= VOLL* varsigma_gen_RT_comp;
	varsigma_gen_RT_plus <= VOLL * (1 - varsigma_gen_RT_comp);

	varsigma_dem_RT_minus <= VOLL* varsigma_dem_RT_comp;
	varsigma_dem_RT_plus <= VOLL * (1 - varsigma_dem_RT_comp);
	
	varsigma_DR_RT_minus <= VOLL* varsigma_DR_RT_comp;
	varsigma_DR_RT_plus <= VOLL * (1 - varsigma_DR_RT_comp);
	
% 	Upsilon_DA_minus <= VOLL* Upsilon_DA_comp;
% 	Upsilon_DA_plus <= VOLL * (1 - Upsilon_DA_comp);
	

cvx_end


for s = 1:nscen
	for k = 1:size(test_vec{s},2)
		tight1(k,s) = ( b*test_vec{s}(:,k) - norm(A*test_vec{s}(:,k))  );
		tight2(k,s) = ( b*test_vec2{s}(:,k) - norm(A*test_vec2{s}(:,k)) );
		if any( tight1(k,s) > 10e-6) ...
				|| any( tight2(k,s) > 10e-6)
			warning(['SOCP relaxation of scenario ' num2str(s) ' and line ' num2str(k) ' is not tight']);
			pause(0.05)
		end
	end
end
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
DA_market_outcome.p_gen_DA_tilde = p_gen_DA_tilde;
DA_market_outcome.p_dem_DA_tilde = p_dem_DA_tilde;


if DA_market_outcome.shed_p > 10e-5
	warning('There is load shedding in the Day-Ahead market');
	pause(3)
end

result_DSO_lookahead.DSO_cost = DSO_cost;
result_DSO_lookahead.cost_DR_DA = cost_DR_DA;
result_DSO_lookahead.p_gen_DA_tilde = p_gen_DA_tilde;
result_DSO_lookahead.p_dem_DA_tilde = p_dem_DA_tilde;
result_DSO_lookahead.p_DR_DA = p_DR_DA;
result_DSO_lookahead.p_DR_RT = p_DR_RT;
result_DSO_lookahead.wind_RT = wind_RT;
result_DSO_lookahead.shed_p_RT = shed_p_RT;
result_DSO_lookahead.shed_q_RT = shed_q_RT;
result_DSO_lookahead.wind_DA_E = wind_DA_E;
result_DSO_lookahead.p_pcc_DA_E = p_pcc_DA_E;
result_DSO_lookahead.p_pcc_RT = p_pcc_rt;
result_DSO_lookahead.p_flow = p_flow;
result_DSO_lookahead.q_flow = q_flow;
result_DSO_lookahead.SOCP_tight1 = tight1;
result_DSO_lookahead.SOCP_tight2 = tight2;


result_DSO_lookahead.fe_up = fe_up;
result_DSO_lookahead.fe_dn = fe_dn;
result_DSO_lookahead.pi_pcc_DA = pi_pcc_DA;
result_DSO_lookahead.pi_pcc_up = pi_pcc_up;
result_DSO_lookahead.pi_pcc_dn = pi_pcc_dn;
end