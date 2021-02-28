function [result_lookahead] = DSO_lookahead(data,fe_up,fe_dn,...
	pi_pcc_DA,pi_pcc_up,pi_pcc_dn,e)



VOLL = 10000;
offer_wind = 0;

cc{1} = double.empty(1,0);
[nb,busL,busN,busType,Pd,Qd,Vmax,Vmin,statusG,activeG,activeD,...
	ng,nd,genL,genN,incidentG,demL,incidentD,incidentW,Qmax_gen,Qmin_gen,...
    Pmax_gen,Pmin_gen,Qmax_dem,Qmin_dem,Pmax_dem,Pmin_dem,statusL,...
	activeL,nl,fbusL,tbusL,SlmMax,fbusN,tbusN,incidentF,incidentT,...
    Yf,Yt,YfP,YtP,Ybus,edges,offer_gen_DA,offer_gen_upreg,...
	offer_gen_downreg,bid_dem_DA,bid_dem_upreg,bid_dem_downreg,...
	nr_vals_cost,p,cost_type,ys,areas,bus_areacodes,....
	incident_singlephase,Y_ft,rBranch,xBranch,narea,prob_wscen,Wmax,Wmin,...
	n_wgen,nscen,Wmax_mean_DA,bus_wgen] = Data_Reader(data,cc);

if ~isempty(data.DR_DSO)
	bus_DR = data.DR_DSO(:,1);
	Pmax_DR = data.DR_DSO(:,2);
	Pmin_DR = data.DR_DSO(:,3);
	bid_DA_DR = data.DR_DSO(:,4);
	bid_up_DR = data.DR_DSO(:,5);
	bid_dn_DR = data.DR_DSO(:,6);

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


[area_codes, n_areas, DSO_codes, fbusL_area, tbusL_area, branch_num_area, ...
	fbusL_ext_area, tbusL_ext_area, branch_num_ext_area, overlap,...
	overlap_t, neigh_area, ext_area, ext_area_singlephase,fbus_local,tbus_local,incidence_area] ...
	= find_overlap(bus_areacodes, fbusL, tbusL, areas);

pcc_global = intersect(overlap{1,e+1},areas{e+1});
pcc_local = areas{e+1} == pcc_global;

ne = n_areas-1;
[~,~,gens_area] = intersect(areas{e+1},genL);
[~,~,dems_area] = intersect(areas{e+1},demL);
[~,~,windg_area] = intersect(areas{e+1},bus_wgen);
[~,~,DR_area] = intersect(areas{e+1},bus_DR);

buses = areas{e+1};
nb = length(areas{e+1});
ng = length(gens_area);
nd = length(dems_area);
nw = length(windg_area);
nDR_d = length(DR_area);
nDR = length(DR_area);

incidentG = incidentG(gens_area,buses);
incidentD = incidentD(dems_area,buses);
incidentW = full(incidentW(windg_area,buses));
incidentDR = incidentDR(DR_area,buses);
incident_pcc = zeros(nb,1);
incident_pcc(pcc_local) = 1;

cvx_begin
	
% 	cvx_precision high
	
	
	variable p_gen_DA(ng)
	variable p_dem_DA(nd) nonnegative
	variable wind_DA nonnegative
	variable shed_p_DA nonnegative
	variable p_pcc_DA
	variable p_DR_DA(nDR)

	variable p_DR_RT(nDR,nscen)
	variable pup_DR(nDR,nscen) nonnegative
	variable pdn_DR(nDR,nscen) nonnegative
	variable pup_g(ng,nscen) nonnegative
	variable pdn_g(ng,nscen) nonnegative
	variable pup_d(nd,nscen) nonnegative
	variable pdn_d(nd,nscen) nonnegative
	variable p_pcc_rt(nscen)
	variable q_pcc_rt(nscen)
	variable pup_pcc(nscen) nonnegative
	variable pdn_pcc(nscen) nonnegative
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

	expression cost_g_RT(ng,nscen)
	expression cost_d_RT(nd,nscen)
	expression cost_RT
	expression cost_DR_RT(nDR,nscen)
	expression cost_DR_DA(nDR)
	
	dual variable lambda_RT
	dual variable lambda_DA
	dual variable lambda_q_RT{nscen}
	
	
	dual variable varsigma_gen_DA_minus
	dual variable varsigma_gen_DA_plus
	dual variable varsigma_dem_DA_minus
	dual variable varsigma_dem_DA_plus
	
	dual variable varsigma_gen_RT_minus{nscen}
	dual variable varsigma_gen_RT_plus{nscen}
	dual variable varsigma_dem_RT_minus{nscen}
	dual variable varsigma_dem_RT_plus{nscen}
	dual variable Upsilon_RT_minus{nscen}
	dual variable Upsilon_RT_plus{nscen}
	dual variable zeta_gen_RT{nscen}
	dual variable zeta_dem_RT{nscen}
	dual variable zeta_pcc_RT{nscen}
	dual variable Upsilon_Q_RT_minus{nscen}
	
	dual variable mu_p_RT{nb,nb,nscen}
	dual variable mu_q_RT{nb,nb,nscen}
	
	dual variable beta_RT{nb,nb,nscen}

	dual variable gamma_RT{nb,nb,nscen}
	
	dual variable eta_RT{nb,nb,nscen}
	
	dual variable nu_RT_plus{nscen}
	dual variable nu_RT_minus{nscen}
	
	dual variable iota_plus_DA
	dual variable iota_minus_DA
	dual variable rho_DA_minus
	dual variable rho_DA_plus
	dual variable sigma_RT_minus{nscen}
	dual variable sigma_RT_plus{nscen}
	dual variable kappa_dem_RT_plus{nscen}
	dual variable kappa_dem_RT_minus{nscen}
	dual variable kappa_gen_RT_plus{nscen}
	dual variable kappa_gen_RT_minus{nscen}
	dual variable rho_RT_minus
	dual variable rho_RT_plus
	
	
	%% Real time objective function
	for m = 1:ng
		cost_g_RT(m,:) = offer_gen_DA(gens_area(m),2).*(p_gen_RT(m,:) - p_gen_DA(m)) ...
			+ offer_gen_upreg(gens_area(m))*pup_g(m,:)...
			+ offer_gen_downreg(gens_area(m)) *pdn_g(m,:);
	end

	for m = 1:nd
		cost_d_RT(m,:) = bid_dem_DA(dems_area(m)).*(p_dem_DA(m) - p_dem_RT(m,:)) ...
			+ bid_dem_upreg(dems_area(m)).*pup_d(m,:) ...
			+ bid_dem_downreg(dems_area(m)) .* pdn_d(m,:);
	end
	cost_shed_RT = VOLL*(shed_p_RT+shed_q_RT);

	cost_wind_RT = wind_RT*offer_wind;
	
	cost_PCC_RT = pi_pcc_DA * (p_pcc_rt - p_pcc_DA) ...
			+ pi_pcc_up * pup_pcc...
			+ pi_pcc_dn * pdn_pcc;
	for m = 1:nDR_d
		cost_DR_RT((m),:) = bid_DA_DR(DR_area(m)).*(p_DR_DA((m)) - p_DR_RT((m),:)) ...
			+ bid_up_DR(DR_area(m)).*pup_DR((m),:) ...
			+ bid_dn_DR(DR_area(m)) .* pdn_DR((m),:);
	end
	
	cost_RT = sum(cost_g_RT,1) + sum(cost_d_RT,1) + sum(cost_shed_RT,1) + sum(cost_wind_RT,1)+cost_PCC_RT' + sum(cost_DR_RT,1);
	
	%% day ahead cost 
	expression cost_gen_DA(ng)
	expression cost_dem_DA(nd)
	for m = 1:ng
		if cost_type(m) == 2
			for k = 1:nr_vals_cost
				cost_gen_DA(m) = cost_gen_DA(m) + sum(offer_gen_DA(gens_area(m),k).*p_gen_DA(m).^(nr_vals_cost-k));
			end
		elseif cost_type(m) == 1
			error('Wrong cost function, piecewise linear is not convex');
		end
	end

	for m = 1:nd
		cost_dem_DA(m) = bid_dem_DA(dems_area(m))*p_dem_DA(m);
	end

	for m = 1:nDR_d
		cost_DR_DA((m)) = bid_DA_DR(DR_area(m))*p_DR_DA(m);
	end

	cost_shed_DA = VOLL*shed_p_DA;
	cost_wind_DA = wind_DA*offer_wind;
	cost_PCC_DA = pi_pcc_DA * p_pcc_DA;
	cost_DA = sum(cost_gen_DA)-sum(cost_dem_DA) + cost_shed_DA + cost_wind_DA + cost_PCC_DA - sum(cost_DR_DA);
	cost = prob_wscen(1,:)*cost_RT' + cost_DA;
%% cost function statement
	minimize(prob_wscen(1,:)*cost_RT' + cost_DA)


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% primal constraints	
	subject to
	%% day-ahead primal constraints
	lambda_DA : sum(p_gen_DA) + wind_DA + shed_p_DA + p_pcc_DA == sum(p_dem_DA) + sum(p_DR_DA);


	varsigma_gen_DA_minus : Pmin_gen(gens_area) <= p_gen_DA <= Pmax_gen(gens_area) : varsigma_gen_DA_plus;
	varsigma_dem_DA_minus : Pmin_dem(dems_area) <= p_dem_DA <= Pmax_dem(dems_area) : varsigma_dem_DA_plus;
	iota_minus_DA : 0 <= wind_DA <= sum(Wmax_mean_DA(windg_area)) : iota_plus_DA;
	0 <= shed_p_DA <= sum(p_dem_DA);
	
	rho_DA_minus : fe_dn <= p_pcc_DA <= fe_up : rho_DA_plus;
	Pmin_DR(DR_area) <= p_DR_DA <= Pmax_DR(DR_area);

	%% real time primal constraints
	
	incident_contra = ones(nb,nb) - incidence_area{e+1};
	incident_contra = logical(incident_contra);

	[n1,n2] = ind2sub([nb nb],find(incident_contra));
	for k = 1:length(n1)
		p_flow(n1(k),n2(k),:) == 0;
		q_flow(n1(k),n2(k),:) == 0;
	end

	branch_num = branch_num_area{e+1};

	fbus = fbus_local{e+1};
	tbus = tbus_local{e+1};

	A = [2	0	0	0;
		0	2	0	0;
		0	0	1	-1];
	b = [0	0	1	1];

	for k = 1:length(fbus)
		for s = 1:nscen
			test_vec{s}(:,k) = [p_flow(fbus(k),tbus(k),s); q_flow(fbus(k),tbus(k),s); ...
				i_sq(fbus(k),tbus(k),s); v_sq(fbus(k),s)];
			gamma_RT{fbus(k),tbus(k),s} : norm(A*test_vec{s}(:,k)) <= b*test_vec{s}(:,k);

			test_vec2{s}(:,k) = [p_flow(tbus(k),fbus(k),s); q_flow(tbus(k),fbus(k),s);...
				i_sq(tbus(k),fbus(k),s); v_sq(tbus(k),s)];
			gamma_RT{tbus(k),fbus(k),s} : norm(A*test_vec2{s}(:,k)) <= b*test_vec2{s}(:,k);
		end
	end

	for s = 1:nscen
		for k = 1:length(fbus)
			mu_p_RT{fbus(k),tbus(k),s} : p_flow(fbus(k),tbus(k),s) + p_flow(tbus(k),fbus(k),s) == rBranch(branch_num(k)).*i_sq(fbus(k),tbus(k),s);
			mu_p_RT{tbus(k),fbus(k),s} : p_flow(tbus(k),fbus(k),s) + p_flow(fbus(k),tbus(k),s) == rBranch(branch_num(k)).*i_sq(tbus(k),fbus(k),s);

			mu_q_RT{fbus(k),tbus(k),s} : q_flow(fbus(k),tbus(k),s) + q_flow(tbus(k),fbus(k),s) == xBranch(branch_num(k)).*i_sq(fbus(k),tbus(k),s);
			mu_q_RT{tbus(k),fbus(k),s} : q_flow(tbus(k),fbus(k),s) + q_flow(fbus(k),tbus(k),s) == xBranch(branch_num(k)).*i_sq(tbus(k),fbus(k),s);

			beta_RT{fbus(k),tbus(k),s} : v_sq(tbus(k),s) == v_sq(fbus(k),s)...
				- 2*(rBranch(branch_num(k)).*p_flow(fbus(k),tbus(k),s)...
				+ xBranch(branch_num(k)).*q_flow(fbus(k),tbus(k),s) ) ...
				+ (rBranch(branch_num(k))^2+xBranch(branch_num(k))^2).*i_sq(fbus(k),tbus(k),s);
		
			beta_RT{tbus(k),fbus(k),s} : v_sq(fbus(k),s) == v_sq(tbus(k),s)...
				- 2*(rBranch(branch_num(k)).*p_flow(tbus(k),fbus(k),s)...
				+ xBranch(branch_num(k)).*q_flow(tbus(k),fbus(k),s) ) ...
				+ (rBranch(branch_num(k))^2+xBranch(branch_num(k))^2).*i_sq(tbus(k),fbus(k),s);
		end
	end
	
	for s = 1:nscen
		sigma_RT_minus{s} : Vmin(buses).^2 <= v_sq(:,s) <= Vmax(buses).^2 : sigma_RT_plus{s};

		lambda_q_RT{s} : sum(q_flow(:,:,s),2) == q_inj_rt(:,s);

		q_inj_rt(:,s) == incidentG' * q_gen_RT(:,s) ...
					- incidentD' * q_dem_RT(:,s) + shed_q_RT(:,s) + incident_pcc * q_pcc_rt(s);


		for k = 1:length(fbus)
			if SlmMax(branch_num(k)) ~= 0
				p_flow(fbus(k),tbus(k),s).^2 + q_flow(fbus(k),tbus(k),s).^2 <= SlmMax(branch_num(k)).^2 : eta_RT{fbus(k),tbus(k),s};
			else
				eta_RT{fbus(k),tbus(k),s} : p_flow(fbus(k),tbus(k),s).^2 + q_flow(fbus(k),tbus(k),s).^2 <= 10000;
			end
		end
	end
	
	p_flow_matrix = reshape(sum(p_flow,2),size(p_inj_rt));
	lambda_RT : 0 == p_flow_matrix - p_inj_rt;
	for s = 1:nscen
		p_DR_RT(:,s) == p_DR_DA - pup_DR(:,s) + pdn_DR(:,s);
		zeta_pcc_RT{s} : p_pcc_rt(s) == p_pcc_DA + pup_pcc(s) - pdn_pcc(s);
		zeta_gen_RT{s} : p_gen_RT(:,s) == p_gen_DA + pup_g(:,s) - pdn_g(:,s);
		zeta_dem_RT{s} : p_dem_RT(:,s) == p_dem_DA - pup_d(:,s) + pdn_d(:,s);

		p_inj_rt(:,s) == incidentG' * p_gen_RT(:,s) - incidentD' * p_dem_RT(:,s) - incidentDR' * p_DR_RT(:,s)...
				+ incidentW' * wind_RT(:,s) + shed_p_RT(:,s) + incident_pcc * p_pcc_rt(s);
			
		Upsilon_RT_minus{s} : 0 <= shed_p_RT(:,s);
		Upsilon_RT_plus{s} : 0 <= (incidentD' * p_dem_RT(:,s)) - shed_p_RT(:,s);
		Upsilon_Q_RT_minus{s} : 0 <= shed_q_RT(:,s) <= (incidentD' * q_dem_RT(:,s));
	end
	

	for s = 1:nscen
		varsigma_gen_RT_minus{s} : Pmin_gen(gens_area) <= p_gen_RT(:,s) <= Pmax_gen(gens_area) : varsigma_gen_RT_plus{s};
		kappa_gen_RT_minus{s} : Qmin_gen(gens_area) <= q_gen_RT(:,s) <= Qmax_gen(gens_area) : kappa_gen_RT_plus{s};

		varsigma_dem_RT_minus{s} : Pmin_dem(dems_area) <= p_dem_RT(:,s) <= Pmax_dem(dems_area) : varsigma_dem_RT_plus{s};
		kappa_dem_RT_minus{s} : Qmin_dem(dems_area) <= q_dem_RT(:,s) <= Qmax_dem(dems_area) : kappa_dem_RT_plus{s};

		nu_RT_minus{s} : Wmin(windg_area,s) <= wind_RT(:,s) <= Wmax(windg_area,s) : nu_RT_plus{s};
		Pmin_DR(DR_area) <= p_DR_RT(:,s) <= Pmax_DR(DR_area);

	end
	
	rho_RT_minus : fe_dn <= p_pcc_rt <= fe_up : rho_RT_plus;
	
	
	
cvx_end
	
dual_objective_DA = sum(Pmin_gen(gens_area)' * varsigma_gen_DA_minus - Pmax_gen(gens_area)' * varsigma_gen_DA_plus) ...
	+ sum(Pmin_dem(dems_area)' * varsigma_dem_DA_minus - Pmax_dem(dems_area)' * varsigma_dem_DA_plus )...
	- iota_plus_DA * sum(Wmax_mean_DA(windg_area)) ...
	+ fe_dn * rho_DA_minus - fe_up * rho_DA_plus;
% 		- eta_RT_s(fbus,tbus,:) .* SlmMax(branch_num) ...
for s = 1:nscen
	dual_objective_RT(s) =  sum(Vmin(buses)'.^2 * sigma_RT_minus{s})  - sum(Vmax(buses)'.^2 * sigma_RT_plus{s}) ...
		- sum(nu_RT_plus{s}' *  Wmax(windg_area,s) )...
		+ sum(Pmin_gen(gens_area)' * varsigma_gen_RT_minus{s} - Pmax_gen(gens_area)' * varsigma_gen_RT_plus{s}) ...
		+ sum(Pmin_dem(dems_area)' * varsigma_dem_RT_minus{s} - Pmax_dem(dems_area)' * varsigma_dem_RT_plus{s}) ...
		+ sum(Qmin_gen(gens_area)' * kappa_gen_RT_minus{s} - Qmax_gen(gens_area)' * kappa_gen_RT_plus{s}) ...
		+ sum(Qmin_dem(dems_area)' * kappa_dem_RT_minus{s} - Qmax_dem(dems_area)' * kappa_dem_RT_plus{s}) ...
		+ fe_dn * rho_RT_minus(s) - fe_up * rho_RT_plus(s) ...
		;
end
dual_objective = dual_objective_DA + sum(dual_objective_RT);

if any(shed_p_RT > 10e-5)
	warning('There is load shedding in the DSO look-ahead Real Time calculation');
end
if any(shed_p_DA > 10e-5)
	warning('There is load shedding in the DSO look-ahead day-ahead calculation');
end
result_lookahead.cost = cost;
result_lookahead.cost_RT = cost_RT;
result_lookahead.p_dem_DA = p_dem_DA;
result_lookahead.p_DR_DA = p_DR_DA;
result_lookahead.p_gen_DA = p_gen_DA;
result_lookahead.wind_DA = wind_DA;
result_lookahead.p_pcc_DA = p_pcc_DA;
result_lookahead.p_pcc_RT = p_pcc_rt;

result_lookahead.p_DR_RT = p_DR_RT;

result_lookahead.lambda_p_RT = lambda_RT;
result_lookahead.lambda_q_RT = cell2mat(lambda_q_RT');
result_lookahead.lambda_DA = lambda_DA;
result_lookahead.zeta_p_dem_RT = cell2mat(zeta_dem_RT');
result_lookahead.zeta_p_gen_RT = cell2mat(zeta_gen_RT');
result_lookahead.zeta_pcc_RT = cell2mat(zeta_pcc_RT');
result_lookahead.varsigma_gen_RT_minus = cell2mat(varsigma_gen_RT_minus');
result_lookahead.varsigma_gen_RT_plus = cell2mat(varsigma_gen_RT_plus');
result_lookahead.varsigma_dem_RT_minus = cell2mat(varsigma_dem_RT_minus');
result_lookahead.varsigma_dem_RT_plus = cell2mat(varsigma_dem_RT_plus');
result_lookahead.rho_RT_minus = rho_RT_minus;
result_lookahead.rho_RT_plus = rho_RT_plus;
result_lookahead.rho_DA_minus = rho_DA_minus;
result_lookahead.rho_DA_plus = rho_DA_plus;

p_gen_DA
p_dem_DA
p_pcc_DA
p_gen_RT
p_dem_RT
p_pcc_rt
wind_DA
wind_RT
v_sq
shed_p_RT
shed_p_DA
i_sq
lambda_DA
lambda_RT
varsigma_gen_DA_minus
varsigma_gen_DA_plus
varsigma_dem_DA_minus
varsigma_dem_DA_plus
Upsilon_RT_minus = cell2mat(Upsilon_RT_minus')
Upsilon_RT_plus = cell2mat(Upsilon_RT_plus')
Upsilon_Q_RT_minus = cell2mat(Upsilon_Q_RT_minus')
lambda_q_RT=cell2mat(lambda_q_RT')
zeta_gen_RT=cell2mat(zeta_gen_RT')
zeta_dem_RT=cell2mat(zeta_dem_RT')
zeta_pcc_RT=cell2mat(zeta_pcc_RT')
varsigma_gen_RT_minus=cell2mat(varsigma_gen_RT_minus')
varsigma_gen_RT_plus=cell2mat(varsigma_gen_RT_plus')
varsigma_dem_RT_minus=cell2mat(varsigma_dem_RT_minus')
varsigma_dem_RT_plus=cell2mat(varsigma_dem_RT_plus')
nu_RT_minus = cell2mat(nu_RT_minus')
nu_RT_plus = cell2mat(nu_RT_plus')
iota_plus_DA
iota_minus_DA
for k = 1:numel(mu_p_RT)
	if isempty(mu_p_RT{k})
		mu_p_RT{k} = 0;
	end
end
mu_p_RT = cell2mat(mu_p_RT);
mu_p_RT
gamma_RT
% for k = 1:numel(beta_RT)
% 	if isempty(beta_RT{k})
% 		beta_RT{k} = 0;
% 	end
% end
% beta_RT = cell2mat(beta_RT);
% beta_RT = beta_RT - permute(beta_RT,[2,1,3])
beta_RT

for k = 1:numel(mu_q_RT)
	if isempty(mu_q_RT{k})
		mu_q_RT{k} = 0;
	end
end
mu_q_RT = cell2mat(mu_q_RT);
mu_q_RT
eta_RT
dual_objective
rho_RT_minus
rho_RT_plus
rho_DA_minus
rho_DA_plus

result_lookahead.dual_objective = dual_objective;
result_lookahead.mu_p_RT = mu_p_RT;
result_lookahead.mu_q_RT = mu_q_RT;
result_lookahead.beta_RT = beta_RT;
result_lookahead.p_flow = p_flow;
result_lookahead.q_flow = q_flow;
result_lookahead.s_flow = sqrt(p_flow.^2 + q_flow.^2);

end