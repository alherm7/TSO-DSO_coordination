function [RT_outcome] = RT_fixed_DA_dispatch(data_tso,p_gen_DA_hat,p_dem_DA_hat,p_DR_DA_hat,wind_DA_hat,s,iter,kk)

VOLL = 10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('results_full_coord.mat');
% fcoo_pg_RT = results_full_coordination.pg_RT;
% fcoo_pd_RT = results_full_coordination.pd_RT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

nw = length(windgL);

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
	Pmax_DR = double.empty(0,1);
	Pmin_DR = double.empty(0,1);
	bid_DA_DR = [];
	bid_up_DR = [];
	bid_dn_DR = [];

	DR_N = [];

	nDR = length(bus_DR);
	incidentDR =  sparse(1:nDR, DR_N, 1 , nDR, nb);
	incidentDR = full(incidentDR);
end

cvx_solver mosek


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

for g = 1:ng
	if abs(p_gen_DA_hat(g)) < 1e-4
		p_gen_DA_hat(g) = 0;
	end
end
for d = 1:nd
	if abs(p_dem_DA_hat(d)) < 1e-4
		p_dem_DA_hat(d) = 0;
	end
end
for w = 1:nw
	if abs(wind_DA_hat(w)) < 1e-4
		wind_DA_hat(w) = 0;
	end
end
%%
cvx_begin quiet
	
% 	cvx_precision high

	variable pup_DR(nDR) nonnegative
	variable pdn_DR(nDR) nonnegative
	variable pup_g(ng) nonnegative
	variable pdn_g(ng) nonnegative
	variable pup_d(nd) nonnegative
	variable pdn_d(nd) nonnegative
	variable shed_p(nb) nonnegative
	variable shed_q(nb) nonnegative
	variable wind_g(n_wgen) nonnegative
	variable p_gen_DA(ng) nonnegative
	variable p_dem_DA(nd) nonnegative
	variable p_DR_DA(nDR) nonnegative
	variable p_wind_DA(nw) nonnegative
	variable pup_wind(nw) nonnegative
	variable pdn_wind(nw) nonnegative

	variable p_inj_rt(nb)
	variable q_inj_rt(nb)
	variable p_flow(nb,nb)
	variable q_flow(nb,nb)
	variable p_DR(nDR)
	variable p_g(ng)
	variable p_d(nd)
	variable q_g(ng)
	variable q_d(nd)

	variable i_sq(nb,nb) nonnegative
	variable v_sq(nb) nonnegative

	variable v_angle(nb)

	expression cost_g(ng)
	expression cost_d(nd)
	expression cost_RT
	expression cost_DR(nDR)
	
	dual variable lambda
	dual variable gamma_price
	dual variable dual_day_ahead_dem
	dual variable dual_day_ahead_gen
	dual variable dual_day_ahead_DR
	dual variable dual_day_ahead_wind
	
% 	cost_g = 0;
% 	cost_d = 0;
% 	clear cost_g cost_d

	for m = 1:ng
		cost_g(m) = offer_gen_DA(m,2).*(p_g(m) - p_gen_DA(m)) ...
			+ (offer_gen_upreg(m))*pup_g(m)...
			+ (offer_gen_downreg(m)) *pdn_g(m);
	end

	for m = 1:nd
		cost_d(m) = bid_dem_DA(m).*(p_dem_DA(m) - p_d(m)) ...
			+ ( bid_dem_upreg(m))*pup_d(m) ...
			+ ( bid_dem_downreg(m)) *pdn_d(m);
	end

	for m = 1:nDR
		cost_DR(m) = bid_DA_DR(m).*(p_DR_DA(m) - p_DR(m)) ...
			+ ( bid_up_DR(m))*pup_DR(m) ...
			+ ( bid_dn_DR(m)) *pdn_DR(m);
	end
	
	for w = 1:nw
		cost_wind(w) = offer_wind_DA(w).*(wind_g(w) - p_wind_DA(w)) ...
			+ offer_wind_up(w) * pup_wind(w) ...
			+ offer_wind_dn(w) * pdn_wind(w);
	end


	cost_shed = VOLL*(shed_p+shed_q);

% 	cost_wind = wind_g*offer_wind;

	cost_RT = sum(cost_g) + sum(cost_d) + sum(cost_DR) + sum(cost_shed) + sum(cost_wind);

	minimize(cost_RT)


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subject to
	
	dual_day_ahead_gen : p_gen_DA == p_gen_DA_hat;
	dual_day_ahead_dem : p_dem_DA == p_dem_DA_hat;
	dual_day_ahead_DR :	p_DR_DA == p_DR_DA_hat;
	dual_day_ahead_wind : p_wind_DA(:) == wind_DA_hat;

	incident_contra = ones(nb,nb) - incident_singlephase;
	incident_contra = logical(incident_contra);

	ind_contra = find(incident_contra);
	p_flow(ind_contra) == 0;
	q_flow(ind_contra) == 0;
%% DSO PF
	for DSO_num = 1:ne

		branch_num = branch_num_ext_area{DSO_num+1};
% 		[~,~,gens_area] = intersect(ext_area{DSO_num+1},genL);
% 		[~,~,dems_area] = intersect(ext_area{DSO_num+1},demL);
% 		[~,~,windg_area] = intersect(ext_area{DSO_num+1},windgL);
		gens_area = find(ismember(genL,ext_area{DSO_num+1}));
		dems_area =  find(ismember(demL,ext_area{DSO_num+1}));
		windg_area = find(ismember(windgL,ext_area{DSO_num+1}));

		buses = ext_area{DSO_num+1};
		fbus = fbusL_ext_area{DSO_num+1};
		tbus = tbusL_ext_area{DSO_num+1};

		A = [2	0	0	0;
			0	2	0	0;
			0	0	1	-1];
		b = [0	0	1	1];
		
		test_vec = 0;
		test_vec2 = 0;
		clear test_vec test_vec2
		
		expression test_vec(4,length(fbus))
		expression test_vec2(4,length(fbus))

		for k = 1:length(fbus)
			test_vec(:,k) = [p_flow(fbus(k),tbus(k)); q_flow(fbus(k),tbus(k)); ...
				i_sq(fbus(k),tbus(k)); v_sq(fbus(k))];
			norm(A*test_vec(:,k)) <= b*test_vec(:,k);

			test_vec2(:,k) = [p_flow(tbus(k),fbus(k)); q_flow(tbus(k),fbus(k));...
				i_sq(tbus(k),fbus(k)); v_sq(tbus(k))];
			norm(A*test_vec2(:,k)) <= b*test_vec2(:,k);
		end




		for k = 1:length(fbus)
			p_flow(fbus(k),tbus(k)) + p_flow(tbus(k),fbus(k)) == rBranch(branch_num(k)).*i_sq(fbus(k),tbus(k));
			p_flow(tbus(k),fbus(k)) + p_flow(fbus(k),tbus(k)) == rBranch(branch_num(k)).*i_sq(tbus(k),fbus(k));

			q_flow(fbus(k),tbus(k)) + q_flow(tbus(k),fbus(k)) == xBranch(branch_num(k)).*i_sq(fbus(k),tbus(k));
			q_flow(tbus(k),fbus(k)) + q_flow(fbus(k),tbus(k)) == xBranch(branch_num(k)).*i_sq(tbus(k),fbus(k));

			v_sq(tbus(k)) == v_sq(fbus(k)) - 2*(rBranch(branch_num(k)).*p_flow(fbus(k),tbus(k))...
				+ xBranch(branch_num(k)).*q_flow(fbus(k),tbus(k))) ...
				+ (rBranch(branch_num(k))^2+xBranch(branch_num(k))^2).*i_sq(fbus(k),tbus(k));
		end

		Vmin(buses).^2/3 <= v_sq(buses) <= Vmax(buses).^2*3;

		sum(q_flow(buses,buses),2) == q_inj_rt(buses);

		q_inj_rt(buses) == incidentG(gens_area,buses)'*q_g(gens_area) ...
					- incidentD(dems_area,buses)'*q_d(dems_area) + shed_q(buses);


% 		for k = 1:length(fbus)
% 			if SlmMax(branch_num(k)) ~= 0
% 				p_flow(fbus(k),tbus(k)).^2 + q_flow(fbus(k),tbus(k)).^2 <= SlmMax(branch_num(k)).^2;
% 				p_flow(tbus(k),fbus(k)).^2 + q_flow(tbus(k),fbus(k)).^2 <= SlmMax(branch_num(k)).^2;
% 			end
% 		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% TSO power flow in Real time with DC-power flow
	branch_num = branch_num_area{1};
	fbus = fbusL_area{1};
	tbus = tbusL_area{1};
	buses =  areas{1};
% 	[~,~,gens_area] = intersect(buses,genL);
% 	[~,~,dems_area] = intersect(buses,demL);
% 	[~,~,windg_area] = intersect(buses,bus_wgen);

	for k = 1:length(fbus)
		p_flow(fbus(k),tbus(k)) == 1/xBranch(branch_num(k)) * (v_angle(fbus(k)) - v_angle(tbus(k)));
		p_flow(tbus(k),fbus(k)) == 1/xBranch(branch_num(k)) * (v_angle(tbus(k)) - v_angle(fbus(k)));
	end

	lambda : sum(p_flow,2) ==  p_inj_rt;

	p_g == p_gen_DA + pup_g - pdn_g;
	p_d == p_dem_DA - pup_d + pdn_d;
	p_DR == p_DR_DA - pup_DR + pdn_DR;
	wind_g == p_wind_DA + pup_wind - pdn_wind;
% 	fcoo_pg_RT(:,s) == p_g;
% 	fcoo_pd_RT(:,s) == p_d;
pup_g == 0;
pdn_g == 0;
pup_d == 0;
pdn_d == 0;
pup_DR == 0;
pdn_DR == 0;
pup_wind == 0;
pdn_wind == 0;


	gamma_price : p_inj_rt == incidentG' * p_g - incidentD' * p_d...
							+ incidentW' * wind_g + shed_p - incidentDR' * p_DR;


	Pmin_gen <= p_g(:) <= Pmax_gen;
	Qmin_gen <= q_g(:) <= Qmax_gen;

	Pmin_dem <= p_d(:) <= Pmax_dem;
	Qmin_dem <= q_d(:) <= Qmax_dem;
	
	Pmin_DR <= p_DR(:) <= Pmax_DR;

	v_angle(1) == 0;
% 	Wmin(:,s) <= wind_g <= Wmax(:,s);
	
% 	for l = 1:length(fbus)
% 		if SlmMax(branch_num(l)) ~= 0
% 			p_flow(fbus(l),tbus(l)) <= SlmMax(branch_num(l));
% 			p_flow(tbus(l),fbus(l)) <= SlmMax(branch_num(l));
% 		end
% 	end
	
	0 <= shed_p <= incidentD' * p_d;



cvx_end
%%
% if any(round(wind_g*10)/10 < round(Wmax(:,s)*10)/10)
% 	disp(['There is Wind Spillage in Scenario ' num2str(s) '. ' num2str(sum(Wmax(:,s)-wind_g)) ' MW of wind is spilled']);
% 	disp(['Wind Farm ' num2str(find(round(wind_g*10)/10 < round(Wmax(:,s)*10)/10)') ' is spilling power']);
% end
congestion = zeros(nb,nb);
s_flow = sqrt(p_flow.^2 + q_flow.^2);
for l = 1:nl
	if any(l == branch_num_area{1}) && ( SlmMax(l) - abs(p_flow(fbusL(l),tbusL(l))) < 0.5 ) && (SlmMax(l) ~= 0) && p_flow(fbusL(l),tbusL(l)) < 0
% 		disp(['Line ' num2str(l) ' (from node ' num2str(fbusL(l)) ' to node ' num2str(tbusL(l)) ' in TSO network ) ' ' in Scenario ' num2str(s) ' is congested'])
		congestion(fbusL(l),tbusL(l)) = -1;
	elseif any(l == branch_num_area{1}) && ( SlmMax(l) - abs(p_flow(fbusL(l),tbusL(l))) < 0.5 ) && (SlmMax(l) ~= 0) && p_flow(fbusL(l),tbusL(l)) > 0
		congestion(fbusL(l),tbusL(l)) = 1;
	elseif any(l == cell2mat(branch_num_ext_area(2:end) ) ) ...
			&& ( SlmMax(l) - abs(s_flow(fbusL(l),tbusL(l))) < 0.5 ) && (SlmMax(l) ~= 0) && s_flow(fbusL(l),tbusL(l)) < 0
% 		disp(['Line ' num2str(l) ' (from node ' num2str(fbusL(l)) ' to node ' num2str(tbusL(l)) ' in DSO feeder) ' ' in Scenario ' num2str(s) ' is congested'])
		congestion(fbusL(l),tbusL(l)) = -1;
	elseif any(l == cell2mat(branch_num_ext_area(2:end) ) ) ...
			&& ( SlmMax(l) - abs(s_flow(fbusL(l),tbusL(l))) < 0.5 ) && (SlmMax(l) ~= 0) && s_flow(fbusL(l),tbusL(l)) > 0
		congestion(fbusL(l),tbusL(l)) = 1;
	end
end
	
	
if any(shed_p > 10e-5)
	warning(['WPP: ' num2str(kk) ', Benders sub-prob. It.: ' num2str(iter) ', Scenario: ' num2str(s) ', There is load shedding in the Real Time market']);
% 		pause(3)
end

% for s = 1:nscen
	for k = 1:size(test_vec,2)
		tight1(k) = ( b*test_vec(:,k) - norm(A*test_vec(:,k))  );
		tight2(k) = ( b*test_vec2(:,k) - norm(A*test_vec2(:,k)) );
% 		if any( tight1(k) > 10e-6) ...
% 				|| any( tight2(k) > 10e-6)
% 			warning(['SOCP relaxation of scenario ' num2str(s) ' and line ' num2str(k) ' is not tight']);
% 			pause(0.01)
% 		end
	end
% end
wind_penetration_actual = sum(wind_g)/sum(p_d);
wind_penetration_total = sum(Wmax(:,s))/sum(p_d);
wind_penetration_offered = sum(Wmax(:,s))/sum(Pmax_dem);

RT_outcome.tight1 = tight1;
RT_outcome.tight2 = tight2;

RT_outcome.cost_RT = cost_RT;
RT_outcome.cost_gen_RT = cost_g;
RT_outcome.cost_dem_RT = cost_d;
RT_outcome.cost_wind_RT = cost_wind;
RT_outcome.p_flow = full(p_flow);
RT_outcome.q_flow = full(q_flow);

RT_outcome.p_d = p_d;
RT_outcome.p_g = p_g;
RT_outcome.p_DR = p_DR;
RT_outcome.pup_g = pup_g;
RT_outcome.pdn_g = pdn_g;
RT_outcome.pup_d = pup_d;
RT_outcome.pdn_d = pdn_d;
RT_outcome.pup_DR = pup_DR;
RT_outcome.pdn_DR = pdn_DR;
RT_outcome.v_sq = v_sq;
RT_outcome.congestion = congestion;
RT_outcome.shed_p = shed_p;
RT_outcome.wind_penetration_actual = wind_penetration_actual;
RT_outcome.wind_penetration_total = wind_penetration_total;
RT_outcome.wind_penetration_offered = wind_penetration_offered;


RT_outcome.dual_day_ahead_gen = dual_day_ahead_gen;
RT_outcome.dual_day_ahead_dem = dual_day_ahead_dem;
RT_outcome.dual_day_ahead_DR = dual_day_ahead_DR;
RT_outcome.dual_day_ahead_wind = dual_day_ahead_wind;
RT_outcome.wind_RT = wind_g;
end

