function [cooptim_outcome] = TSO_DSO_cooptim(data_tso)

%% load parameters
VOLL_T = 10000;
VOLL_E = 10000;
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
	nw,nscen,Wmax_mean_DA,windgL,...
	offer_wind_DA,offer_wind_up,offer_wind_dn] = Data_Reader(data_tso,cc);


cvx_solver mosek


[area_codes, n_areas, DSO_codes, fbusL_area, tbusL_area, branch_num_area, ...
	fbusL_ext_area, tbusL_ext_area, branch_num_ext_area, overlap,...
	overlap_t, neigh_area, ext_area, ext_area_singlephase,PCC_branch_id] ...
	= find_overlap(bus_areacodes, fbusL, tbusL, areas);

windg_area_T = find(ismember(windgL,areas{1}));

% nw = length(windgL);

ne = n_areas-1; %#ok<NASGU>


gens_DSO = [];
dems_DSO = [];
windg_DSO = [];
fbus_DSO = [];
tbus_DSO = [];
buses_DSO = [];
branch_num_DSO = [];

for DSO_num = 1:length(DSO_codes)
	gens_area{DSO_num} = find(ismember(genL,areas{DSO_num+1}));
	dems_area{DSO_num} =  find(ismember(demL,areas{DSO_num+1}));
	windg_area{DSO_num} = find(ismember(windgL,areas{DSO_num+1}));
	
	gens_DSO = [gens_DSO; gens_area{DSO_num}];
	dems_DSO = [dems_DSO; dems_area{DSO_num}];
	windg_DSO = [windg_DSO; windg_area{DSO_num}];
	
	buses = ext_area{DSO_num+1};
	fbus = fbusL_ext_area{DSO_num+1};
	tbus = tbusL_ext_area{DSO_num+1};
	branch_num = branch_num_ext_area{DSO_num+1};

	
	fbus_DSO = [fbus_DSO; fbus];
	tbus_DSO = [tbus_DSO; tbus];
	buses_DSO = [buses_DSO; buses];
	branch_num_DSO = [branch_num_DSO; branch_num(:)];

end


ng_DSO = length(gens_DSO);
nd_DSO = length(dems_DSO);

%% begin with cvx and declare variables
build_time = tic;
cvx_begin quiet

% cvx_precision high

variable pup_g(ng,nscen) nonnegative
variable pdn_g(ng,nscen) nonnegative
variable pup_d(nd,nscen) nonnegative
variable pdn_d(nd,nscen) nonnegative
variable shed_p_RT(nb,nscen) nonnegative
variable shed_q_RT(nb,nscen) nonnegative
variable wind_g_RT(nw,nscen) nonnegative
variable pup_wind(nw,nscen) nonnegative
variable pdn_wind(nw,nscen) nonnegative
	
variable p_inj_rt(nb,nscen)
variable q_inj_rt(nb,nscen)
variable p_flow(nb,nb,nscen)
variable q_flow(nb,nb,nscen)
variable pg_RT(ng,nscen)
variable pd_RT(nd,nscen)
variable qg_RT(ng,nscen)
variable qd_RT(nd,nscen)
variable i_sq(nb,nb,nscen) nonnegative
variable v_sq(nb,nscen) nonnegative
variable v_angle(nb,nscen)

dual variable lambda_RT
dual variable gamma_RT

expression cost_g_RT(ng,nscen)
expression cost_d_RT(nd,nscen)
expression cost_RT(nscen)

variable p_gen_DA(ng)
variable p_dem_DA(nd) nonnegative
% variable p_wind_T_DA nonnegative
% variable p_wind_E_DA(ne) nonnegative
variable p_wind_DA(nw) nonnegative
variable s_shed_E_DA(ne) nonnegative
variable s_shed_T_DA nonnegative

dual variable lambda_DA

expression cost_gen_DA(ng)
expression cost_dem_DA(nd)

% variable pf_lin(nb,nb,t,s)
% variable qf_lin(nb,nb,t,s)
%% Day ahead objective function
for m = 1:ng
	if cost_type(m) == 2
		for k = 1:nr_vals_cost
			cost_gen_DA(m) = cost_gen_DA(m) + sum(offer_gen_DA(m,k).*p_gen_DA(m).^(nr_vals_cost-k)); %#ok<AGROW>
		end
	elseif cost_type(m) == 1
		error('Wrong cost function, piecewise linear is not convex');
	end
end

for m = 1:nd
	cost_dem_DA(m) = bid_dem_DA(m)*p_dem_DA(m); %#ok<AGROW>
end

cost_shed_T_DA = s_shed_T_DA*VOLL_T;
cost_shed_E_DA = sum(s_shed_E_DA*VOLL_E);
cost_shed_DA = cost_shed_E_DA+cost_shed_T_DA;


cost_wind = p_wind_DA.*offer_wind_DA;
cost_wind_e_DA = sum(p_wind_DA(windg_DSO).*offer_wind_DA(windg_DSO));
cost_wind_T_DA = p_wind_DA(windg_area_T).*offer_wind_DA(windg_area_T);

cost_DA = sum(cost_gen_DA)-sum(cost_dem_DA) + cost_shed_DA + sum(cost_wind);

%% Real time objective function
for m = 1:ng
	cost_g_RT(m,:) = offer_gen_DA(m,2).*(pg_RT(m,:) - p_gen_DA(m)) ...
		+ (offer_gen_upreg(m))*pup_g(m,:)...
		+ (offer_gen_downreg(m)) *pdn_g(m,:); %#ok<IDISVAR,NODEF,AGROW>
end

for m = 1:nd
	cost_d_RT(m,:) = bid_dem_DA(m).*(p_dem_DA(m) - pd_RT(m,:)) ...
		+ ( bid_dem_upreg(m))*pup_d(m,:) ...
		+ ( bid_dem_downreg(m)) *pdn_d(m,:); %#ok<IDISVAR,NODEF,AGROW>
end

cost_shed_RT = VOLL*(shed_p_RT+shed_q_RT);

% cost_wind_RT = wind_g_RT*offer_wind_DA;
for w = 1:nw
	cost_wind_RT(w,:) = offer_wind_DA(w).*(wind_g_RT(w,:) - p_wind_DA(w)) ...
		+  offer_wind_up(w) * pup_wind(w,:) ...
		+  offer_wind_dn(w) * pdn_wind(w,:);
end

cost_RT = sum(cost_g_RT) + sum(cost_d_RT) + sum(cost_shed_RT) + sum(cost_wind_RT,1);
cost_total = cost_DA+prob_wscen(1,:)*cost_RT';
%% objective function statement
minimize(cost_DA+prob_wscen(1,:)*cost_RT')
%% constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject to
%% Day ahead market
lambda_DA : sum(p_gen_DA) + sum(p_wind_DA) + sum(s_shed_E_DA) + s_shed_T_DA == sum(p_dem_DA); %#ok<EQEFF>

% for DSO_num = 1:length(DSO_codes)
	
% 	[~,~,gens_area] = intersect(areas{DSO_num+1},genL);
% 	[~,~,dems_area] = intersect(areas{DSO_num+1},demL);
% 	[~,~,windg_area] = intersect(areas{DSO_num+1},bus_wgen);
	
% 	for k = 1:ng
% 		if k == gens_area
			Pmin_gen <= p_gen_DA <= Pmax_gen; %#ok<VUNUS>
% 		else
% 			Pmin_gen(k) <= p_gen_DA(k) <= Pmax_gen(k); %#ok<VUNUS>
% 		end
% 	end
	
% 	for k = 1:nd
% 		for m = 1:length(dems_area)
% 			if k == dems_area(m)
% 				Pmin_dem(k) <= p_dem_DA(k) <= Pmax_dem(k); %#ok<VUNUS>
% 			else
				Pmin_dem <= p_dem_DA <= Pmax_dem; %#ok<VUNUS>
% 			end
% 		end
% 	end
	

% 	0 <= p_wind_E_DA(DSO_num) <= sum(Wmax_mean_DA(windg_DSO)); %#ok<VUNUS,CHAIN>
% end



0 <= p_wind_DA <= Wmax_mean_DA; %#ok<VUNUS,CHAIN>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% expected outcome of the real time market:
incident_contra = zeros(nb,nb,nscen);
for s = 1:nscen
	incident_contra(:,:,s) = ones(nb,nb) - incident_singlephase;
end
incident_contra = logical(incident_contra);


ind_contra = find(incident_contra);
p_flow(ind_contra) == 0;
q_flow(ind_contra) == 0;
%% DSO real time
% for DSO_num = 1:length(DSO_codes)

	branch_num = branch_num_DSO;
% 	[~,~,gens_area] = intersect(ext_area{DSO_num+1},genL);
% 	[~,~,dems_area] = intersect(ext_area{DSO_num+1},demL);
% 	[~,~,windg_area] = intersect(ext_area{DSO_num+1},bus_wgen);
	gens_area = gens_DSO;
	dems_area = dems_DSO;
	windg_area = windg_DSO;
	
% 	buses = ext_area{DSO_num+1};
% 	fbus = fbusL_ext_area{DSO_num+1};
% 	tbus = tbusL_ext_area{DSO_num+1};
	buses = buses_DSO;
	tbus = tbus_DSO;
	fbus = fbus_DSO;

	A = [2	0	0	0;
		0	2	0	0;
		0	0	1	-1];
	b = [0	0	1	1];
% 	expression test_vec(4,length(fbus))
% 	expression test_vec2(4,length(fbus))

	for k = 1:length(fbus)
		for s = 1:nscen
			test_vec{s}(:,k) = [p_flow(fbus(k),tbus(k),s); q_flow(fbus(k),tbus(k),s); ...
				i_sq(fbus(k),tbus(k),s); v_sq(fbus(k),s)];
			norm(A*test_vec{s}(:,k)) <= b*test_vec{s}(:,k);

			test_vec2{s}(:,k) = [p_flow(tbus(k),fbus(k),s); q_flow(tbus(k),fbus(k),s);...
				i_sq(tbus(k),fbus(k),s); v_sq(tbus(k),s)];
			norm(A*test_vec2{s}(:,k)) <= b*test_vec2{s}(:,k);
		end
	end



% 	for s = 1:nscen
		for k = 1:length(fbus)
			p_flow(fbus(k),tbus(k),:) + p_flow(tbus(k),fbus(k),:) == rBranch(branch_num(k)).*i_sq(fbus(k),tbus(k),:);
			p_flow(tbus(k),fbus(k),:) + p_flow(fbus(k),tbus(k),:) == rBranch(branch_num(k)).*i_sq(tbus(k),fbus(k),:);

			q_flow(fbus(k),tbus(k),:) + q_flow(tbus(k),fbus(k),:) == xBranch(branch_num(k)).*i_sq(fbus(k),tbus(k),:);
			q_flow(tbus(k),fbus(k),:) + q_flow(fbus(k),tbus(k),:) == xBranch(branch_num(k)).*i_sq(tbus(k),fbus(k),:);

			v_sq(tbus(k),:) == v_sq(fbus(k),:) - reshape( 2*(rBranch(branch_num(k)).*p_flow(fbus(k),tbus(k),:)...
				+ xBranch(branch_num(k)).*q_flow(fbus(k),tbus(k),:) ) ...
				+ (rBranch(branch_num(k))^2+xBranch(branch_num(k))^2).*i_sq(fbus(k),tbus(k),:), [1 nscen]);
		end
% 	end
	
	repmat(Vmin(buses),1,nscen).^2 <= v_sq(buses,:) <= repmat(Vmax(buses),1,nscen).^2;
	
	q_flow_matrix = reshape(sum(q_flow(buses,buses,:),2), size(q_inj_rt(buses,:)) );

	q_flow_matrix == q_inj_rt(buses,:);
	q_inj_rt(buses,:) == incidentG(gens_area,buses)'*qg_RT(gens_area,:) ...
					- incidentD(dems_area,buses)'*qd_RT(dems_area,:) + shed_q_RT(buses,:);
	for k = 1:length(fbus)
		if SlmMax(branch_num(k)) ~= 0
			p_flow(fbus(k),tbus(k),:).^2 + q_flow(fbus(k),tbus(k),:).^2 <= SlmMax(branch_num(k)).^2;
			p_flow(tbus(k),fbus(k),:).^2 + q_flow(tbus(k),fbus(k),:).^2 <= SlmMax(branch_num(k)).^2;
		end
	end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TSO power flow in Real time with DC-power flow
branch_num = branch_num_area{1};
fbus = fbusL_area{1};
tbus = tbusL_area{1};
% buses =  areas{1};
% [~,~,gens_area] = intersect(buses,genL);
% [~,~,dems_area] = intersect(buses,demL);
% [~,~,windg_area] = intersect(buses,bus_wgen);

for k = 1:length(fbus)
% 	for s = 1:nscen
		reshape( p_flow(fbus(k),tbus(k),:), [1 nscen] ) == 1/xBranch(branch_num(k)) * (v_angle(fbus(k),:) - v_angle(tbus(k),:) );
		reshape( p_flow(tbus(k),fbus(k),:), [1 nscen] ) == 1/xBranch(branch_num(k)) * (v_angle(tbus(k),:) - v_angle(fbus(k),:) );
% 	end
end

p_flow_matrix = reshape(sum(p_flow,2),size(p_inj_rt));

lambda_RT : p_flow_matrix ==  p_inj_rt;
for s = 1:nscen
	pg_RT(:,s) == p_gen_DA + pup_g(:,s) - pdn_g(:,s);
	pd_RT(:,s) == p_dem_DA - pup_d(:,s) + pdn_d(:,s);
	wind_g_RT(:,s) == p_wind_DA + pup_wind(:,s) - pdn_wind(:,s);
end

p_inj_rt == incidentG' * pg_RT - incidentD' * pd_RT...
		+ incidentW' * wind_g_RT + shed_p_RT;

% for s = 1:nscen
	repmat( Pmin_gen,1 ,nscen ) <= pg_RT <= repmat( Pmax_gen, 1, nscen);
	repmat( Qmin_gen,1 ,nscen ) <= qg_RT <= repmat( Qmax_gen,1 ,nscen );

	repmat( Pmin_dem,1 ,nscen ) <= pd_RT <= repmat( Pmax_dem,1 ,nscen );
	repmat( Qmin_dem,1 ,nscen ) <= qd_RT <= repmat( Qmax_dem,1 ,nscen );

	v_angle(1,:) == 0;
	Wmin <= wind_g_RT <= Wmax;
	
% end
% for s = 1:nscen
	for l = 1:length(fbus)
		if SlmMax(branch_num(l)) ~= 0
			p_flow(fbus(l),tbus(l),:) <= SlmMax(branch_num(l));
			p_flow(tbus(l),fbus(l),:) <= SlmMax(branch_num(l));
		end
	end
% end
t_overhead = toc(build_time);
cvx_end
%% output variables
voltage = sqrt(v_sq);
wind_penetration_actual = sum(wind_g_RT)./sum(pd_RT);
wind_penetration_total = sum(Wmax)./sum(pd_RT);
wind_penetration_offered = sum(Wmax)/sum(Pmax_dem);
congestion = cell(nscen,1);
s_flow = sqrt(p_flow.^2 + q_flow.^2);
for s = 1:nscen
	congestion{s} = zeros(nb,nb);
	for l = 1:nl
		if any(l == branch_num_area{1}) && ( SlmMax(l) - abs(p_flow(fbusL(l),tbusL(l),s)) < 0.5 ) && (SlmMax(l) ~= 0) && p_flow(fbusL(l),tbusL(l),s) < 0
	% 		disp(['Line ' num2str(l) ' (from node ' num2str(fbusL(l)) ' to node ' num2str(tbusL(l)) ' in TSO network ) ' ' in Scenario ' num2str(s) ' is congested'])
			congestion{s}(fbusL(l),tbusL(l)) = -1;
		elseif any(l == branch_num_area{1}) && ( SlmMax(l) - abs(p_flow(fbusL(l),tbusL(l),s)) < 0.5 ) && (SlmMax(l) ~= 0) && p_flow(fbusL(l),tbusL(l),s) > 0
			congestion{s}(fbusL(l),tbusL(l)) = 1;
		elseif any(l == cell2mat(branch_num_ext_area(2:end) ) ) ...
				&& ( SlmMax(l) - abs(s_flow(fbusL(l),tbusL(l),s)) < 0.5 ) && (SlmMax(l) ~= 0) && s_flow(fbusL(l),tbusL(l),s) < 0
	% 		disp(['Line ' num2str(l) ' (from node ' num2str(fbusL(l)) ' to node ' num2str(tbusL(l)) ' in DSO feeder) ' ' in Scenario ' num2str(s) ' is congested'])
			congestion{s}(fbusL(l),tbusL(l)) = -1;
		elseif any(l == cell2mat(branch_num_ext_area(2:end) ) ) ...
				&& ( SlmMax(l) - abs(s_flow(fbusL(l),tbusL(l),s)) < 0.5 ) && (SlmMax(l) ~= 0) && s_flow(fbusL(l),tbusL(l),s) > 0
			congestion{s}(fbusL(l),tbusL(l)) = 1;
		end
	end
end

cooptim_outcome.pg_RT = pg_RT;
cooptim_outcome.pd_RT = pd_RT;
cooptim_outcome.cost_DA = cost_DA;
cooptim_outcome.cost_RT = cost_RT;
cooptim_outcome.cost_wind_RT = cost_wind_RT;
cooptim_outcome.p_gen_DA = p_gen_DA;
cooptim_outcome.p_dem_DA = p_dem_DA;
cooptim_outcome.p_flow = p_flow;
cooptim_outcome.v_sq = v_sq;
cooptim_outcome.voltage = voltage;
cooptim_outcome.s_shed_E_DA = s_shed_E_DA;
cooptim_outcome.s_shed_T_DA = s_shed_T_DA;
cooptim_outcome.shed_p_RT = shed_p_RT;
cooptim_outcome.wind_g_RT = wind_g_RT;
% cooptim_outcome.p_wind_T_DA = p_wind_T_DA;
% cooptim_outcome.p_wind_E_DA = p_wind_E_DA;
cooptim_outcome.p_wind_DA = p_wind_DA;
cooptim_outcome.cost_total = cost_total;
cooptim_outcome.wind_penetration_actual = wind_penetration_actual;
cooptim_outcome.wind_penetration_total = wind_penetration_total;
cooptim_outcome.wind_penetration_offered = wind_penetration_offered;
cooptim_outcome.congestion = congestion;
cooptim_outcome.s_flow = s_flow;

for s = 1:nscen
	for k = 1:size(test_vec{s},2)
		tight1(k,s) = ( b*test_vec{s}(:,k) - norm(A*test_vec{s}(:,k))  );
		tight2(k,s) = ( b*test_vec2{s}(:,k) - norm(A*test_vec2{s}(:,k)) );
% 		if any( tight1(k,s) > 10e-6) ...
% 				|| any( tight2(k,s) > 10e-6)
% 			warning(['SOCP relaxation of scenario ' num2str(s) ' and line ' num2str(k) ' is not tight']);
% 		end
	end
end

cooptim_outcome.tight1 = tight1;
cooptim_outcome.tight2 = tight2;
