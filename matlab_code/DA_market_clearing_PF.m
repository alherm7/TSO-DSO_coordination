function [DA_market_outcome] = DA_market_clearing_PF(data_tso,dual_DA_gen,dual_DA_dem,...
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
gens_T = find(ismember(genL,areas{1}));
dems_T = find(ismember(demL,areas{1}));

nw = length(windgL);

%% define variables
cvx_begin quiet

% cvx_precision high

variable p_gen_DA(ng)
variable p_dem_DA(nd) nonnegative
variable wind_DA(nw) nonnegative
variable shed_p(nb) nonnegative
variable shed_q(nb) nonnegative
variable p_inj_rt(nb)
variable q_inj_rt(nb)


variable alpha_cut(nscen)
variable p_flow(nb,nb)
variable q_flow(nb,nb)
variable i_sq(nb,nb) nonnegative
variable v_sq(nb) nonnegative
variable q_g(ng)
variable q_d(nd)
    
dual variable lambda
dual variable lambda_2
dual variable gamma_price
variable v_angle(nb)

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
cost_shed = VOLL_DA*sum(shed_p);

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
% lambda : sum(p_gen_DA) + sum(wind_DA) + sum(shed_p) == sum(p_dem_DA);

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

%% DSO PF
incident_contra = ones(nb,nb) - incident_singlephase;
incident_contra = logical(incident_contra);

ind_contra = find(incident_contra);
p_flow(ind_contra) == 0;
q_flow(ind_contra) == 0;

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

    Vmin(buses).^2 <= v_sq(buses) <= Vmax(buses).^2;

    sum(q_flow(buses,buses),2) == q_inj_rt(buses);

    q_inj_rt(buses) == incidentG(gens_area,buses)'*q_g(gens_area) ...
                - incidentD(dems_area,buses)'*q_d(dems_area) + shed_q(buses);


    for k = 1:length(fbus)
        if SlmMax(branch_num(k)) ~= 0
            p_flow(fbus(k),tbus(k)).^2 + q_flow(fbus(k),tbus(k)).^2 <= SlmMax(branch_num(k)).^2;
            p_flow(tbus(k),fbus(k)).^2 + q_flow(tbus(k),fbus(k)).^2 <= SlmMax(branch_num(k)).^2;
        end
    end
end

%%%%%%%%%%%%%%%%%
%% DC power flow

branch_num = branch_num_area{1};
fbus = fbusL_area{1};
tbus = tbusL_area{1};


for k = 1:length(fbus)
    p_flow(fbus(k),tbus(k)) == 1/xBranch(branch_num(k)) * (v_angle(fbus(k)) - v_angle(tbus(k)));
    p_flow(tbus(k),fbus(k)) == 1/xBranch(branch_num(k)) * (v_angle(tbus(k)) - v_angle(fbus(k)));
end

lambda_2 : sum(p_flow,2) ==  p_inj_rt;


% 	fcoo_pg_RT(:,s) == p_g;
% 	fcoo_pd_RT(:,s) == p_d;

gamma_price : p_inj_rt == incidentG' * p_gen_DA - incidentD' * p_dem_DA...
                        + incidentW' * wind_DA + shed_p;


Pmin_gen <= p_gen_DA(:) <= Pmax_gen;

Pmin_dem <= p_dem_DA(:) <= Pmax_dem;


v_angle(1) == 0;

for l = 1:length(fbus)
    if SlmMax(branch_num(l)) ~= 0
        p_flow(fbus(l),tbus(l)) <= SlmMax(branch_num(l));
        p_flow(tbus(l),fbus(l)) <= SlmMax(branch_num(l));
    end
end

0 <= shed_p <= incidentD' * p_dem_DA;

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

