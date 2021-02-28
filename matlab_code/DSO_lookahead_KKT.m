function [result_lookahead] = DSO_lookahead_KKT(data,fe_up,fe_dn,...
	pi_pcc_DA,pi_pcc_up,pi_pcc_dn,e,npart_in,direc,relax)



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

buses = areas{e+1};
nb = length(areas{e+1});
ng = length(gens_area);
nd = length(dems_area);
nw = length(windg_area);
nl = length(fbus_local{e+1});


incidentG = incidentG(gens_area,buses);
incidentD = incidentD(dems_area,buses);
incidentW = full(incidentW(windg_area,buses));
incident_pcc = zeros(nb,1);
incident_pcc(pcc_local) = 1;

cvx_begin %quiet
	cvx_solver gurobi

% 	cvx_precision high
	
	%% primal variables
% 	wind_DA = result_lookahad_s.wind_DA;
		
	variable p_gen_DA(ng)
	variable p_dem_DA(nd) nonnegative
	variable wind_DA nonnegative
	variable shed_p_DA nonnegative
	variable p_pcc_DA

	
	variable pup_g(ng,nscen) nonnegative
	variable pdn_g(ng,nscen) nonnegative
	variable pup_d(nd,nscen) nonnegative
	variable pdn_d(nd,nscen) nonnegative
	variable p_pcc_rt(ne,nscen)
	variable q_pcc_rt(nscen)
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

	expression cost_g_RT(ng,nscen)
	expression cost_d_RT(nd,nscen)
	expression cost_RT
    
    variable p_pcc_comp(ne,nscen) binary
	variable p_gen_comp(ng,nscen) binary
	variable p_dem_comp(nd,nscen) binary
	
	%% dual variables
% 	lambda_q_RT = -result_lookahad_s.lambda_q_RT;
% 	lambda_p_RT = result_lookahad_s.lambda_p_RT;
% 	zeta_p_gen_RT = result_lookahad_s.zeta_p_gen_RT;
% 	zeta_p_dem_RT = result_lookahad_s.zeta_p_dem_RT;
% 	mu_p_RT = result_lookahad_s.mu_p_RT;
% 	mu_q_RT = result_lookahad_s.mu_q_RT;
% 	beta_RT = -result_lookahad_s.beta_RT;
	
	variable lambda_DA_e
	variable varsigma_gen_DA_minus_e(ng) nonnegative
	variable varsigma_gen_DA_plus_e(ng) nonnegative
	variable varsigma_dem_DA_minus_e(nd) nonnegative
	variable varsigma_dem_DA_plus_e(nd) nonnegative
	variable iota_minus_DA nonnegative
	variable iota_plus_DA nonnegative
	variable rho_DA_minus nonnegative
	variable rho_DA_plus nonnegative
	variable Upsilon_DA_minus nonnegative
	variable Upsilon_DA_plus nonnegative
	
	
	variable zeta_p_gen_RT(ng,nscen)
	variable zeta_p_dem_RT(nd,nscen)
	variable zeta_pcc_RT(nscen)
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
	variable epsilon_pcc_up_RT(nscen) nonnegative
	variable epsilon_pcc_dn_RT(nscen) nonnegative
	
	variable Upsilon_RT_minus(nb,nscen) nonnegative
	variable Upsilon_RT_plus(nb,nscen) nonnegative
	
	variable Upsilon_Q_RT_minus(nb,nscen) nonnegative
	variable Upsilon_Q_RT_plus(nb,nscen) nonnegative
	
	%% McCormick Variables for the linear envelopes
	variable w_piDA_pRT(ne,nscen)
	variable w_piDA_pDA(ne)
	variable w_piup_pup(ne,nscen)
	variable w_pidn_pdn(ne,nscen)

	variable w_fdn_rhominus_DA(ne)
	variable w_fup_rhoplus_DA(ne)
	variable w_fdn_rhominus_RT(ne,nscen)
	variable w_fup_rhoplus_RT(ne,nscen)
	
	%% binary variables for the complimentarity constraints
% 	variable bin_gen(ng,2) binary
% 	variable bin_dem(nd,2) binary
% 	variable bin_windlim_DA(2) binary
% 	variable bin_pcclim_DA(2) binary
% 	variable bin_socpflow(4,nb,nb,nscen) binary
% 	variable bin_flowlim(2,nb,nb,nscen) binary
% 	variable bin_voltlim(2,nb,nscen) binary
% 	variable bin_windlim_RT(2,nw,nscen) binary
% 	variable bin_pgenlim_RT(2,ng,nscen) binary
% 	variable bin_pdemlim_RT(2,nd,nscen) binary
% 	variable bin_qgenlim_RT(2,ng,nscen) binary
% 	variable bin_qdemlim_RT(2,nd,nscen) binary
% 	variable bin_pcclim_RT(2,nscen) binary
% 	variable bin_gen_RT(2,ng,nscen) binary
% 	variable bin_dem_RT(2,nd,nscen) binary
% 	variable bin_pcc_RT(2,nscen) binary
% 	variable bin_shed_RT(2,nb,nscen) binary
% 	variable bin_shed_Q_RT(2,nb,nscen) binary
% 	variable bin_shed_DA(2) binary
	
	
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
	
    if strcmp(relax,'none')
        cost_PCC_RT = pi_pcc_DA * (p_pcc_rt - p_pcc_DA) ...
                + pi_pcc_up * pup_pcc...
                + pi_pcc_dn * pdn_pcc;
    elseif strcmp(relax,'McCormick')
    	cost_PCC_RT = w_piDA_pRT(e,:) - w_piDA_pDA(e) ...
    			+ w_piup_pup(e,:) + w_pidn_pdn(e,:);
    end
	cost_RT = sum(cost_g_RT,1) + sum(cost_d_RT,1) + sum(cost_shed_RT,1) + sum(cost_wind_RT,1)+cost_PCC_RT;
	cost_RT_expec = prob_wscen(1,:)*cost_RT';
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

	cost_shed_DA = VOLL*shed_p_DA;
	cost_wind_DA = wind_DA*offer_wind;
	if strcmp(relax,'none')
        cost_PCC_DA = pi_pcc_DA * p_pcc_DA;
    elseif strcmp(relax,'McCormick')
        cost_PCC_DA = w_piDA_pDA(e);
    end
	cost_DA = sum(cost_gen_DA)-sum(cost_dem_DA) + cost_shed_DA + cost_wind_DA + cost_PCC_DA;
%% cost function statement
	cost = cost_RT_expec + cost_DA;
	if strcmp(direc,'min')
		minimize(sum(p_dem_DA(2)))
	elseif strcmp(direc,'max')
		maximize(sum(p_dem_DA(2)))
    elseif strcmp(direc,'dummy')
        minimize(1)
    end
	
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% primal constraints	
	subject to
	%% McCormick expansion of bi-linear terms
		npart = npart_in;
		
		variable p_pcc_rt_part(ne,npart,nscen) binary
		variable pi_pcc_DA_part(ne,npart) binary
		variable pi_pcc_up_part(ne,npart,nscen) binary
		variable pdn_pcc_part(ne,npart,nscen) binary

		
		for s = 1:nscen
			sum(p_pcc_rt_part(e,:,s)) == 1;
			sum(pi_pcc_up_part(e,:,s)) == 1;
			sum(pdn_pcc_part(e,:,s)) == 1;
		end
		sum(pi_pcc_DA_part(e,:)) == 1;
		
		
		pi_pcc_DA_lo = 1;
		pi_pcc_DA_up = 35;
		
		p_pcc_rt_lo = -10;
		p_pcc_rt_up = 10;
		
		p_pcc_DA_E_lo = -10;
		p_pcc_DA_E_up = 10;
		
		pi_pcc_up_lo = 1;
		pi_pcc_up_up = 14;
		
		pup_pcc_lo = 0;
		pup_pcc_up = 8;
		
		pi_pcc_dn_lo = 1;
		pi_pcc_dn_up = 14;
		
		pdn_pcc_lo = 0;
		pdn_pcc_up = 8;
		
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

		variable pi_pcc_DA_hat(ne,npart,nscen)
		variable p_pcc_DA_hat(ne,npart)
		variable pup_pcc_hat(ne,npart,nscen)
		variable pi_pcc_dn_hat(ne,npart,nscen)
		
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
		p_pcc_DA(e) == sum(p_pcc_DA_hat(e,:));
			
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
		npart2 = npart_in;
		variable fdn_part(ne,npart2) binary
		variable rho_DA_plus_part(ne,npart2) binary
		variable rho_RT_minus_part(ne,npart2,nscen) binary
		variable fup_part(ne,npart2,nscen) binary
		
		
		sum(fdn_part(e,:)) == 1;
		sum(rho_DA_plus_part(e,:)) == 1;
		for s = 1:nscen
			sum(rho_RT_minus_part(e,:,s)) == 1;
			sum(fup_part(e,:,s)) == 1;
		end
		
		
		fdn_lo = -10;
		fdn_up = 5;
		
		rhominus_lo = 0;
		rhominus_up = 30;
		
		fup_lo = -5;
		fup_up = 10;
		
		rhoplus_lo = 0;
		rhoplus_up = 30;
		
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
		variable rho_DA_minus_hat(ne,npart2)
		variable fup_hat(ne,npart2)
		variable fdn_hat(ne,npart2,nscen)
		variable rho_RT_plus_hat(ne,npart2,nscen)
		
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
	
	%% day-ahead primal constraints
	sum(p_gen_DA) + wind_DA + shed_p_DA + p_pcc_DA == sum(p_dem_DA);


	Pmin_gen(gens_area) <= p_gen_DA <= Pmax_gen(gens_area);
	Pmin_dem(dems_area) <= p_dem_DA <= Pmax_dem(dems_area);
	0 <= wind_DA <= sum(Wmax_mean_DA(windg_area));
	0 <= shed_p_DA <= sum(p_dem_DA);
	
	fe_dn <= p_pcc_DA <= fe_up;
	
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

	for k = 1:nl
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
		Vmin(buses).^2 <= v_sq(:,s) <= Vmax(buses).^2;

		sum(q_flow(:,:,s),2) == q_inj_rt(:,s);

		q_inj_rt(:,s) == incidentG' * q_gen_RT(:,s) ...
					- incidentD' * q_dem_RT(:,s) + shed_q_RT(:,s) + incident_pcc * q_pcc_rt(s);


		for k = 1:nl
			if SlmMax(branch_num(k)) ~= 0
				p_flow(fbus(k),tbus(k),s).^2 + q_flow(fbus(k),tbus(k),s).^2 <= SlmMax(branch_num(k)).^2;
			end
		end
	end
	
	p_flow_matrix = reshape(sum(p_flow,2),size(p_inj_rt));
	p_flow_matrix ==  p_inj_rt;
	for s = 1:nscen
		p_pcc_rt(e,s) == p_pcc_DA + pup_pcc(e,s) - pdn_pcc(e,s);
		p_gen_RT(:,s) == p_gen_DA + pup_g(:,s) - pdn_g(:,s);
		p_dem_RT(:,s) == p_dem_DA - pup_d(:,s) + pdn_d(:,s);

		p_inj_rt(:,s) == incidentG' * p_gen_RT(:,s) - incidentD' * p_dem_RT(:,s)...
				+ incidentW' * wind_RT(:,s) + shed_p_RT(:,s) + incident_pcc * p_pcc_rt(e,s);
			
		0 <= shed_p_RT(:,s) <= incidentD' * p_dem_RT(:,s);
	end
	

	for s = 1:nscen
		Pmin_gen(gens_area) <= p_gen_RT(:,s) <= Pmax_gen(gens_area);
		Qmin_gen(gens_area) <= q_gen_RT(:,s) <= Qmax_gen(gens_area);

		Pmin_dem(dems_area) <= p_dem_RT(:,s) <= Pmax_dem(dems_area);
		Qmin_dem(dems_area) <= q_dem_RT(:,s) <= Qmax_dem(dems_area);

		Wmin(windg_area,s) <= wind_RT(:,s) <= Wmax(windg_area,s);
	end
	
	fe_dn <= p_pcc_rt(e,:) <= fe_up;
% 	0 <= shed_p_RT <= incidentD' * p_dem_RT;
    for s = 1:nscen
		pup_pcc(e,s) <= VOLL * p_pcc_comp(e,s);
		pdn_pcc(e,s) <= VOLL * (1 - p_pcc_comp(e,s));
		
		pup_g(:,s) <= VOLL * p_gen_comp(:,s);
		pdn_g(:,s) <= VOLL * (1 - p_gen_comp(:,s));
		
		pup_d(:,s) <= VOLL * p_dem_comp(:,s);
		pdn_d(:,s) <= VOLL * (1 - p_dem_comp(:,s));
	end
	%% KKT stationary conditions
	for g = 1:ng
		offer_gen_DA(gens_area(g),2) - lambda_DA_e - varsigma_gen_DA_minus_e(g) + varsigma_gen_DA_plus_e(g) ...
				- sum(prob_wscen(1,:)*(offer_gen_DA(gens_area(g),2))) + sum(zeta_p_gen_RT(g,:))  == 0; % p_gen_DA
		for s = 1:nscen
			prob_wscen(1,s) * offer_gen_upreg(gens_area(g)) + zeta_p_gen_RT(g,s) - epsilon_up_RT(g,s) == 0; % pup_g
			prob_wscen(1,s) * offer_gen_downreg(gens_area(g)) - zeta_p_gen_RT(g,s) - epsilon_dn_RT(g,s) == 0; % pdn_g
			
			prob_wscen(1,s)*offer_gen_DA(gens_area(g),2) - zeta_p_gen_RT(g,s) - varsigma_gen_RT_minus(g,s) ...
				+ varsigma_gen_RT_plus(g,s) - incidentG(g,:) *lambda_p_RT(:,s) == 0; % p_gen_RT
			
			- kappa_gen_RT_minus(g,s) ...
				+ kappa_gen_RT_plus(g,s) - incidentG(g,:) *lambda_q_RT(:,s) == 0; % q_gen_RT
		end
	end
	
	for d = 1:nd
		-bid_dem_DA(dems_area(d)) + lambda_DA_e - varsigma_dem_DA_minus_e(d) + varsigma_dem_DA_plus_e(d) ...
				+  sum(prob_wscen(1,:)*bid_dem_DA(dems_area(d)) ) + sum(zeta_p_dem_RT(d,:)) - Upsilon_DA_plus == 0; % p_dem_DA
		for s = 1:nscen
			prob_wscen(1,s) * bid_dem_upreg(dems_area(d)) - zeta_p_dem_RT(d,s) - varepsilon_up_RT(d,s) == 0; % pup_d
			prob_wscen(1,s) * bid_dem_downreg(dems_area(d)) + zeta_p_dem_RT(d,s) - varepsilon_dn_RT(d,s) == 0; % pdn_d
			
			-prob_wscen(1,s)*bid_dem_DA(dems_area(d)) - zeta_p_dem_RT(d,s) - varsigma_dem_RT_minus(d,s) ...
				+ varsigma_dem_RT_plus(d,s) + incidentD(d,:) * (lambda_p_RT(:,s) - Upsilon_RT_plus(:,s)) == 0; % p_dem_RT
			
			- kappa_dem_RT_minus(d,s) + kappa_dem_RT_plus(d,s) ...
				+ incidentD(d,:) *lambda_q_RT(:,s) - incidentD(d,:)*Upsilon_Q_RT_plus(:,s) == 0; % q_dem_RT
		end
	end
	
	for n = 1:nb
		for s = 1:nscen
			VOLL - lambda_p_RT(n,s) - Upsilon_RT_minus(n,s) + Upsilon_RT_plus(n,s) == 0; % shed_p_RT
			VOLL - lambda_q_RT(n,s) - Upsilon_Q_RT_minus(n,s) + Upsilon_Q_RT_plus(n,s) == 0;
		end
	end
	
	for w = 1:nw
		for s = 1:nscen
			offer_wind - incidentW(w,:) * lambda_p_RT(:,s) - nu_RT_minus(w,s) + nu_RT_plus(w,s) == 0; % wind_RT
		end
	end
	
	
	for l = 1:nl
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
			
			
% 			beta_RT(fbus(l),tbus(l),s) == - beta_RT(tbus(l),fbus(l),s);
		end
	end
% 	gamma_RT_w <= 2;
% % 	gamma_RT_l <= 1;
% 	eta_RT_p == 0;
% 	eta_RT_q == 0;
% 	sigma_RT_plus(2:3,:) == 0;
% 	sigma_RT_minus == 0;
% eta_RT_s == 0;
	for s = 1:nscen
		for k = 1:length(n1)
			gamma_RT_l(n1(k),n2(k),s) == 0;
			gamma_RT_w(n1(k),n2(k),s) == 0;
			beta_RT(n1(k),n2(k),s) == 0;
			eta_RT_s(n1(k),n2(k),s) == 0;
		end
		
		prob_wscen(1,s)*pi_pcc_DA - incident_pcc' * lambda_p_RT(:,s) - zeta_pcc_RT(s) - rho_RT_minus(e,s)...
			+ rho_RT_plus(e,s) == 0; % p_pcc_RT
		
		incident_pcc' * lambda_q_RT(:,s) == 0; % q_pcc_RT
		
		prob_wscen(1,s)*pi_pcc_up + zeta_pcc_RT(s) - epsilon_pcc_up_RT(s) == 0; % pup_pcc
		prob_wscen(1,s)*pi_pcc_dn - zeta_pcc_RT(s) - epsilon_pcc_dn_RT(s) == 0; % pdn_pcc
	end
	pi_pcc_DA - sum(prob_wscen(1,:)*pi_pcc_DA) - lambda_DA_e - rho_DA_minus...
		+ rho_DA_plus + sum(zeta_pcc_RT) == 0; % p_pcc_DA
	
	VOLL - lambda_DA_e - Upsilon_DA_minus + Upsilon_DA_plus == 0; % shed_p_DA
	offer_wind - lambda_DA_e - iota_minus_DA + iota_plus_DA == 0; % wind_DA
	
	
	%% KKT complimentarity constraints
	
% 	for g = 1:ng
% 		varsigma_gen_DA_minus_e(g) <= VOLL * bin_gen(g,1);
% 		p_gen_DA(g) - Pmin_gen(gens_area(g)) <= VOLL *  (1 - bin_gen(g,1));
% 
% 		varsigma_gen_DA_plus_e(g) <= VOLL * bin_gen(g,2);
% 		Pmax_gen(gens_area(g)) - p_gen_DA(g) <= VOLL *  (1 - bin_gen(g,2));
% 		
% 	end
% 	
% 	
% 	for d = 1:nd
% 		varsigma_dem_DA_minus_e(d) <= VOLL * bin_dem(d,1);
% 		p_dem_DA(d) - Pmin_dem(dems_area(d)) <= VOLL *  (1 - bin_dem(d,1));
% 
% 		varsigma_dem_DA_plus_e(d) <= VOLL * bin_dem(d,2);
% 		Pmax_dem(dems_area(d)) - p_dem_DA(d) <= VOLL *  (1 - bin_dem(d,2));
% 	end
%  	
% 	iota_minus_DA <= VOLL * bin_windlim_DA(1);
% 	wind_DA <= VOLL * (1 - bin_windlim_DA(1));
% 	
% 	iota_plus_DA <= VOLL * bin_windlim_DA(2);
% 	sum(Wmax_mean_DA(windg_area)) - wind_DA <= VOLL * (1 - bin_windlim_DA(2));
% 	
% 	rho_DA_minus <= VOLL * bin_pcclim_DA(1);
% 	p_pcc_DA - fe_dn <= VOLL * (1 - bin_pcclim_DA(1));
% 	rho_DA_plus <= VOLL * bin_pcclim_DA(2);
% 	fe_up - p_pcc_DA <= VOLL * (1 - bin_pcclim_DA(2));
% 	
% % 	for l = 1:nl
% % 		for s = 1:nscen
% % 			
% % % 			i_sq(fbus(l),tbus(l),s) * v_sq(fbus(l),s) - (p_flow(fbus(l),tbus(l),s)^2 + q_flow(fbus(l),tbus(l),s)^2) ...
% % % 				<= VOLL * (1 - bin_socpflow(l,s));
% % 			gamma_RT_p(fbus(l),tbus(l),s) <= VOLL * bin_socpflow(1,fbus(l),tbus(l),s);
% % 			square(p_flow(fbus(l),tbus(l),s))  <= VOLL * (1 - bin_socpflow(1,fbus(l),tbus(l),s));
% % 			
% % 			gamma_RT_q(fbus(l),tbus(l),s) <= VOLL * bin_socpflow(2,fbus(l),tbus(l),s);
% % 			square(q_flow(fbus(l),tbus(l),s)) <= VOLL * (1 - bin_socpflow(2,fbus(l),tbus(l),s));
% % 			
% % 			gamma_RT_l(fbus(l),tbus(l),s) <= VOLL * bin_socpflow(3,fbus(l),tbus(l),s);
% % 			square(i_sq(fbus(l),tbus(l),s) - v_sq(fbus(l),s)) <= VOLL * (1 - bin_socpflow(3,fbus(l),tbus(l),s));
% % 			
% % 			gamma_RT_w(fbus(l),tbus(l),s) <= VOLL * bin_socpflow(4,fbus(l),tbus(l),s);
% % 			square(i_sq(fbus(l),tbus(l),s) + v_sq(fbus(l),s)) <= VOLL * (1 - bin_socpflow(4,fbus(l),tbus(l),s));
% % 			
% % 			
% % 			
% % % 			gamma_RT_p(tbus(l),fbus(l),s) <= VOLL * bin_socpflow(1,tbus(l),fbus(l),s);
% % % 			square(p_flow(tbus(l),fbus(l),s)) <= VOLL * (1 - bin_socpflow(1,tbus(l),fbus(l),s));
% % % 			
% % % 			gamma_RT_q(tbus(l),fbus(l),s) <= VOLL * bin_socpflow(2,tbus(l),fbus(l),s);
% % % 			square(q_flow(tbus(l),fbus(l),s)) <= VOLL * (1 - bin_socpflow(2,tbus(l),fbus(l),s));
% % % 			
% % % 			gamma_RT_l(tbus(l),fbus(l),s) <= VOLL * bin_socpflow(3,tbus(l),fbus(l),s);
% % % 			square(i_sq(tbus(l),fbus(l),s) - v_sq(tbus(l),s)) <= VOLL * (1 - bin_socpflow(3,tbus(l),fbus(l),s));
% % % 			
% % % 			gamma_RT_w(tbus(l),fbus(l),s) <= VOLL * bin_socpflow(4,tbus(l),fbus(l),s);
% % % 			square(i_sq(tbus(l),fbus(l),s) + v_sq(tbus(l),s)) <= VOLL * (1 - bin_socpflow(4,tbus(l),fbus(l),s));
% % 
% % 
% % % 			if SlmMax(k) ~= 0
% % % 				eta_RT(l,s) <= VOLL *bin_flowlim(l,s);
% % % 				SlmMax(l)^2 - p_flow(fbus(l),tbus(l),s)^2 - q_flow(fbus(l),tbus(l),s)^2 <= VOLL * (1 - bin_flowlim(l,s));
% % % 			else
% % % 				eta_RT(l,s) <= VOLL *bin_flowlim(l,s);
% % % 				VOLL - p_flow(fbus(l),tbus(l),s)^2 - q_flow(fbus(l),tbus(l),s)^2 <= VOLL * (1 - bin_flowlim(l,s));
% % % 			end
% % 
% % % 			eta_RT_p(fbus(l),tbus(l),s) <= VOLL * bin_flowlim(1,fbus(l),tbus(l),s);
% % % 			(p_flow(fbus(l),tbus(l),s))^2 <= VOLL * (1 - bin_flowlim(1,fbus(l),tbus(l),s));
% % % 			
% % % 			eta_RT_q(fbus(l),tbus(l),s) <= VOLL * bin_flowlim(2,fbus(l),tbus(l),s);
% % % 			(q_flow(fbus(l),tbus(l),s))^2 <= VOLL * (1 - bin_flowlim(2,fbus(l),tbus(l),s));
% % % 			
% % % 			eta_RT_p(tbus(l),fbus(l),s) <= VOLL * bin_flowlim(1,tbus(l),fbus(l),s);
% % % 			(p_flow(tbus(l),fbus(l),s))^2 <= VOLL * (1 - bin_flowlim(1,tbus(l),fbus(l),s));
% % % 			
% % % 			eta_RT_q(tbus(l),fbus(l),s) <= VOLL * bin_flowlim(2,tbus(l),fbus(l),s);
% % % 			(q_flow(tbus(l),fbus(l),s))^2 <= VOLL * (1 - bin_flowlim(2,tbus(l),fbus(l),s));
% % 		end
% % 	end
% 	
% 	for n = 1:nb
% 		for s = 1:nscen
% 			sigma_RT_minus(n,s) <= VOLL * bin_voltlim(1,n,s);
% 			v_sq(n,s) - Vmin(buses(n)).^2 <= VOLL * (1 - bin_voltlim(1,n,s));
% 			
% 			sigma_RT_plus(n,s) <= VOLL * bin_voltlim(2,n,s);
% 			Vmax(buses(n)).^2 - v_sq(n,s) - 1e-6 <= VOLL * (1 - bin_voltlim(2,n,s));
% 		
% 		end
% 	end
% 	
% 	for w = 1:nw
% 		for s = 1:nscen
% 			nu_RT_minus(w,s) <= VOLL * bin_windlim_RT(1,w,s);
% 			wind_RT(w,s) <= VOLL * (1 - bin_windlim_RT(1,w,s));
% 			
% 			nu_RT_plus(w,s) <= VOLL * bin_windlim_RT(2,w,s);
% 			Wmax(windg_area(w),s) - wind_RT(w,s) <= VOLL * (1 - bin_windlim_RT(2,w,s));
% 		end
% 	end
% 	
% 	for g = 1:ng
% 		for s = 1:nscen
% 			varsigma_gen_RT_minus(g,s) <= VOLL * bin_pgenlim_RT(1,g,s);
% 			p_gen_RT(g,s) - Pmin_gen(gens_area(g)) <= VOLL * (1 - bin_pgenlim_RT(1,g,s));
% 			
% 			varsigma_gen_RT_plus(g,s) <= VOLL * bin_pgenlim_RT(2,g,s);
% 			Pmax_gen(gens_area(g)) - p_gen_RT(g,s) <= VOLL * (1 - bin_pgenlim_RT(2,g,s));
% 			
% 			kappa_gen_RT_minus(g,s) <= VOLL * bin_qgenlim_RT(1,g,s);
% 			q_gen_RT(g,s) - Qmin_gen(gens_area(g)) <= VOLL * (1 - bin_qgenlim_RT(1,g,s));
% 			
% 			kappa_gen_RT_plus(g,s) <= VOLL * bin_qgenlim_RT(2,g,s);
% 			Qmax_gen(gens_area(g)) - q_gen_RT(g,s) <= VOLL * (1 - bin_qgenlim_RT(2,g,s));
% 		end
% 	end
% 	
% 	for d = 1:nd
% 		for s = 1:nscen
% 			varsigma_dem_RT_minus(d,s) <= VOLL * bin_pdemlim_RT(1,d,s);
% 			p_dem_RT(d,s) - Pmin_dem(dems_area(d)) <= VOLL * (1 - bin_pdemlim_RT(1,d,s));
% 			
% 			varsigma_dem_RT_plus(d,s) <= VOLL * bin_pdemlim_RT(2,d,s);
% 			Pmax_dem(dems_area(d)) - p_dem_RT(d,s) <= VOLL * (1 - bin_pdemlim_RT(2,d,s));
% 			
% 			kappa_dem_RT_minus(d,s) <= VOLL * bin_qdemlim_RT(1,d,s);
% 			q_dem_RT(d,s) - Qmin_dem(dems_area(d)) <= VOLL * (1 - bin_qdemlim_RT(1,d,s));
% 			
% 			kappa_dem_RT_plus(d,s) <= VOLL * bin_qdemlim_RT(2,d,s);
% 			Qmax_dem(dems_area(d)) - q_dem_RT(d,s) <= VOLL * (1 - bin_qdemlim_RT(2,d,s));
% 		end
% 	end
% 	
% 	for s = 1:nscen
% 		rho_RT_minus(s) - 1e-6 <= VOLL * bin_pcclim_RT(1,s);
% 		p_pcc_rt(s) - fe_dn <= VOLL *(1 - bin_pcclim_RT(1,s));
% 		
% 		rho_RT_plus(s) - 1e-7 <= VOLL * bin_pcclim_RT(2,s);
% 		fe_up - p_pcc_rt(s) <= VOLL *(1 - bin_pcclim_RT(2,s));
% 	end
% 	
% 	for g = 1:ng
% 		for s = 1:nscen
% 			epsilon_up_RT(g,s) <= VOLL * bin_gen_RT(1,g,s);
% 			pup_g(g,s) <= VOLL * (1 - bin_gen_RT(1,g,s));
% 			
% 			epsilon_dn_RT(g,s) <= VOLL * bin_gen_RT(2,g,s);
% 			pdn_g(g,s) <= VOLL * (1 - bin_gen_RT(2,g,s));
% 			
% 		end
% 	end
% 	
% 	for d = 1:nd
% 		for s = 1:nscen
% 			varepsilon_up_RT(d,s) <= VOLL * bin_dem_RT(1,d,s);
% 			pup_d(d,s) - 5e-8 <= VOLL * (1 - bin_dem_RT(1,d,s));
% 			
% 			varepsilon_dn_RT(d,s) <= VOLL * bin_dem_RT(2,d,s);
% 			pdn_d(d,s) - 4e-8 <= VOLL * (1 - bin_dem_RT(2,d,s));
% 			
% 		end
% 	end
% 	
% 	for s = 1:nscen
% 		epsilon_pcc_up_RT(s) <= VOLL * bin_pcc_RT(1,s);
% 		pup_pcc(s) <= VOLL * (1 - bin_pcc_RT(1,s));
% 
% 		epsilon_pcc_dn_RT(s) <= VOLL * bin_pcc_RT(2,s);
% 		pdn_pcc(s) <= VOLL * (1 - bin_pcc_RT(2,s));
% 	end
% 	
% 	for n = 1:nb
% 		for s = 1:nscen
% 			Upsilon_RT_minus(n,s) <= 1.5*VOLL *bin_shed_RT(1,n,s);
% 			shed_p_RT(n,s) <= VOLL * (1 - bin_shed_RT(1,n,s));
% 			
% 			Upsilon_RT_plus(n,s) <= 1.5*VOLL *bin_shed_RT(2,n,s);
% 			incidentD(:,n)' * p_dem_RT(:,s) - shed_p_RT(n,s) <= VOLL * (1 - bin_shed_RT(2,n,s));
% 		end
% 	end
% 	
% 	for n = 1:nb
% 		for s = 1:nscen
% 			Upsilon_Q_RT_minus(n,s) <= 1.5*VOLL *bin_shed_Q_RT(1,n,s);
% 			shed_q_RT(n,s) <= VOLL * (1 - bin_shed_Q_RT(1,n,s));
% 			
% 			Upsilon_Q_RT_plus(n,s) <= 1.5*VOLL *bin_shed_Q_RT(2,n,s);
% 			incidentD(:,n)' * q_dem_RT(:,s) - shed_q_RT(n,s) <= VOLL * (1 - bin_shed_Q_RT(2,n,s));
% 		end
% 	end
% 	
% 	
% 	Upsilon_DA_minus <= 1.5 * VOLL * bin_shed_DA(1);
% 	shed_p_DA <= VOLL * (1 - bin_shed_DA(1));
% 	
% 	Upsilon_DA_plus <= 1.5 * VOLL * bin_shed_DA(2);
% 	sum(p_dem_DA) - shed_p_DA <= VOLL * (1 - bin_shed_DA(2));

	%% dual constraints
		for s = 1:nscen
			for l = 1:nl
				
				ui{s}(:,l) = A'*[gamma_RT_p(fbus(l),tbus(l),s); gamma_RT_q(fbus(l),tbus(l),s); gamma_RT_l(fbus(l),tbus(l),s)];
				bi{s}(:,l) = b' * gamma_RT_w(fbus(l),tbus(l),s);
				
% 				sum( ui{s}(:,l) - bi{s}(:,l) ) == 0;
			end
        end
        
       %%
    variable rho_DA_comp(ne) binary

    rho_DA_minus(e) <= VOLL* rho_DA_comp(e);
    rho_DA_plus(e) <= VOLL * (1 - rho_DA_comp(e));

    variable rho_RT_comp(ne,nscen) binary

    rho_RT_minus(e,:) <= VOLL* rho_RT_comp(e,:);
    rho_RT_plus(e,:) <= VOLL * (1 - rho_RT_comp(e,:));
	%% strong duality objective
	expression dual_obj_linecap(nscen)
    if strcmp(relax,'none')
        dual_objective_DA = sum(Pmin_gen(gens_area)' * varsigma_gen_DA_minus_e - Pmax_gen(gens_area)' * varsigma_gen_DA_plus_e) ...
            + sum(Pmin_dem(dems_area)' * varsigma_dem_DA_minus_e - Pmax_dem(dems_area)' * varsigma_dem_DA_plus_e )...
            - iota_plus_DA * sum(Wmax_mean_DA(windg_area)) ...
            + fe_dn * rho_DA_minus - fe_up * rho_DA_plus;
    elseif strcmp(relax,'McCormick')
        dual_objective_DA = sum(Pmin_gen(gens_area)' * varsigma_gen_DA_minus_e - Pmax_gen(gens_area)' * varsigma_gen_DA_plus_e) ...
            + sum(Pmin_dem(dems_area)' * varsigma_dem_DA_minus_e - Pmax_dem(dems_area)' * varsigma_dem_DA_plus_e )...
            - iota_plus_DA * sum(Wmax_mean_DA(windg_area)) ...
            + w_fdn_rhominus_DA(e) - w_fup_rhoplus_DA(e);
    end

	for s = 1:nscen
        if strcmp(relax,'none')
            dual_objective_RT(s) =  sum(Vmin(buses)'.^2 * sigma_RT_minus(:,s))  - sum(Vmax(buses)'.^2 * sigma_RT_plus(:,s)) ...
                - sum(nu_RT_plus(:,s)' *  Wmax(windg_area,s) )...
                + sum(Pmin_gen(gens_area)' * varsigma_gen_RT_minus(:,s) - Pmax_gen(gens_area)' * varsigma_gen_RT_plus(:,s)) ...
                + sum(Pmin_dem(dems_area)' * varsigma_dem_RT_minus(:,s) - Pmax_dem(dems_area)' * varsigma_dem_RT_plus(:,s)) ...
                + sum(Qmin_gen(gens_area)' * kappa_gen_RT_minus(:,s) - Qmax_gen(gens_area)' * kappa_gen_RT_plus(:,s)) ...
                + sum(Qmin_dem(dems_area)' * kappa_dem_RT_minus(:,s) - Qmax_dem(dems_area)' * kappa_dem_RT_plus(:,s)) ...
                + fe_dn * rho_RT_minus(s) - fe_up * rho_RT_plus(s);
        elseif strcmp(relax,'McCormick')
            dual_objective_RT(s) =  sum(Vmin(buses)'.^2 * sigma_RT_minus(:,s))  - sum(Vmax(buses)'.^2 * sigma_RT_plus(:,s)) ...
                - sum(nu_RT_plus(:,s)' *  Wmax(windg_area,s) )...
                + sum(Pmin_gen(gens_area)' * varsigma_gen_RT_minus(:,s) - Pmax_gen(gens_area)' * varsigma_gen_RT_plus(:,s)) ...
                + sum(Pmin_dem(dems_area)' * varsigma_dem_RT_minus(:,s) - Pmax_dem(dems_area)' * varsigma_dem_RT_plus(:,s)) ...
                + sum(Qmin_gen(gens_area)' * kappa_gen_RT_minus(:,s) - Qmax_gen(gens_area)' * kappa_gen_RT_plus(:,s)) ...
                + sum(Qmin_dem(dems_area)' * kappa_dem_RT_minus(:,s) - Qmax_dem(dems_area)' * kappa_dem_RT_plus(:,s)) ...
                + w_fdn_rhominus_RT(e,s) - w_fup_rhoplus_RT(e,s);
        end
		
		for l = 1:nl
			if SlmMax(branch_num(l)) ~= 0
				dual_obj_linecap(s) = dual_obj_linecap(s) - sum(  (eta_RT_s(fbus(l),tbus(l),s)+eta_RT_s(tbus(l),fbus(l),s)) * SlmMax(branch_num(l)) ) ;
			else
				eta_RT_s(fbus(l),tbus(l),s) + eta_RT_s(tbus(l),fbus(l),s) == 0;
			end
		end
	end
	dual_objective = dual_objective_DA + sum(dual_objective_RT) + sum(dual_obj_linecap);
	
	cost == dual_objective;
	
cvx_end
	
% cost
	
% if any(shed_p_RT > 10e-5)
% 	warning('There is load shedding in the DSO look-ahead Real Time calculation');
% end
% if any(shed_p_DA > 10e-5)
% 	warning('There is load shedding in the DSO look-ahead day-ahead calculation');
% end
	
result_lookahead.cost_RT = cost_RT;
result_lookahead.p_dem_RT = p_dem_RT;
result_lookahead.p_gen_RT = p_gen_RT;
result_lookahead.p_gen_DA = p_gen_DA;
result_lookahead.p_dem_DA = p_dem_DA;
result_lookahead.wind_DA = wind_DA;
result_lookahead.wind_RT = wind_RT;
result_lookahead.p_pcc_DA = p_pcc_DA;
result_lookahead.p_pcc_RT = p_pcc_rt;
result_lookahead.v_sq = v_sq;
result_lookahead.i_sq = i_sq;

end