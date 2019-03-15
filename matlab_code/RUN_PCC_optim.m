clc
clear all
close all


%% in case the script is submitted to a cluster
% configCluster

% c = parcluster('dcc');
% c.AdditionalProperties.EmailAddress = ;
% c.AdditionalProperties.MemUsage = '4GB';
% c.AdditionalProperties.WallTime = '05:00';
% c.AdditionalProperties.ProcsPerNode = 20;
% c.NumWorkers = 60;
% c.AdditionalProperties.QueueName = 'elektro';
% % clust.AdditionalProperties.Queue = 'elektro';
% c.saveProfile
% nw_avail = str2num(getenv('LSB_DJOB_NUMPROC'));

% disp(['NW Available: ' num2str(nw_avail)])

% pool=parpool(c,60);

% c.AdditionalProperties
%% Setup CVX in case this has not been installed yet

%cd cvx
% cvx_setup

%cd ..
%% setup the necessary variables and the case data

test_case = 'case_24TSO_3DSO'; %name of the file with case data
mpc = eval(test_case);
warning('off','all')
simpoints = 24; % how many different wind penetration points do you want to simulate
pcc_points = 1;
wind_factor = linspace(0.2,4.6,simpoints);
PCC_lim_factor = 1; %linspace(0.1,2.5,pcc_points);
line_factor = 1; %linspace(0,2,35);
line_var_cap = 35; % line ID number to vary

	[nb,busL,busN,busType,Pd,Qd,Vmax,Vmin,statusG,activeG,activeD,...
	ng,nd,genL,genN,incidentG,demL,incidentD,incidentW,Qmax_gen,Qmin_gen,...
    Pmax_gen,Pmin_gen,Qmax_dem,Qmin_dem,Pmax_dem,Pmin_dem,statusL,...
	activeL,nl,fbusL,tbusL,SlmMax,fbusN,tbusN,incidentF,incidentT,...
    Yf,Yt,YfP,YtP,Ybus,edges,offer_gen_DA,offer_gen_upreg,...
	offer_gen_downreg,bid_dem_DA,bid_dem_upreg,bid_dem_downreg,...
	nr_vals_cost,p,cost_type,ys,areas,bus_areacodes,....
	incident_singlephase,Y_ft,rBranch,xBranch,narea,prob_wscen,Wmax,Wmin,...
	n_wgen,nscen,Wmax_mean_DA,windgL,...
	offer_wind_DA,offer_wind_up,offer_wind_dn] = Data_Reader(mpc,double.empty(1,0));

[area_codes, n_areas, DSO_codes, fbusL_area, tbusL_area, branch_num_area, ...
	fbusL_ext_area, tbusL_ext_area, branch_num_ext_area, overlap,...
	overlap_t, neigh_area, ext_area, ext_area_singlephase,fbus_local,tbus_local,incidence_area,PCC_branch_id] ...
	= find_overlap(bus_areacodes, fbusL, tbusL, areas);

TSO_lines = sub2ind([nb nb],fbusL_area{1},tbusL_area{1});

[PCC_lim_factor_grided,wind_factor_grided,line_factor_gridded] = ndgrid(PCC_lim_factor,wind_factor,line_factor);

base_PCC_capacity = mpc.branch(PCC_branch_id,6);
base_varLine_capacity = mpc.branch(line_var_cap,6);

variance_base =	[750	740	760	300	300	300	300];
mean_base = [200	200	200	40	40	40	10];

%%%%%%%%%%%%%%%%
nscen_IS=4; % set the number of in-sample scenarios
in_sample_iter = numel(PCC_lim_factor_grided);
[ii,jj,pp] = ind2sub(size(PCC_lim_factor_grided),1:numel(PCC_lim_factor_grided));

nscen_OOS = 200;
outofsample_iter = numel(PCC_lim_factor_grided);

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

%% Run the in sample optimization

% 
t_sim_IS = tic;
for kk = 1:in_sample_iter
	[results_PCCoptim(kk), results_FC(kk), results_CM(kk)] ...
		= run_insample_optim(mpc,base_PCC_capacity,PCC_lim_factor_grided,...
			wind_factor_grided,line_factor_gridded,line_var_cap,variance_base,mean_base,PCC_branch_id,nscen_IS,kk);
end
t_fin_IS = toc(t_sim_IS);
disp(['Time to finish IN-sample optimization: ' num2str(t_fin_IS)]);
% save('./solutions/results_inceasing_windpenetration.mat');
% load('./solutions/results_inceasing_windpenetration_run7.mat')
%% Run the out-of-sample validation

for k = 1:in_sample_iter
	DA_outcome_PCC(k) = results_PCCoptim(k).DA_outcome;
	DA_outcome_CM(k) = results_CM(k).DA_outcome; 
end
t_sim_OOS = tic;
parfor kk = 1:outofsample_iter
	[results_PCCoptim_OOS(kk), results_FC_OOS(kk), results_CM_OOS(kk)] ...
		= run_outofsample_optim(mpc,base_PCC_capacity,PCC_lim_factor_grided,...
			wind_factor_grided,line_factor_gridded,line_var_cap,variance_base,mean_base,PCC_branch_id,nscen_IS,nscen_OOS,...
			DA_outcome_PCC(kk),DA_outcome_CM(kk),results_FC(kk),kk);
end
t_fin_OOS = toc(t_sim_OOS);
disp(['Time to finish OUT-OF-sample validation: ' num2str(t_fin_OOS)]);

% save('./solutions/results_inceasing_windpenetration.mat');

%% Find the power flow incurred by the DA dispatches
% The RT simulation is done without redispatch and without line flow and
% voltage limits.
for k = 1:in_sample_iter
	DA_outcome_PCC(k) = results_PCCoptim(k).DA_outcome;
	DA_outcome_CM(k) = results_CM(k).DA_outcome; 
end
t_sim_OOS = tic;
parfor kk = 1:outofsample_iter
	[results_PCCoptim_DAflow(kk), results_FC_DAflow(kk), results_CM_DAflow(kk)] ...
		= run_DA_powerflow(mpc,base_PCC_capacity,PCC_lim_factor_grided,...
			wind_factor_grided,line_factor_gridded,line_var_cap,variance_base,mean_base,PCC_branch_id,nscen_IS,nscen_OOS,...
			DA_outcome_PCC(kk),DA_outcome_CM(kk),results_FC(kk),kk);
end
t_fin_OOS = toc(t_sim_OOS);
disp(['Time to finish OUT-OF-sample validation: ' num2str(t_fin_OOS)]);


save('./solutions/results_inceasing_windpenetration.mat');%,...
% 	'results_PCCoptim','results_FC','results_CM',...
% 	'results_PCCoptim_OOS','results_FC_OOS','results_CM_OOS',...
% 	'results_PCCoptim_DAflow','results_FC_DAflow','results_CM_DAflow');
%% extract the data to plot in-sample results

% load('.\solutions\results_inceasing_windpenetration_run6.mat')
% 

cost_RT_PCCoptim = zeros(size(wind_factor_grided));
cost_RT_conv = zeros(size(wind_factor_grided));
co_optim_cost =  zeros(size(wind_factor_grided));
pcc_optim_cost =  zeros(size(wind_factor_grided));
conv_cost =  zeros(size(wind_factor_grided));
cost_RT_cooptim =  zeros(size(wind_factor_grided));
DA_co_optim_cost =  zeros(size(wind_factor_grided));
DA_pcc_optim_cost =  zeros(size(wind_factor_grided));
DA_conv_cost =  zeros(size(wind_factor_grided));
wind_penetration_PCC =  zeros(size(wind_factor_grided));
wind_penetration_FC =  zeros(size(wind_factor_grided));
wind_penetration_CM =  zeros(size(wind_factor_grided));
wind_penetration_IS =  zeros(size(wind_factor_grided));


for kk = 1:in_sample_iter
	
	co_optim_cost(kk) = results_FC(kk).cost_total;
	pcc_optim_cost(kk) = results_PCCoptim(kk).DA_outcome.cost;
	conv_cost(kk) = results_CM(kk).DA_outcome.cost;
	
	
	cost_RT_cooptim(kk) = sum(results_FC(kk).cost_RT)/nscen_IS;
	for s = 1:nscen_IS
		cost_RT_PCCoptim(kk) = cost_RT_PCCoptim(kk)...
			+ (results_PCCoptim(kk).RT_outcome(s).cost_RT)/nscen_IS;
		cost_RT_conv(kk) = cost_RT_conv(kk)...
			+ (results_CM(kk).RT_outcome(s).cost_RT)/nscen_IS;
	end
	
	DA_co_optim_cost(kk) = results_FC(kk).cost_DA;
	DA_pcc_optim_cost(kk) = results_PCCoptim(kk).DA_outcome.cost_DA;
	DA_conv_cost(kk) = results_CM(kk).DA_outcome.cost_DA;


	for s = 1:nscen_IS
		wind_penetration_PCC(kk) = wind_penetration_PCC(kk) + results_PCCoptim(kk).RT_outcome(s).wind_penetration_total;
		wind_penetration_CM(kk) = wind_penetration_CM(kk) + results_CM(kk).RT_outcome(s).wind_penetration_total;
		
		wind_penetration_IS(kk) = wind_penetration_IS(kk) + results_PCCoptim(kk).RT_outcome(s).wind_penetration_offered;
	end
	wind_penetration_FC(kk) = wind_penetration_FC(kk) + sum(results_FC(kk).wind_penetration_total);

	wind_penetration_PCC(kk) = wind_penetration_PCC(kk)/nscen_IS;	
	wind_penetration_CM(kk) = wind_penetration_CM(kk)/nscen_IS;			
	wind_penetration_FC(kk) = wind_penetration_FC(kk)/nscen_IS;	
	wind_penetration_IS(kk) = wind_penetration_IS(kk)/nscen_IS;

	
end

%% extract the data to plot out-of-sample results

co_optim_cost_OOS =  zeros(size(wind_factor_grided));
pcc_optim_cost_OOS =  zeros(size(wind_factor_grided));
conv_cost_OOS =  zeros(size(wind_factor_grided));
wind_penetration_OOS =  zeros(size(wind_factor_grided));
FC_costSTD_OOS =  zeros(size(wind_factor_grided));
PCCO_costSTD_OOS =  zeros(size(wind_factor_grided));
CM_costSTD_OOS =  zeros(size(wind_factor_grided));

FC_cost_OOS_quant09 =  zeros(size(wind_factor_grided));
PCCO_cost_OOS_quant09 =  zeros(size(wind_factor_grided));
CM_cost_OOS_quant09 =  zeros(size(wind_factor_grided));
FC_cost_OOS_quant01 =  zeros(size(wind_factor_grided));
PCCO_cost_OOS_quant01 =  zeros(size(wind_factor_grided));
CM_cost_OOS_quant01 =  zeros(size(wind_factor_grided));

p_flow_PCC_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
p_flow_CM_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
p_flow_FC_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
q_flow_PCC_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
q_flow_CM_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
q_flow_FC_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
s_flow_PCC_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
s_flow_CM_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);
s_flow_FC_RT_expec =  zeros(length(PCC_branch_id),outofsample_iter);


PCC_subind = sub2ind([nb nb],fbusL(PCC_branch_id),tbusL(PCC_branch_id));
p_flow_PCC_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);
q_flow_PCC_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);
s_flow_PCC_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);


p_flow_FC_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);
q_flow_FC_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);
s_flow_FC_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);


p_flow_CM_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);
q_flow_CM_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);
s_flow_CM_DAflow_expec = zeros(length(PCC_branch_id),outofsample_iter);


PCC_TSO_congestion = zeros(size(wind_factor_grided));
PCC_TSO_line_cong = zeros([size(wind_factor_grided,1) size(wind_factor_grided,2) size(wind_factor_grided,3) length(TSO_lines) ]);
FC_TSO_congestion = zeros(size(wind_factor_grided));
FC_TSO_line_cong = zeros([size(wind_factor_grided,1) size(wind_factor_grided,2) size(wind_factor_grided,3) length(TSO_lines) ]);

CM_TSO_congestion = zeros(size(wind_factor_grided));
CM_TSO_line_cong = zeros([size(wind_factor_grided,1) size(wind_factor_grided,2) size(wind_factor_grided,3) length(TSO_lines) ]);
for kk = 1:outofsample_iter

	co_optim_cost_OOS(kk) = results_FC_OOS(kk).DA_outcome.cost;
	pcc_optim_cost_OOS(kk) = results_PCCoptim_OOS(kk).DA_outcome.cost;
	conv_cost_OOS(kk) = results_CM_OOS(kk).DA_outcome.cost;
	FC_costSTD_OOS(kk) = results_FC_OOS(kk).DA_outcome.cost_STD;
	PCCO_costSTD_OOS(kk) = results_PCCoptim_OOS(kk).DA_outcome.cost_STD;
	CM_costSTD_OOS(kk) = results_CM_OOS(kk).DA_outcome.cost_STD;
	
	FC_cost_OOS_quant09(kk) = results_FC_OOS(kk).DA_outcome.quantile08;
	PCCO_cost_OOS_quant09(kk) = results_PCCoptim_OOS(kk).DA_outcome.quantile08;
	CM_cost_OOS_quant09(kk) = results_CM_OOS(kk).DA_outcome.quantile08;
	FC_cost_OOS_quant01(kk) = results_FC_OOS(kk).DA_outcome.quantile02;
	PCCO_cost_OOS_quant01(kk) = results_PCCoptim_OOS(kk).DA_outcome.quantile02;
	CM_cost_OOS_quant01(kk) = results_CM_OOS(kk).DA_outcome.quantile02;
	
	
	for s = 1:nscen_OOS
		wind_penetration_OOS(kk) = wind_penetration_OOS(kk) + results_PCCoptim_OOS(kk).RT_outcome(s).wind_penetration_offered;
	end
	wind_penetration_OOS(kk) = wind_penetration_OOS(kk)/nscen_OOS;
% 	wind_penetration_OOS(kk) = sum(Wmax_mean_DA)/sum(Pmax_dem);
	
	for s = 1:nscen_OOS
		p_flow_PCCO_RT{kk,s} = results_PCCoptim_OOS(kk).RT_outcome(s).p_flow;
		p_flow_FC_RT{kk,s} = results_FC_OOS(kk).RT_outcome(s).p_flow;
		p_flow_CM_RT{kk,s} = results_CM_OOS(kk).RT_outcome(s).p_flow;

		q_flow_PCC_RT{kk,s} = results_PCCoptim_OOS(kk).RT_outcome(s).q_flow;
		q_flow_FC_RT{kk,s} = results_FC_OOS(kk).RT_outcome(s).q_flow;
		q_flow_CM_RT{kk,s} = results_CM_OOS(kk).RT_outcome(s).q_flow;
	end
	
	for s = 1:nscen_OOS
		p_flow_PCC_RT_expec(:,kk) = p_flow_PCC_RT_expec(:,kk) + p_flow_PCCO_RT{kk,s}(PCC_subind)/nscen_OOS;
		q_flow_PCC_RT_expec(:,kk) = q_flow_PCC_RT_expec(:,kk) + q_flow_PCC_RT{kk,s}(PCC_subind)/nscen_OOS; 

		p_flow_FC_RT_expec(:,kk) = p_flow_FC_RT_expec(:,kk) + p_flow_FC_RT{kk,s}(PCC_subind)/nscen_OOS;
		q_flow_FC_RT_expec(:,kk) = q_flow_FC_RT_expec(:,kk) + q_flow_FC_RT{kk,s}(PCC_subind)/nscen_OOS;

		p_flow_CM_RT_expec(:,kk) = p_flow_CM_RT_expec(:,kk) + p_flow_CM_RT{kk,s}(PCC_subind)/nscen_OOS;
		q_flow_CM_RT_expec(:,kk) = q_flow_CM_RT_expec(:,kk) + q_flow_CM_RT{kk,s}(PCC_subind)/nscen_OOS;
		
		pcc_flow_PCCO_RT(:,kk,s) = sqrt(p_flow_PCCO_RT{kk,s}(PCC_subind).^2 + q_flow_PCC_RT{kk,s}(PCC_subind).^2);
		pcc_flow_FC_RT(:,kk,s) = sqrt(p_flow_FC_RT{kk,s}(PCC_subind).^2 + q_flow_FC_RT{kk,s}(PCC_subind).^2);
		pcc_flow_CM_RT(:,kk,s) = sqrt(p_flow_CM_RT{kk,s}(PCC_subind).^2 + q_flow_CM_RT{kk,s}(PCC_subind).^2);
	end
	
	STD_pcc_flow_PCCO_RT(:,kk) = std(pcc_flow_PCCO_RT(:,kk,:),0,3);
	STD_pcc_flow_FC_RT(:,kk) = std(pcc_flow_FC_RT(:,kk,:),0,3);
	STD_pcc_flow_CM_RT(:,kk) = std(pcc_flow_CM_RT(:,kk,:),0,3);
	
	quant09_pcc_flow_PCCO_RT(:,kk,:) = quantile(pcc_flow_PCCO_RT(:,kk,:),0.9,3);
	quant09_pcc_flow_FC_RT(:,kk,:) = quantile(pcc_flow_FC_RT(:,kk,:),0.9,3);
	quant09_pcc_flow_CM_RT(:,kk,:) = quantile(pcc_flow_CM_RT(:,kk,:),0.9,3);
	
	quant01_pcc_flow_PCCO_RT(:,kk,:) = quantile(pcc_flow_PCCO_RT(:,kk,:),0.1,3);
	quant01_pcc_flow_FC_RT(:,kk,:) = quantile(pcc_flow_FC_RT(:,kk,:),0.1,3);
	quant01_pcc_flow_CM_RT(:,kk,:) = quantile(pcc_flow_CM_RT(:,kk,:),0.1,3);

	s_flow_PCC_RT_expec(:,kk) = sqrt( p_flow_PCC_RT_expec(:,kk).^2 + q_flow_PCC_RT_expec(:,kk).^2);
	s_flow_FC_RT_expec(:,kk) = sqrt( p_flow_FC_RT_expec(:,kk).^2 + q_flow_FC_RT_expec(:,kk).^2);
	s_flow_CM_RT_expec(:,kk) = sqrt( p_flow_CM_RT_expec(:,kk).^2 + q_flow_CM_RT_expec(:,kk).^2);
	%% extract the results fromt the fixed DA-dispatch power flow
	p_flow_PCC_DAflow{kk} = results_PCCoptim_DAflow(kk).RT_outcome.p_flow;
	p_flow_FC_DAflow{kk} = results_FC_DAflow(kk).RT_outcome.p_flow;
	p_flow_CM_DAflow{kk} = results_CM_DAflow(kk).RT_outcome.p_flow;

	q_flow_PCC_DAflow{kk} = results_PCCoptim_DAflow(kk).RT_outcome.q_flow;
	q_flow_FC_DAflow{kk} = results_FC_DAflow(kk).RT_outcome.q_flow;
	q_flow_CM_DAflow{kk} = results_CM_DAflow(kk).RT_outcome.q_flow;


	p_flow_PCC_DAflow_expec(:,kk) = p_flow_PCC_DAflow_expec(:,kk) + p_flow_PCC_DAflow{kk}(PCC_subind);
	q_flow_PCC_DAflow_expec(:,kk) = q_flow_PCC_DAflow_expec(:,kk) + q_flow_PCC_DAflow{kk}(PCC_subind); 

	p_flow_FC_DAflow_expec(:,kk) = p_flow_FC_DAflow_expec(:,kk) + p_flow_FC_DAflow{kk}(PCC_subind);
	q_flow_FC_DAflow_expec(:,kk) = q_flow_FC_DAflow_expec(:,kk) + q_flow_FC_DAflow{kk}(PCC_subind);

	p_flow_CM_DAflow_expec(:,kk) = p_flow_CM_DAflow_expec(:,kk) + p_flow_CM_DAflow{kk}(PCC_subind);
	q_flow_CM_DAflow_expec(:,kk) = q_flow_CM_DAflow_expec(:,kk) + q_flow_CM_DAflow{kk}(PCC_subind);
	
	
	s_flow_PCC_DAflow_expec(:,kk) = sqrt( p_flow_PCC_DAflow_expec(:,kk).^2 + q_flow_PCC_DAflow_expec(:,kk).^2);
	s_flow_FC_DAflow_expec(:,kk) = sqrt( p_flow_FC_DAflow_expec(:,kk).^2 + q_flow_FC_DAflow_expec(:,kk).^2);
	s_flow_CM_DAflow_expec(:,kk) = sqrt( p_flow_CM_DAflow_expec(:,kk).^2 + q_flow_CM_DAflow_expec(:,kk).^2);
	
	%% the congestion in the TSO network
	congestion_test_limit = 3;
	for s = 1:nscen_OOS
% 		try
		PCC_TSO_num_lines_cong(kk,s) = sum(sum(abs(results_PCCoptim_OOS(kk).RT_outcome(s).congestion(areas{1},areas{1})))) ;
		FC_TSO_num_lines_cong(kk,s) = sum(sum(abs(results_FC_OOS(kk).RT_outcome(s).congestion(areas{1},areas{1}))));
		CM_TSO_num_lines_cong(kk,s) =  sum(sum(abs(results_CM_OOS(kk).RT_outcome(s).congestion(areas{1},areas{1}))));

		
		PCC_TSO_congestion(kk) = PCC_TSO_congestion(kk) + (PCC_TSO_num_lines_cong(kk,s))/nscen_OOS;
		PCC_TSO_line_cong(ii(kk),jj(kk),pp(kk),:) = PCC_TSO_line_cong(ii(kk),jj(kk),pp(kk),:) + reshape( results_PCCoptim_OOS(kk).RT_outcome(s).congestion(TSO_lines) /nscen_OOS,[1 1 1 length(TSO_lines) ]) ;
		
		FC_TSO_congestion(kk) = FC_TSO_congestion(kk) + (FC_TSO_num_lines_cong(kk,s))/nscen_OOS;
		FC_TSO_line_cong(ii(kk),jj(kk),pp(kk),:) = FC_TSO_line_cong(ii(kk),jj(kk),pp(kk),:) + reshape( results_FC_OOS(kk).RT_outcome(s).congestion(TSO_lines) /nscen_OOS,[1 1 1 length(TSO_lines) ]) ;
		
		CM_TSO_congestion(kk) = CM_TSO_congestion(kk) + (CM_TSO_num_lines_cong(kk,s))/nscen_OOS;
		CM_TSO_line_cong(ii(kk),jj(kk),pp(kk),:) = CM_TSO_line_cong(ii(kk),jj(kk),pp(kk),:) + reshape( results_CM_OOS(kk).RT_outcome(s).congestion(TSO_lines) /nscen_OOS,[1 1 1 length(TSO_lines) ]) ;
% 		catch mes
% 			PCC_TSO_congestion(kk) = PCC_TSO_congestion(kk) + sum(sum(abs(results_PCCoptim_OOS(kk).RT_outcome(s).congestion(areas{1},areas{1}))));
% 			PCC_TSO_line_cong(ii(kk),jj(kk),:) = PCC_TSO_line_cong(ii(kk),jj(kk),:) + reshape( results_PCCoptim_OOS(kk).RT_outcome(s).congestion(TSO_lines) /nscen_OOS,[1 1 length(TSO_lines) ]) ;
% 		end
	end
% 	PCC_TSO_congestion(kk) = PCC_TSO_congestion(kk)/nscen_OOS;
	
	for s = 1:nscen_OOS
		p_flow_PCC{kk,s} = results_PCCoptim_OOS(kk).RT_outcome(s).p_flow;
		p_flow_FC{kk,s} = results_FC_OOS(kk).RT_outcome(s).p_flow;
		p_flow_CM{kk,s} = results_CM_OOS(kk).RT_outcome(s).p_flow;
	end

	
end

PCC_TSO_line_cong = squeeze(PCC_TSO_line_cong);
FC_TSO_line_cong = squeeze(FC_TSO_line_cong);
CM_TSO_line_cong = squeeze(CM_TSO_line_cong);

%% Start plotting the results
close all
for k = 1:pcc_points
IS_results_actualWP = figure;
plot(wind_penetration_FC(k,:),-co_optim_cost(k,:))
hold on
plot(wind_penetration_PCC(k,:),-pcc_optim_cost(k,:))
plot(wind_penetration_CM(k,:),-conv_cost(k,:))
xlabel('Actual Wind Penetration')
ylabel('Social Welfare')
legend('SW with full co-optimization','SW with PCC optimizer','SW with sequential market clearing','Location','best')
grid on
saveas(IS_results_actualWP,['./solutions/IS_results_actualWP' num2str(k)],'png')
% distri_fig

ISfig = figure;
plot(wind_penetration_IS(k,:),-co_optim_cost(k,:),'Marker','x')
hold on
plot(wind_penetration_IS(k,:),-pcc_optim_cost(k,:),'Marker','x')
plot(wind_penetration_IS(k,:),-conv_cost(k,:),'Marker','x')
title('In-Sample Results')
xlabel('Wind Penetration')
ylabel('Expected Social Welfare [$]')
legend('SW with full co-optimization','SW with PCC optimizer','SW with sequential market clearing','Location','best')
grid on
% xlim([min(wind_penetration_IS),max(wind_penetration_IS)])
saveas(ISfig,['./solutions/IS_results' num2str(k)],'png')

OOSfig = figure;
plot(wind_penetration_OOS(k,:),-co_optim_cost_OOS(k,:),'Marker','x')
hold on
plot(wind_penetration_OOS(k,:),-pcc_optim_cost_OOS(k,:),'Marker','x')
plot(wind_penetration_OOS(k,:),-conv_cost_OOS(k,:),'Marker','x')
title('Out-Of-Sample validation')
xlabel('Wind Penetration')
ylabel('Expected Social Welfare [$]')
legend('SW with full co-optimization','SW with PCC optimizer','SW with sequential market clearing','Location','best')
grid on
% xlim([min(wind_penetration_OOS),max(wind_penetration_OOS)])
saveas(OOSfig,['./solutions/OOS_results' num2str(k)],'png')

% figure
% plot(wind_penetration_IS(k,:),cost_RT_cooptim(k,:))
% hold on
% plot(wind_penetration_IS(k,:),cost_RT_PCCoptim(k,:))
% plot(wind_penetration_IS(k,:),cost_RT_conv(k,:))
% xlabel('Wind Penetration + Spilled Wind')
% ylabel('RT Re-dispatch cost')
% legend('SW with full co-optimization','SW with PCC optimizer','SW with sequential market clearing','Location','best')
% grid on

% figure
% plot(wind_penetration_IS,-DA_co_optim_cost)
% hold on
% plot(wind_penetration_IS,-DA_pcc_optim_cost)
% plot(wind_penetration_IS,-DA_conv_cost)
% xlabel('Wind Penetration + Spilled Wind')
% ylabel('DA SW')
% legend('SW with full co-optimization','SW with PCC optimizer','SW with sequential market clearing','Location','best')
% grid on
ind = 0;
clear line_names
for l = 1:length(TSO_lines)
	if sum(PCC_TSO_line_cong(:,l)) == 0
		PCC_TSO_line_cong(:,l) = NaN;
% 		line_names{l} = [];
	else
		ind = ind+1;
		line_names{ind} = ['Line' num2str(l)];
	end
end

PCC_TSO_line_cong_plot = PCC_TSO_line_cong(:,all(~isnan(PCC_TSO_line_cong)));

congestion_plot_PCC = figure;
plot(wind_penetration_IS(k,:),PCC_TSO_line_cong_plot,'Marker','x')
grid on
title('Congestion plot PCC')
xlabel('Wind Penetration')
ylabel('Congestion Probability in RT')
legend(line_names,'Location','best')
saveas(congestion_plot_PCC,['./solutions/PCCcongestion_plot' num2str(k)],'png')


ind = 0;
clear line_names
for l = 1:length(TSO_lines)
	if sum(FC_TSO_line_cong(:,l)) == 0
		FC_TSO_line_cong(:,l) = NaN;
% 		line_names{l} = [];
	else
		ind = ind+1;
		line_names{ind} = ['Line' num2str(l)];
	end
end

FC_TSO_line_cong_plot = FC_TSO_line_cong(:,all(~isnan(FC_TSO_line_cong)));

try
	congestion_plot_FC = figure;
	plot(wind_penetration_IS(k,:),FC_TSO_line_cong_plot,'Marker','x')
	grid on
	title('Congestion plot FC')
	xlabel('Wind Penetration')
	ylabel('Congestion Probability in RT')
	legend(line_names,'Location','best')
	saveas(congestion_plot_FC,['./solutions/FCcongestion_plot' num2str(k)],'png')
catch congestion_plot_FC_error
	disp('no data for FC line congestion plot')
	close(congestion_plot_FC)
end


ind = 0;
clear line_names
for l = 1:length(TSO_lines)
	if sum(CM_TSO_line_cong(:,l)) == 0
		CM_TSO_line_cong(:,l) = NaN;
% 		line_names{l} = [];
	else
		ind = ind+1;
		line_names{ind} = ['Line' num2str(l)];
	end
end

CM_TSO_line_cong_plot = CM_TSO_line_cong(:,all(~isnan(CM_TSO_line_cong)));
try
	congestion_plot_CM = figure;
	plot(wind_penetration_IS(k,:),CM_TSO_line_cong_plot,'Marker','x')
	grid on
	title('Congestion plot CM')
	xlabel('Wind Penetration')
	ylabel('Congestion Probability in RT')
	legend(line_names,'Location','best')
	saveas(congestion_plot_CM,['./solutions/CMcongestion_plot' num2str(k)],'png')
catch congestion_plot_CM_error
	disp('no data for CM line congestion plot')
	close(congestion_plot_CM)
end
% for k = 1:3
p_flow_mean = zeros(simpoints,4);
for p = 1:simpoints
	for s = 1:nscen_IS
		p_flow_mean(p,1) = p_flow_mean(p,1) + sum(p_flow_PCC{p,s}(3,24))/nscen_IS;
		p_flow_mean(p,2) = p_flow_mean(p,2) + sum((p_flow_PCC{p,s}(9,11)))/nscen_IS;
		p_flow_mean(p,3) = p_flow_mean(p,3) + sum((p_flow_PCC{p,s}(10,12)))/nscen_IS;
		p_flow_mean(p,4) = p_flow_mean(p,4) + sum((p_flow_PCC{p,s}(15,21)))/nscen_IS;
	end
end
tranformer_flow_plot = figure;
plot(wind_penetration_IS(k,:),p_flow_mean)
grid on
title('Line flow in transformers PCC')
xlabel('Wind Penetration')
ylabel('Expected Active Power Flow [MW]')
legend('From 3 to 24','From 9 to 11','From 10 to 12','From 15 to 21','Location','best')
saveas(tranformer_flow_plot,['./solutions/p_flow_transformers_PCC' num2str(k)],'png')
% end
p_flow_mean = zeros(simpoints,4);
for p = 1:simpoints
	for s = 1:nscen_IS
		p_flow_mean(p,1) = p_flow_mean(p,1) + sum(p_flow_CM{p,s}(3,24))/nscen_IS;
		p_flow_mean(p,2) = p_flow_mean(p,2) + sum((p_flow_CM{p,s}(9,11)))/nscen_IS;
		p_flow_mean(p,3) = p_flow_mean(p,3) + sum((p_flow_CM{p,s}(10,12)))/nscen_IS;
		p_flow_mean(p,4) = p_flow_mean(p,4) + sum((p_flow_CM{p,s}(15,21)))/nscen_IS;
	end
end
tranformer_flow_plot = figure;
plot(wind_penetration_IS(k,:),p_flow_mean)
grid on
title('Line flow in transformers CM')
xlabel('Wind Penetration')
ylabel('Expected Active Power Flow [MW]')
legend('From 3 to 24','From 9 to 11','From 10 to 12','From 15 to 21','Location','best')
saveas(tranformer_flow_plot,['./solutions/p_flow_transformers_CM' num2str(k)],'png')

p_flow_mean = zeros(simpoints,4);
for p = 1:simpoints
	for s = 1:nscen_IS
		p_flow_mean(p,1) = p_flow_mean(p,1) + sum(p_flow_FC{p,s}(3,24))/nscen_IS;
		p_flow_mean(p,2) = p_flow_mean(p,2) + sum((p_flow_FC{p,s}(9,11)))/nscen_IS;
		p_flow_mean(p,3) = p_flow_mean(p,3) + sum((p_flow_FC{p,s}(10,12)))/nscen_IS;
		p_flow_mean(p,4) = p_flow_mean(p,4) + sum((p_flow_FC{p,s}(15,21)))/nscen_IS;
	end
end
tranformer_flow_plot = figure;
plot(wind_penetration_IS(k,:),p_flow_mean)
grid on
title('Line flow in transformers FC')
xlabel('Wind Penetration')
ylabel('Expected Active Power Flow [MW]')
legend('From 3 to 24','From 9 to 11','From 10 to 12','From 15 to 21','Location','best')
saveas(tranformer_flow_plot,['./solutions/p_flow_transformers_FC' num2str(k)],'png')

end


%% plot the SW together with the standard deviation
clc
set(groot, 'defaultLegendInterpreter','latex');

z = wind_penetration_OOS(1,:)'*100;

% red_color = [ones(simpoints,1) zeros(simpoints,1) zeros(simpoints,1)];

max_congestion_level = max(max([PCC_TSO_congestion; FC_TSO_congestion; CM_TSO_congestion]));

PCC_TSO_congestion_plot = PCC_TSO_congestion./max(max_congestion_level);
FC_TSO_congestion_plot = FC_TSO_congestion./max(max_congestion_level);
CM_TSO_congestion_plot = CM_TSO_congestion./max(max_congestion_level);


% generate some data to plot
normal = wind_penetration_OOS(1,:)*100;
Parameter1_mean = -co_optim_cost_OOS(1,:);
Parameter1_plus_std = FC_cost_OOS_quant09;
Parameter1_minus_std = FC_cost_OOS_quant01;
Parameter2_mean = -pcc_optim_cost_OOS(1,:);
Parameter2_plus_std = PCCO_cost_OOS_quant09;
Parameter2_minus_std = PCCO_cost_OOS_quant01;
Parameter3_mean = -conv_cost_OOS(1,:);
Parameter3_plus_std = CM_cost_OOS_quant09;
Parameter3_minus_std = CM_cost_OOS_quant01;

% create figure
OOS_STD_fig = figure('Units', 'normalized'); %, 'outerposition', [0 0 1 1]);
colormap(OOS_STD_fig,'hot');
% CT=cbrewer('div', 'RdYlGn', 13);
% colormap(OOS_STD_fig,flipud(CT(2:end-2,:)));
hold on

% plot the two solid fillings without border
fi1 = fill([z',fliplr(z')], [Parameter1_plus_std,fliplr(Parameter1_minus_std)], 0.9*[1 1 1], 'EdgeColor','none');
fi2 = fill([z',fliplr(z')], [Parameter2_plus_std,fliplr(Parameter2_minus_std)], 0.5*[1 1 1], 'EdgeColor','none');
fi3 = fill([z',fliplr(z')], [Parameter3_plus_std,fliplr(Parameter3_minus_std)], 0.1*[1 1 1], 'EdgeColor','none');
set(fi1,'facealpha',.6)
set(fi2,'facealpha',.6)
set(fi3,'facealpha',.6)


% plot all the lines
plot(normal,Parameter1_mean,'k-','LineWidth',2,'Marker','x');
plot(normal,Parameter1_plus_std,'k-')
plot(normal,Parameter1_minus_std,'k-')
scatter(normal,Parameter1_mean,30,FC_TSO_congestion_plot,'filled')


plot(normal,Parameter2_mean,'k--','LineWidth',2,'Marker','x');
plot(normal,Parameter2_plus_std,'k--')
plot(normal,Parameter2_minus_std,'k--')
scatter(normal,Parameter2_mean,30,PCC_TSO_congestion_plot,'filled')

plot(normal,Parameter3_mean,'k:','LineWidth',2,'Marker','x');
plot(normal,Parameter3_plus_std,'k:')
plot(normal,Parameter3_minus_std,'k:')
scatter(normal,Parameter3_mean,30,CM_TSO_congestion_plot,'filled')

c = colorbar;
c.TickLabelInterpreter = 'latex';
% c.TickLabelFormat = '\\textbf{%g}';
c.Label.String = '\textbf{Congestion Level [p.u.]}';
c.Label.Interpreter = 'latex';
% c.Label.FontWeight = 'bold';

% some tweaking
xlim([min(z),max(z)])
axes1 = gca;
grid on
xlabel('\textbf{Wind Penetration [\%]}','Interpreter','latex')
ylabel('\textbf{Expected Social Welfare [\$]}','Interpreter','latex')
% legend('\textbf{Full Coordination (Ideal Benchmark)}','\textbf{PCC Optimizer (proposed coord.)}','\textbf{No Coordination}','Location','best','Interpreter','latex','FontSize',10)
legend('Full Coordination (Ideal Benchmark)','PCC Optimizer (proposed coord.)','No Coordination','Location','best')
set(gca,'FontWeight','bold','FontSize',10)
set(gca,'TickLabelInterpreter','latex')
axes1.XAxis.TickLabelFormat = '\\textbf{%g}';
axes1.YAxis.TickLabelFormat = '\\textbf{%g}';


set(OOS_STD_fig,'PaperSize',[15 11.5]);
saveas(OOS_STD_fig,'./solutions/OOS_results_STD','pdf')

%% plot the congestion level
clc
congestion_level = figure;
plot(wind_penetration_OOS(1,:)*100,PCC_TSO_congestion_plot,'LineWidth',2)
hold on
plot(wind_penetration_OOS(1,:)*100,FC_TSO_congestion_plot,'LineWidth',2)
plot(wind_penetration_OOS(1,:)*100,CM_TSO_congestion_plot,'LineWidth',2)
legend('PCC-Optimizer','Full Coordination','No Coordination','location','best')
xlabel('\textbf{Wind Penetration [\%]}','Interpreter','latex')
ylabel('\textbf{Probability of congestion}','Interpreter','latex')
set(gcf, 'PaperUnits', 'points');
set(congestion_level,'PaperSize',[385 180]);

axes1 = gca;
axes1.XAxis.TickLabelFormat = '\\textbf{%g}';
axes1.YAxis.TickLabelFormat = '\\textbf{%g}';
set(gca,'FontWeight','bold','FontSize',10)
set(gca,'TickLabelInterpreter','latex')
grid on
set(gcf, 'Position',  [100, 100, 500, 220])
saveas(congestion_level,'./solutions/congestion_level','pdf')

%% plot the power flow over the PCC with power flow from the DA-market
k = 1;
% wind_penetration_OOS = wind_penetration_OOS*100;
color_range = 'kbrgmy';
marker_range = 'o+pxsd';
plot_PCC_cap = repmat(base_PCC_capacity,[1 length(wind_penetration_OOS(k,:))]);
DAflow_plot_PCC = figure;
pl_sf = plot(wind_penetration_OOS(k,:)*100,s_flow_PCC_DAflow_expec','LineWidth',2);
for m = 1:length(pl_sf)
	pl_sf(m).Color = color_range(m);
	pl_sf(m).Marker = marker_range(m);
end
hold on
pl_lim = plot(wind_penetration_OOS(k,:)*100,plot_PCC_cap,':','LineWidth',2);
for m = 1:length(pl_lim)
	pl_lim(m).Color = color_range(m);
	pl_lim(m).Marker = marker_range(m);
end
title('Flow in PCCs with DA-dispatch, PCC optimizer')
xlabel('Wind Penetration')
ylabel('Apparent Power Flow [MW]')
legend('PCC1','PCC2','PCC3','PCC4','PCC5')
grid on
set(gca,'FontWeight','bold','FontSize',13)
saveas(DAflow_plot_PCC,['./solutions/DAflow_PCC' num2str(k)],'png')


DAflow_plot_FC = figure;
pl_sf = plot(wind_penetration_OOS(k,:)*100,s_flow_FC_DAflow_expec','LineWidth',2);
for m = 1:length(pl_sf)
	pl_sf(m).Color = color_range(m);
	pl_sf(m).Marker = marker_range(m);
end
hold on
pl_lim = plot(wind_penetration_OOS(k,:)*100,plot_PCC_cap,':','LineWidth',2);
for m = 1:length(pl_lim)
	pl_lim(m).Color = color_range(m);
	pl_lim(m).Marker = marker_range(m);
end
title('Apparent power flow in PCCs with DA-dispatch, FC')
xlabel('Wind Penetration')
ylabel('Apparent Power Flow [MW]')
legend('PCC1','PCC2','PCC3','PCC4','PCC5')
grid on
set(gca,'FontWeight','bold','FontSize',13)
saveas(DAflow_plot_FC,['./solutions/DAflow_FC' num2str(k)],'png')


DAflow_plot_CM = figure;
pl_sf = plot(wind_penetration_OOS(k,:)*100,s_flow_CM_DAflow_expec','LineWidth',2);
for m = 1:length(pl_sf)
	pl_sf(m).Color = color_range(m);
	pl_sf(m).Marker = marker_range(m);
end
hold on
pl_lim = plot(wind_penetration_OOS(k,:)*100,plot_PCC_cap,':','LineWidth',2);
for m = 1:length(pl_lim)
	pl_lim(m).Color = color_range(m);
	pl_lim(m).Marker = marker_range(m);
end
title('Apparent power flow in PCCs with DA-dispatch, CM')
xlabel('Wind Penetration')
ylabel('Apparent Power Flow [MW]')
legend('PCC1','PCC2','PCC3','PCC4','PCC5')
grid on
set(gca,'FontWeight','bold','FontSize',13)
saveas(DAflow_plot_CM,['./solutions/DAflow_CM' num2str(k)],'png')

%% plot the optimal caps for gens and dems
k=1;
for kk = 1:in_sample_iter
	p_gen_tilde_PCC(:,kk) = results_PCCoptim(kk).DA_outcome.p_gen_DA_tilde;
	p_dem_tilde_PCC(:,kk) = results_PCCoptim(kk).DA_outcome.p_dem_DA_tilde;
end
for kk = 1:length(dems_area)
	feeder_name{kk,1} = ['DSO feeder ' num2str(kk)];
end

p_dem_tilde_plot = figure;
hold on
for kk = 1:length(dems_area)
plot(wind_penetration_IS(k,:),sum(p_dem_tilde_PCC(dems_area{kk},:),1),'Marker','x')
end
legend(feeder_name,'Location','best')
saveas(p_dem_tilde_plot,['./solutions/p_dem_tilde_plot' num2str(k)],'png')


p_gen_tilde_plot = figure;
hold on
for kk = 1:length(dems_area)
plot(wind_penetration_IS(k,:),sum(p_gen_tilde_PCC(gens_area{kk},:),1),'Marker','x')
end
legend(feeder_name,'Location','best')
saveas(p_gen_tilde_plot,['./solutions/p_gen_tilde_plot' num2str(k)],'png')

%% Plot the bar plot for the power flow in the PCCs
clc
% for p = 1:length(PCC_branch_id)
% 	for kk = 1:size(s_flow_PCC_DAflow_expec,2)
% 		if s_flow_PCC_DAflow_expec(p,kk) - s_flow_PCC_RT_expec(p,kk) > 0
% 			s_flow_PCC_diff(p,kk) = s_flow_PCC_DAflow_expec(p,kk) - s_flow_PCC_RT_expec(p,kk);
% 		elseif s_flow_PCC_DAflow_expec - s_flow_PCC_RT_expec < 0
% 			s_flow_PCC_diff(p,kk) = s_flow_PCC_RT_expec(p,kk) - s_flow_PCC_DAflow_expec(p,kk);
% 		end
% 	end
% end
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
colors_avail = 'rgb';%[1 0 0; 0.8 0 0; 0.5 0 0];
% set(groot, 'defaultLabelInterpreter','latex');
for pccn = 1:length(PCC_branch_id)
	
	
	
	n_plot_st = 6;
	n_plot = 4;
	n_plot_en = 0;
	
	wpp = wind_penetration_OOS(k,n_plot_st:n_plot:end-n_plot_en);
	s_pcc = s_flow_PCC_RT_expec(pccn,n_plot_st:n_plot:end-n_plot_en);
	s_FC = s_flow_FC_RT_expec(pccn,n_plot_st:n_plot:end-n_plot_en);
	s_CM = s_flow_CM_RT_expec(pccn,n_plot_st:n_plot:end-n_plot_en);
	
	s_pcc_DA = s_flow_PCC_DAflow_expec(pccn,n_plot_st:n_plot:end-n_plot_en);
	s_FC_DA = s_flow_FC_DAflow_expec(pccn,n_plot_st:n_plot:end-n_plot_en);
	s_CM_DA = s_flow_CM_DAflow_expec(pccn,n_plot_st:n_plot:end-n_plot_en);
	
	x_axis_tick = wind_penetration_OOS(k,n_plot_st:n_plot:end-n_plot_en);
	bar_plot_PCC_flow = figure;
	hold on

	w2 = .8;
	bl = bar(wpp,  [s_pcc_DA;s_FC_DA; s_CM_DA]',w2);
	set(bl(1), 'FaceColor',[0.1 0.3 0.7])
	set(bl(2), 'FaceColor',[0.1 0.6 0.9])
	set(bl(3), 'FaceColor',[0.6 0.9 0.99])
	
	hb = bar(wpp,  [s_pcc; s_FC; s_CM]',0.25);
	set(hb(1), 'FaceColor','k')
	set(hb(2), 'FaceColor','k')
	set(hb(3), 'FaceColor','k')
	set(hb(1),'HandleVisibility','off')
	set(hb(2),'HandleVisibility','off')
	
	hold off
	pause(0.1)
	clear ctr ydt
	for k1 = 1:length(hb)
		ctr(k1,:) = bsxfun(@plus, hb(1).XData, [hb(k1).XOffset]')';
		ydt(k1,:) = hb(k1).YData;
	end
	hold on
	clear hb
	std_pcc_flow = STD_pcc_flow_PCCO_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	std_FC_flow = STD_pcc_flow_FC_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	std_CM_flow = STD_pcc_flow_CM_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	
	
	quant09_pcc_flow = quant09_pcc_flow_PCCO_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	quant09_FC_flow = quant09_pcc_flow_FC_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	quant09_CM_flow = quant09_pcc_flow_CM_RT(pccn,n_plot_st:n_plot:end-n_plot_en);

	quant01_pcc_flow = quant01_pcc_flow_PCCO_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	quant01_FC_flow = quant01_pcc_flow_FC_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	quant01_CM_flow = quant01_pcc_flow_CM_RT(pccn,n_plot_st:n_plot:end-n_plot_en);
	for kk = 1:length(std_pcc_flow)	
		p1 = plot([ctr(1,kk) ctr(1,kk)],[quant09_pcc_flow(kk) quant01_pcc_flow(kk)],'Marker','square','LineWidth',2,'Color','r',...
			'MarkerIndices',[1 2]);
		p2 = plot([ctr(2,kk) ctr(2,kk)],[quant09_FC_flow(kk) quant01_FC_flow(kk)],'Marker','square','LineWidth',2,'Color','r');
		p3 = plot([ctr(3,kk) ctr(3,kk)],[quant09_CM_flow(kk) quant01_CM_flow(kk)],'Marker','square','LineWidth',2,'Color','r');
	end
	
%  	errorbar(ctr,ydt,[std_pcc_flow; std_FC_flow; std_CM_flow],'.','LineWidth',1.5,'Color','r')
	
	
	
	
	% pl_lim = plot(wind_penetration_OOS(k,:),plot_PCC_cap(3,:),':','LineWidth',2);
	hline(base_PCC_capacity(pccn),'m:','Physical Limit');
	xlabel('Wind Penetration [\%]','Interpreter','latex','FontSize',12)
	ylabel('Apparent Power Flow in PCC [MVA]','Interpreter','latex','FontSize',12)
	leg = legend('DA - PCC Optimizer','DA - Full Coordination','DA - No Coordination','Expect. RT Power Flow','0.9 and 0.1 quantile','location','northeast');
	ax = gca;
	grid on
	ax.XTick = x_axis_tick; 
	ax.XTickLabels = round(x_axis_tick*1000)/1000*100;
	set(ax,'FontSize',12)


	set(bar_plot_PCC_flow,'PaperSize',[15 11.5]);
	saveas(bar_plot_PCC_flow,['./solutions/bar_flow_PCC' num2str(pccn)],'pdf')
	hold off

end

%% plot the prices
lambda_RT_PCCO = zeros([length(wind_factor) nb]);
for kk = 1:outofsample_iter
	lambda_DA_PCCO(ii(kk),jj(kk),pp(kk)) = results_PCCoptim_OOS(kk).DA_outcome.lambda;
	for s = 1:nscen_OOS
		lambda_RT_PCCO(kk,:) = lambda_RT_PCCO(kk,:) + results_PCCoptim_OOS(kk).RT_outcome(s).lambda_RT'/nscen_OOS;
	end
end

lambda_RT_PCCO = squeeze(lambda_RT_PCCO);

price_plot = figure;
plot(wind_penetration_OOS,lambda_DA_PCCO)
hold on
plot(wind_penetration_OOS,-lambda_RT_PCCO(:,[15 25 26 27 28 29]))

legend('DA-price','RT Node 15','RT node 25','RT node 26','RT node 27','RT node 28','RT node 29')
set(price_plot,'PaperSize',[15 11.5]);
saveas(price_plot,['./solutions/price_plots' num2str(k)],'pdf')

%% plot the loads and gens cleard in in DA

for kk = 1:in_sample_iter
	p_dem_DA(kk,:) = results_PCCoptim(kk).DA_outcome.p_dem_DA;
	p_gen_DA(kk,:) = results_PCCoptim(kk).DA_outcome.p_gen_DA;
end
for d = 1:nd
	dem_names{d} = ['Demand ' num2str(d)];
end

dem_plot = figure;

plot(wind_penetration_IS,p_dem_DA)
legend(dem_names)


for g = 1:ng
	gen_names{g} = ['Gen. ' num2str(g)];
end

gen_plot = figure;

plot(wind_penetration_IS,p_gen_DA)
legend(gen_names)

%% plot the aggregated caps as a function of the varying line capacity
clear p_gen_tilde_PCC p_dem_tilde_PCC
clc
close all

k=1;
for kk = 1:in_sample_iter
	p_gen_tilde_PCC(:,ii(kk),jj(kk),pp(kk)) = results_PCCoptim(kk).DA_outcome.p_gen_DA_tilde;
	p_dem_tilde_PCC(:,ii(kk),jj(kk),pp(kk)) = results_PCCoptim(kk).DA_outcome.p_dem_DA_tilde;
end
p_gen_tilde_PCC = squeeze(p_gen_tilde_PCC);
p_dem_tilde_PCC = squeeze(p_dem_tilde_PCC);
for kk = 1:length(dems_area)
	feeder_name{kk,1} = ['DSO feeder ' num2str(kk)];
end
feeder_name{kk+1,1} = ['Social Welfare'];
for k = 1:max(jj)
	p_dem_tilde_plot = figure;
	hold on
	for kk = [1 2 3 4 5]%1:length(dems_area)
		plot(base_varLine_capacity*line_factor,squeeze(sum(p_dem_tilde_PCC(dems_area{kk},k,:),1)),'LineWidth',2)
	end
	xlabel('\textbf{Physical Capacity of PCC 1 [MW]}','Interpreter','latex')
	ylabel('\textbf{Accumulated optimal caps on demands [MW]}','Interpreter','latex')
	set(p_dem_tilde_plot,'PaperSize',[15.5 11.5]);
	
	
	yyaxis right
	plot(base_varLine_capacity*line_factor,-squeeze(pcc_optim_cost(1,k,:)),'LineStyle','--','LineWidth',2)
	ylabel('\textbf{Social Welfare [\$]}','Interpreter','latex')
	legend(feeder_name{[1 2 3 4 5 6]},'Location','best')

	
	axes1 = gca;
	axes1.XAxis.TickLabelFormat = '\\textbf{%g}';
	axes1.YAxis(1).TickLabelFormat = '\\textbf{%g}';
	axes1.YAxis(2).TickLabelFormat = '\\textbf{%g}';
	axes1.YAxis(2).Color = 'k';
	set(gca,'FontWeight','bold','FontSize',10)
	set(gca,'TickLabelInterpreter','latex')
	
	title(['\textbf{Wind Penetration ' num2str(round(wind_penetration_IS(1,k,1)*1000)/10) '\%}'],'Interpreter','latex');
	grid on
	saveas(p_dem_tilde_plot,['./solutions/var_line_cap_DSO_dem_caps' num2str(k)],'pdf')
	

% 	p_gen_tilde_plot = figure;
% 	hold on
% 	for kk = 1:length(dems_area)
% 		plot(base_varLine_capacity*line_factor,squeeze(sum(p_gen_tilde_PCC(gens_area{kk},k,:),1)),'Marker','x')
% 	end
% 	legend(feeder_name,'Location','best')
% 	saveas(p_gen_tilde_plot,['./solutions/var_line_cap_DSO_gen_caps' num2str(k)],'png')


end
