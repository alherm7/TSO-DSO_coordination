function [result_pccoptim] = market_clearing_DA_PF(data,kk)

tolerance = 0.001;
cc{1} = double.empty(1,0);
[~,~,~,~,~,~,~,~,~,~,~,...
	ng,nd,~,~,~,~,~,~,~,~,...
    Pmax_gen,~,~,~,Pmax_dem,~,~,...
	~,~,fbusL,tbusL,~,~,~,~,~,...
    ~,~,~,~,~,~,~,~,...
	~,~,~,~,...
	~,~,~,~,~,~,....
	~,~,~,~,~,prob_wscen,~,~,...
	~,nscen,~,bus_wgen,...
	offer_wind_DA,offer_wind_up,offer_wind_dn] = Data_Reader(data,cc);

if ~isempty(data.DR_DSO)
	bus_DR = data.DR_DSO(:,1);
	nDR = length(bus_DR);
else
	bus_DR = [];
	nDR = length(bus_DR);
end

nw = length(bus_wgen);

dual_DA_gen{1} = zeros(ng,nscen);
dual_DA_dem{1} = zeros(nd,nscen);
dual_DA_DR{1} = zeros(nDR,nscen);
dual_day_ahead_wind{1} = zeros(nw,nscen);

p_gen_DA_hat{1} = zeros(ng,1);
p_dem_DA_hat{1} = zeros(nd,1);
p_DR_DA_hat{1} = zeros(nDR,1);
wind_DA_hat{1} = zeros(nw,1);

cost_RT{1} = -100000 * ones(nscen,1);
time_all = tic;

p_gen_DA_tilde = Pmax_gen; %[200; 3; 2.5];
p_dem_DA_tilde = Pmax_dem; %[200; 200; 5; 6];
for iter = 1:1
	disp(['WPP PF market: ' num2str(kk) ', PF Market, Iteration ' num2str(iter)]);

	t_master = tic;
% 	[DA_outcome] = DA_market_KKT(data, dual_DA_gen, dual_DA_dem, dual_day_ahead_wind, cost_RT,... 
% 		iter, p_gen_DA_hat, p_dem_DA_hat, wind_DA_hat, p_gen_DA_tilde, p_dem_DA_tilde,kk);
    [DA_outcome] = DA_market_clearing_PF(data, dual_DA_gen, dual_DA_dem, dual_day_ahead_wind, cost_RT,... 
		iter, p_gen_DA_hat, p_dem_DA_hat, wind_DA_hat, p_gen_DA_tilde, p_dem_DA_tilde,kk);

	time_master = toc(t_master);
	disp(['WPP PF market: ' num2str(kk) ', PF Market, Benders It.: ' num2str(iter) ', Master Objective: ' num2str(DA_outcome.cost) '; Solver Time : ' num2str(time_master)]);


	p_gen_DA_hat{iter+1} = DA_outcome.p_gen_DA;
	p_dem_DA_hat{iter+1} = DA_outcome.p_dem_DA;
	p_DR_DA_hat{iter+1} = p_DR_DA_hat{iter};
	wind_DA_hat{iter+1} = DA_outcome.wind_DA;

		
	p_gen_hat_single = p_gen_DA_hat{iter+1};
	p_dem_hat_single = p_dem_DA_hat{iter+1};
	p_DR_hat_single = p_DR_DA_hat{iter+1};
	wind_DA_hat_single = wind_DA_hat{iter+1};

	t_sub = tic;
	parfor s = 1:nscen
		RT_outcome(s) = RT_TSO_DSO(data,p_gen_hat_single,p_dem_hat_single,p_DR_hat_single,wind_DA_hat_single,s,iter,kk);
	end
	time_sub = toc(t_sub);
	disp(['Solution Time for Sub-Problem: ' num2str(time_sub)]);
	disp(['WPP: ' num2str(kk) ', Conventional Market, Benders It.: ' num2str(iter) ', Solution Time for Sub-Problem: ' num2str(time_sub)]);
	for s = 1:nscen
		dual_DA_gen{iter+1}(:,s) = RT_outcome(s).dual_day_ahead_gen;
		dual_DA_dem{iter+1}(:,s) = RT_outcome(s).dual_day_ahead_dem;
		dual_DA_DR{iter+1}(:,s) = RT_outcome(s).dual_day_ahead_DR;
		dual_day_ahead_wind{iter+1}(:,s) = RT_outcome(s).dual_day_ahead_wind;
		cost_RT{iter+1}(s,1) = RT_outcome(s).cost_RT;
	end
	
	
	
	objective_master(iter) = DA_outcome.cost;
	cost_DA(iter) = DA_outcome.cost_DA;
	objective_sub(:,iter) = cost_RT{iter+1};
	
	cuts(:,iter) = DA_outcome.alpha_cut;
	


	time_benders = toc(time_all);
    % prob_wscen(1,:)*cuts(:,iter) - prob_wscen(1,:)*objective_sub(:,iter)
	relative_gap = abs(objective_master(iter) - ( prob_wscen(1,:)*objective_sub(:,iter) + cost_DA(iter)))...
		/abs(objective_master(iter));
	disp(['WPP: ' num2str(kk) ', Conv. Market, Benders It.: ' num2str(iter)  ', Relative Gap: ' num2str(relative_gap)]);
% 	if  relative_gap < tolerance
% 		disp(['WPP: ' num2str(kk) ', Conventional Market, Benders algorithm has converged to within tolerance after ' num2str(iter) ' iterations and ' num2str(time_benders) ' seconds.'])
% 		break;
%     end

    cost = prob_wscen(1,:)*objective_sub(:,iter) + cost_DA(iter);
    for k = 1: nscen
        p_flow(:,:,k) = RT_outcome(k).p_flow;
    end
	p_flow_mean= mean([p_flow],3);
    p_flow_ft = p_flow_mean(sub2ind(size(mean(p_flow,3)),fbusL,tbusL));
end
result_pccoptim.DA_outcome = DA_outcome;
result_pccoptim.RT_outcome = RT_outcome;
result_pccoptim.cost = cost;
result_pccoptim.p_flow_ft = p_flow_ft;



