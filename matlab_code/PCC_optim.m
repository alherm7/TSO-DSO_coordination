function [result_pccoptim] = PCC_optim(data,kk)

tolerance = 0.001;
cc{1} = double.empty(1,0);
[~,~,~,~,~,~,~,~,~,~,~,...
	ng,nd,~,~,~,~,~,~,~,~,...
    Pmax_gen,~,~,~,Pmax_dem,~,~,...
	~,~,~,~,~,~,~,~,~,...
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


time_master_solve = [];
time_subprob_solve = [];
for iter = 1:100 % Benders Iterations, set your desired max iteration number here
	disp(['WPP: ' num2str(kk) ', PCC Optimizer, Iteration ' num2str(iter)]);

	t_master = tic;
	[DA_outcome] = PCC_optimizer_DA_KKT_simple(data,dual_DA_gen,dual_DA_dem,...
		dual_day_ahead_wind,cost_RT,iter,p_gen_DA_hat,p_dem_DA_hat,wind_DA_hat,kk);
	time_master = toc(t_master);
	disp(['WPP: ' num2str(kk) ', PCC Optimizer, Benders It.: ' num2str(iter) ', Master Objective: ' num2str(DA_outcome.cost) '; Solver Time : ' num2str(time_master)]);

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
	disp(['WPP: ' num2str(kk) ', PCC Optimizer, Benders It.: ' num2str(iter) ', Solution Time for Sub-Problem: ' num2str(time_sub)]);
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
	
	benders_it1 = figure(5);
	plot(1:iter,[prob_wscen(1,:)*cuts ; prob_wscen(1,:)*objective_sub])
	xlabel('Iterations')
	legend('Mean Benders Cut','Subproblem Objective','Location','best')
	grid on
% 	saveas(benders_it1,['./solutions/benders_it/bender_it1_iteration' num2str(kk)],'fig')

    
	benders_it2 = figure(6);
	plot(1:iter,[objective_master; prob_wscen(1,:)*objective_sub + cost_DA])
    xlabel('Iterations')
	legend('Master Objective','Subproblem Objective','Location','best')
	grid on
% 	saveas(benders_it2,['./solutions/benders_it/bender_it2_iteration' num2str(kk)],'fig')

    time_master_solve(iter) = DA_outcome.fin_time_master_solve;
    for s = 1:nscen
       time_subprob_solve(iter,s) = RT_outcome(s).tfin_subprob;
    end
% 	distri_fig;
	time_benders = toc(time_all);
    % prob_wscen(1,:)*cuts(:,iter) - prob_wscen(1,:)*objective_sub(:,iter)
	relative_gap = abs(objective_master(iter) - ( prob_wscen(1,:)*objective_sub(:,iter) + cost_DA(iter)))...
		/abs(objective_master(iter));
	disp(['WPP: ' num2str(kk) ', PCC Optimizer, Benders It.: ' num2str(iter)  ', Relative Gap: ' num2str(relative_gap)]);
	if  relative_gap < tolerance
		disp(['WPP: ' num2str(kk) ', PCC Optimizer, Benders algorithm has converged to within tolerance after ' num2str(iter) ' iterations and ' num2str(time_benders) ' seconds.'])
		break;
	end
% 	DA_outcome.p_gen_DA - DA_outcome_KKT.p_gen_DA
% 	DA_outcome.p_dem_DA - DA_outcome_KKT.p_dem_DA
	
end
t_solver_total = sum(time_master_solve) + sum(sum(time_subprob_solve));


result_pccoptim.DA_outcome = DA_outcome;
result_pccoptim.RT_outcome = RT_outcome;
result_pccoptim.t_solver_total = t_solver_total;
result_pccoptim.time_master_solve = time_master_solve;
result_pccoptim.time_subprob_solve = time_subprob_solve;
% result_pccoptim.result_DSO_lookahead = result_DSO_lookahead;
% result_pccoptim.result_lookahead_single = result_lookahead_single;



