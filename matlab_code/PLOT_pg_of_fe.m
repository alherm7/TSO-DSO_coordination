%% plot p_g_tilde and p_d_tilde with respect to the pcc limits and prices.
clear all
close all
clc

test_case = 'casedata_6bus_1scen';
mpc = eval(test_case);

fe_up = -7:0.3:1;
fe_dn = [-10];
pi_pcc_DA = 5:1:7;
pi_pcc_up = [3]; % 0.1:0.4:3;
pi_pcc_dn = [4]; %0.1:0.4:3;
e = 1;

direc = 'min';

[x1,x2] = ndgrid(fe_up,pi_pcc_DA);

% [xq,yq] = ndgrid(linspace(20,49,200),linspace(10,33,200));

% vq = F(xq,yq);

% figure(2)
% [C,h] = contourf(xq,yq,vq,40,'edgecolor','none');
% k=1;
% [result_lookahead] = DSO_lookahead_mccormcick_cvx_obj(mpc,x1(k),fe_dn,...
% 	x2(k),pi_pcc_up,pi_pcc_dn,e,[]);
npart = 2;
tic
for k = 1:numel(x1)
	display(['Iteration ' num2str(k)]);
	result_lookahad(k) = DSO_lookahead_KKT(mpc,x1(k),fe_dn,...
		x2(k),pi_pcc_up,pi_pcc_dn,e,npart,direc);
end
t_sim = toc
for k = 1:numel(x1)
	p_dem_DA_1(k) = result_lookahad(k).p_dem_DA(1);
	p_dem_DA_2(k) = result_lookahad(k).p_dem_DA(2);
	p_gen_DA(k) = result_lookahad(k).p_gen_DA;
end
p_dem_DA_1 = reshape(p_dem_DA_1,size(x1));
p_dem_DA_2 = reshape(p_dem_DA_2,size(x1));
p_gen_DA = reshape(p_gen_DA,size(x1));

% p_dem_1 = griddedInterpolant(x1,x2,p_dem_DA_1,'linear');
% p_dem_2 = griddedInterpolant(x1,x2,p_dem_DA_2,'linear');

figure(1)
surf(x1,x2,p_dem_DA_1)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d1}^{DA}')

figure(2)
surf(x1,x2,p_dem_DA_2)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d2}^{DA}')

figure(3)
surf(x1,x2,p_gen_DA)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{g}^{DA}')


%%
npart = 20;
tic
for k = 1:numel(x1)
	display(['Iteration ' num2str(k)]);
	result_lookahad(k) = DSO_lookahead_KKT(mpc,x1(k),fe_dn,...
		x2(k),pi_pcc_up,pi_pcc_dn,e,npart,direc);
end
t_sim = toc
for k = 1:numel(x1)
	p_dem_DA_1(k) = result_lookahad(k).p_dem_DA(1);
	p_dem_DA_2(k) = result_lookahad(k).p_dem_DA(2);
	p_gen_DA(k) = result_lookahad(k).p_gen_DA;
end
p_dem_DA_1 = reshape(p_dem_DA_1,size(x1));
p_dem_DA_2 = reshape(p_dem_DA_2,size(x1));
p_gen_DA = reshape(p_gen_DA,size(x1));

% p_dem_1 = griddedInterpolant(x1,x2,p_dem_DA_1,'linear');
% p_dem_2 = griddedInterpolant(x1,x2,p_dem_DA_2,'linear');

figure(4)
surf(x1,x2,p_dem_DA_1)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d1}^{DA}')

figure(5)
surf(x1,x2,p_dem_DA_2)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d2}^{DA}')

figure(6)
surf(x1,x2,p_gen_DA)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{g}^{DA}')

%%
direc = 'max';
npart = 2;
tic
for k = 1:numel(x1)
	display(['Iteration ' num2str(k)]);
	result_lookahad(k) = DSO_lookahead_KKT(mpc,x1(k),fe_dn,...
		x2(k),pi_pcc_up,pi_pcc_dn,e,npart,direc);
end
t_sim = toc
for k = 1:numel(x1)
	p_dem_DA_1(k) = result_lookahad(k).p_dem_DA(1);
	p_dem_DA_2(k) = result_lookahad(k).p_dem_DA(2);
	p_gen_DA(k) = result_lookahad(k).p_gen_DA;
end
p_dem_DA_1 = reshape(p_dem_DA_1,size(x1));
p_dem_DA_2 = reshape(p_dem_DA_2,size(x1));
p_gen_DA = reshape(p_gen_DA,size(x1));

% p_dem_1 = griddedInterpolant(x1,x2,p_dem_DA_1,'linear');
% p_dem_2 = griddedInterpolant(x1,x2,p_dem_DA_2,'linear');

figure(7)
surf(x1,x2,p_dem_DA_1)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d1}^{DA}')

figure(8)
surf(x1,x2,p_dem_DA_2)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d2}^{DA}')

figure(9)
surf(x1,x2,p_gen_DA)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{g}^{DA}')


%%
direc = 'max';
npart = 20;
tic
for k = 1:numel(x1)
	display(['Iteration ' num2str(k)]);
	result_lookahad(k) = DSO_lookahead_KKT(mpc,x1(k),fe_dn,...
		x2(k),pi_pcc_up,pi_pcc_dn,e,npart,direc);
end
t_sim = toc
for k = 1:numel(x1)
	p_dem_DA_1(k) = result_lookahad(k).p_dem_DA(1);
	p_dem_DA_2(k) = result_lookahad(k).p_dem_DA(2);
	p_gen_DA(k) = result_lookahad(k).p_gen_DA;
end
p_dem_DA_1 = reshape(p_dem_DA_1,size(x1));
p_dem_DA_2 = reshape(p_dem_DA_2,size(x1));
p_gen_DA = reshape(p_gen_DA,size(x1));

% p_dem_1 = griddedInterpolant(x1,x2,p_dem_DA_1,'linear');
% p_dem_2 = griddedInterpolant(x1,x2,p_dem_DA_2,'linear');

figure(10)
surf(x1,x2,p_dem_DA_1)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d1}^{DA}')

figure(11)
surf(x1,x2,p_dem_DA_2)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{d2}^{DA}')

figure(12)
surf(x1,x2,p_gen_DA)
xlabel('Max Import at PCC')
ylabel('Day ahead price at PCC')
zlabel('P_{g}^{DA}')
%%

fup_test = 0.7;
pi_DA_test = 4;
p_dem_2 = 3.2;
test_if_feas = DSO_lookahead_KKT(mpc,fup_test,fe_dn,...
		pi_DA_test,pi_pcc_up,pi_pcc_dn,e,npart,direc);