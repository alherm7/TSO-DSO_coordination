function [nb,busL,busN,busType,Pd,Qd,Vmax,Vmin,statusG,activeG,activeD,...
	ng,nd,genL,genN,incidentG,demL,incidentD,incidentW,Qmax_gen,Qmin_gen,...
    Pmax_gen,Pmin_gen,Qmax_dem,Qmin_dem,Pmax_dem,Pmin_dem,statusL,...
	activeL,nl,fbusL,tbusL,SlmMax,fbusN,tbusN,incidentF,incidentT,...
    Yf,Yt,YfP,YtP,Ybus,edges,offer_gen_DA,offer_gen_upreg,...
	offer_gen_downreg,bid_dem_DA,bid_dem_upreg,bid_dem_downreg,...
	nr_vals_cost,p,cost_type,ys,areas,bus_area,....
	incident_singlephase,Y_ft,rBranch,xBranch,narea,prob_wscen,Wmax,Wmin,...
	n_wgen,nscen,Wmax_mean_DA,bus_wgen,...
	offer_wind_DA,offer_wind_up,offer_wind_dn] = Data_Reader(mpc,cc)
     
	version = mpc.version;
	
    nc = length(cc);
% 	Z_base = mpc.basekV^2/(mpc.baseMVA*1000);
    %Importing bus parameres
    nb = size(mpc.bus,1);
    busL = mpc.bus(:,1);
    busN = sparse(busL, ones(nb,1), (1:nb) );
    busType = mpc.bus(:,2);
    Pd = mpc.bus(:,3)/mpc.baseMVA;
    Qd = mpc.bus(:,4)/mpc.baseMVA;

    Gs = mpc.bus(:,5)/mpc.baseMVA;
    Bs = mpc.bus(:,6)/mpc.baseMVA;
    Vmax = mpc.bus(:,12);
    Vmin = mpc.bus(:,13);

    Ysh = Gs + Bs*1i;

    %Importing generator parameres
    statusG = mpc.gen(:,8);
    activeG = find(statusG);

	statusD = mpc.dem(:,8);
    activeD = find(statusD);
	
    ng = size(mpc.gen,1);
	nd = size(mpc.dem,1);
	
    genL = mpc.gen(:,1);
    genN = busN(genL);
    incidentG = diag(statusG)*sparse(1:ng, genN, 1 , ng, nb);
	
	demL = mpc.dem(:,1);
    demN = busN(demL);
    incidentD = diag(statusD)*sparse(1:nd, demN, 1 , nd, nb);
	
    Qmax_gen = (mpc.gen(:,4).*statusG)/mpc.baseMVA;
    Qmin_gen = (mpc.gen(:,5).*statusG)/mpc.baseMVA;
    Pmax_gen = (mpc.gen(:,9).*statusG)/mpc.baseMVA;
    Pmin_gen = (mpc.gen(:,10).*statusG)/mpc.baseMVA;
	
	Qmax_dem = (mpc.dem(:,4).*statusD)/mpc.baseMVA;
    Qmin_dem = (mpc.dem(:,5).*statusD)/mpc.baseMVA;
	if version == '2'
		Pmax_dem = (mpc.dem(:,9).*statusD)/mpc.baseMVA;
		Pmin_dem = (mpc.dem(:,10).*statusD)/mpc.baseMVA;
	elseif version == '3'
		max_syst_load = mpc.max_syst_load;
		Pmax_dem = (mpc.dem(:,9).*statusD)/mpc.baseMVA .* max_syst_load/100;
		Pmin_dem = (mpc.dem(:,10).*statusD)/mpc.baseMVA .* max_syst_load/100;
	end
	

    %Importing branch parameres
    statusL = mpc.branch(:,9);
    activeL = find(statusL);

    nl = size(mpc.branch,1);

    fbusL = mpc.branch(:,1);
    tbusL = mpc.branch(:,2);
    rBranch = mpc.branch(:,3);
    xBranch = mpc.branch(:,4);
    bBranch = mpc.branch(:,5);
    SlmMax = mpc.branch(:,6)/mpc.baseMVA;

    ratio = mpc.branch(:,7);
    phase = mpc.branch(:,8);

    ratio(ratio == 0) = 1;

    ys = 1 ./ (rBranch + xBranch*1i);% + 2*10^(-7)*ones(nl,1);
    sh = exp(phase*(pi/180)*1i);

    fbusN = busN(fbusL);
    tbusN = busN(tbusL);

    incidentF = diag(statusL) * sparse(1 : nl, fbusN, 1, nl, nb);
    incidentT = diag(statusL) * sparse(1 : nl, tbusN, 1, nl, nb);

    Yff = (ys + (bBranch/2)*1i) ./ (ratio.^2);
    Yft = -(ys ./ ratio) .* sh;
    Ytf = -(ys ./ ratio) ./ sh;
    Ytt = (ys + (bBranch/2) * 1i);

    Yf = sparse(diag(Yff) * incidentF + diag(Yft) * incidentT);
    Yt = sparse(diag(Ytf) * incidentF + diag(Ytt) * incidentT);

    YffP = ys ./ (ratio .^ 2);             YftP = -ys ./ (ratio .* sh);
    YtfP = -ys ./ (ratio ./ sh);           YttP = ys;
    
    YfP = sparse(diag(YffP) * incidentF + diag(YftP) * incidentT);
    YtP = sparse(diag(YtfP) * incidentF + diag(YttP) * incidentT);

    Ybus = cell(nc,1);
    edges = zeros(nl,nc);
    for rr = 1 : nc
        ex = setdiff(1 : nl, cc{rr});
        ex2 = setdiff(activeL, cc{rr});
        edges(ex2,rr) = 1;
        Ybus{rr} = incidentF(ex, :)' * Yf(ex, :) + incidentT(ex, :)' * Yt(ex, :) + diag(Ysh);
	end

	for k = 1:size(mpc.gencost,1)
		nr_vals_cost = mpc.gencost(k,4);
		if size(mpc.gencost,1) ~= size(mpc.gen,1)
			error('Gencost and gen are not the same size')
		end
		
		if mpc.gencost(k,1) == 2
			for m = 1:nr_vals_cost
				offer_gen_DA(k,m) = (mpc.gencost(k,4+m).*statusG(k))*mpc.baseMVA^(nr_vals_cost-m);
			end
			p = double.empty(1,0);
		end
		
		if mpc.gencost(k,1) == 1
			for m = 1:nr_vals_cost
				offer_gen_DA(k,m) = (mpc.gencost(k,4+m*2-1).*statusG(k))*mpc.baseMVA;
				p(k,m) = (mpc.gencost(k,4+m*2).*statusG(k))*mpc.baseMVA;
			end
		end
	end
	cost_type = mpc.gencost(:,1);
    %Importing cost parameres

    offer_gen_upreg = mpc.gencost_RT(:,2);
	offer_gen_downreg = mpc.gencost_RT(:,3);
	
	bid_dem_DA = mpc.demcost(:,6);
	
	bid_dem_upreg = mpc.demcost_RT(:,2);
	bid_dem_downreg = mpc.demcost_RT(:,3);
	
    ident = eye(nb);
    for k=1:length(ys)
        ylm{k} = (ys(k)'+ys(k))*ident(:,fbusL(k))*ident(fbusL(k),:) - ys(k)*ident(:,fbusL(k))*ident(tbusL(k),:);
        
        Ylm{k} = 1/2*[real(ylm{k}+ylm{k}.'), imag(ylm{k}.' - ylm{k}); imag(ylm{k} - ylm{k}.'), real(ylm{k} + ylm{k}.')];
    end
    
    bus_area = mpc.bus(:,7);
	narea = length(unique(bus_area));
	areacode = unique(bus_area);
	
	for k = 1:narea
		areas{k} = busL(bus_area == areacode(k));
	end
	
	incident_singlephase = zeros(nb,nb);
	for k = 1:nl
		incident_singlephase(fbusL(k),tbusL(k)) = true;
	end
	incident_singlephase = incident_singlephase + incident_singlephase';
	
	for k = 1:nl
		Y_ft{fbusL(k),tbusL(k)} = ys(k);
		Y_ft{tbusL(k),fbusL(k)} = ys(k);
	end
	
	n_wgen = length(mpc.RTscen(:,1));
	nscen = mpc.RTscen(1,2);
	for k = 1:n_wgen
		if nscen ~= mpc.RTscen(k,2)
			error('Number of real-time scenarios must be equal for all wind generators');
		end
		for m = 1:nscen
			prob_wscen(k,m) = mpc.RTscen(k,2+3*m);
			Wmax(k,m) = mpc.RTscen(k,3*m);
			Wmin(k,m) = mpc.RTscen(k,1+3*m);
		end
	end
	
	Wmax_mean_DA = mean(prob_wscen(1,:)*Wmax',1)';
	bus_wgen = mpc.RTscen(:,1);
	
	winN = busN(bus_wgen);
    incidentW = sparse(1:n_wgen, winN, 1 , n_wgen, nb);
	
	offer_wind_DA = mpc.wind_offer(:,1);
	offer_wind_up = mpc.wind_offer(:,2);
	offer_wind_dn = mpc.wind_offer(:,3);
	
end
