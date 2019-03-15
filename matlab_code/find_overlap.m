function [area_codes, n_areas, DSO_codes, fbusL_area, tbusL_area,...
	branch_num_area, fbusL_ext_area, tbusL_ext_area, branch_num_ext_area, overlap,...
	overlap_t, neigh_area, ext_area,...
	ext_area_singlephase,fbus_local,tbus_local,incidence_area,...
	PCC_branch_id] = find_overlap(bus_areacodes, fbusL, tbusL, areas_buses)


ftbus = [fbusL,tbusL];

area_codes = unique(bus_areacodes);
n_areas = length(area_codes);
DSO_codes = area_codes(area_codes > 0);

incidence_ftbus = sparse(fbusL,tbusL,1,length(fbusL),length(fbusL)) + sparse(tbusL,fbusL,1,length(fbusL),length(fbusL));


for k = 1:length(areas_buses)
	ind = 1;
	for m = 1:length(areas_buses{k})
		[line_num,~] = find(ftbus == areas_buses{k}(m));
		lines_zone{k}(ind:ind+length(line_num)-1) = line_num;
		ind = ind + length(line_num);
	end
	
	fb_ind=ismember(fbusL,areas_buses{k});
	tb_ind=ismember(tbusL,areas_buses{k});
	
	branch_ind = fb_ind & tb_ind;
	
	branch_num_area{k} = find(branch_ind);
	fbusL_area{k} = fbusL(branch_num_area{k});
	tbusL_area{k} = tbusL(branch_num_area{k});
	
	lines_zone{k} = unique(lines_zone{k});
	fbusL_ext_area{k} = fbusL(lines_zone{k});
	tbusL_ext_area{k} = tbusL(lines_zone{k});
	branch_num_ext_area{k} = lines_zone{k};
	
	clear ind
	
	
	incidence_area{k} = full(incidence_ftbus(areas_buses{k},areas_buses{k}));
	[fbus_local{k},tbus_local{k}] = find(triu(incidence_area{k}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:n_areas
    for k = 1:length(areas_buses{m})
        if isempty(find( fbusL == areas_buses{m}(k), 1 ) )
            fbus_area{m,k} = [];
        else
            fbus_area{m,k} = find( fbusL == areas_buses{m}(k) );
            fbus_area{m,k} = fbus_area{m,k}(:)';
        end
        
        if isempty(find( tbusL == areas_buses{m}(k), 1 ) )
            tbus_area{m,k} = [];
        else
            tbus_area{m,k} = find( tbusL == areas_buses{m}(k) );
            tbus_area{m,k} = tbus_area{m,k}(:)';
        end
    end
end
ind = 0;
for m = 1:n_areas-1
    for k = m+1:n_areas
        uni = intersect(cell2mat(fbus_area(m,:)),cell2mat(tbus_area(k,:)));
        
        if ~isempty(uni)
            ind = ind +1;
            overlap{m,k} = [fbusL(uni); tbusL(uni)];
            overlap{k,m} = overlap{m,k};
        end
    end
end

for k = 1:n_areas
    neigh_area{k} = [];
	if exist('overlap')
		for m = 1:n_areas
			neigh_area{k} = [neigh_area{k} ~isempty(overlap{k,m})];

		end
		neigh_area{k} = find(neigh_area{k});
	else
		overlap = [];
	end
end

PCC_line_busnum = cell2mat(overlap(1,2:end))';
PCC_branch_id = find(ismember([fbusL tbusL], PCC_line_busnum,'rows'));

for m = 1:n_areas
    for k = 1:length(areas_buses{m})
        if isempty(find( fbusL == areas_buses{m}(k), 1 ) )
            fbus_as{m,k} = [];
        else
            fbus_as{m,k} = find( fbusL == areas_buses{m}(k) );
            fbus_as{m,k} = fbus_as{m,k}(:)';
        end
        
        if isempty(find( tbusL == areas_buses{m}(k), 1 ) )
            tbus_as{m,k} = [];
        else
            tbus_as{m,k} = find( tbusL == areas_buses{m}(k) );
            tbus_as{m,k} = tbus_as{m,k}(:)';
        end
    end
end
ind = 0;
for m = 1:n_areas-1
    for k = m+1:n_areas
        uni = intersect(cell2mat(fbus_as(m,:)),cell2mat(tbus_as(k,:)));
        
        if ~isempty(uni)
            ind = ind +1;
            overlap_s{m,k} = [fbusL(uni); tbusL(uni)];
            overlap_s{k,m} = overlap_s{m,k};
        end
    end
end

ext_area_singlephase = cell(n_areas,1);
for k = 1:n_areas
     ext_area{k} = [];
    for m = 1:length(neigh_area{k})
        ext_area_s{k} = union(overlap_s{k,neigh_area{k}(m)},areas_buses{k});
        ext_area_singlephase{k} = union(ext_area_singlephase{k},ext_area_s{k});
    end
end
if isempty(ext_area_singlephase{1})
	ext_area_singlephase{1} = areas_buses{1};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ext_area = cell(n_areas,1);
for k = 1:n_areas
     ext_area{k} = [];
    for m = 1:length(neigh_area{k})
        ext_area_t{k} = union(overlap{k,neigh_area{k}(m)},areas_buses{k});
        ext_area{k} = union(ext_area{k},ext_area_t{k});
    end
end
if isempty(ext_area{1})
	ext_area{1} = areas_buses{1};
end

for m = 1:n_areas
    [C,IA,IB] = intersect(fbusL,areas_buses{m},'stable');
    [c,ia,ib] = intersect(tbusL,areas_buses{m},'stable');
    
    if length(C) > length(c)
        fbus_a{m} = fbusL(IA);
        tbus_a{m} = tbusL(IA);
        lines{m} = IA;
    else
        fbus_a{m} = fbusL(ia);
        tbus_a{m} = tbusL(ia);
        lines{m} = ia;
    end
end

for m = 1:n_areas
    for k = m+1:n_areas
        [~,overlap_t{k,m},~] = intersect(overlap{m,k},ext_area{m});
		[~,~,overlap_t{m,k}] = intersect(overlap{k,m},ext_area{m});
		
		if isempty(overlap_t{m,k})
			overlap_t{m,k} = [];
			overlap_t{k,m} = [];
		end
	end
end

if ~exist('overlap_t')
	overlap_t = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for m = 1:length(area_codes)
% 	incident_zones{m} = zeros(length(ext_area_singlephase{k}),length(ext_area_singlephase{k}));
% 	for k = 1:nl
% 		incident_zones{m}(fbusL_area{m}(k),tbusL_area{m}(k)) = true;
% 	end
% 	incident_zones{m} = incident_zones{m} + incident_zones{m}';
% end