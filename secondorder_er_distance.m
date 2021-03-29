function [div_found,con_found,chain_found,rec_found,doublerec_found] = secondorder_er_distance(tcluster,tcell,bins,dist_out_prob,dist_in_prob,repeats,central_celltype,partner_celltype,layer)


%preallocate
div_found = zeros(repeats,size(tcluster,1),8);
con_found = zeros(repeats,size(tcluster,1),8);
chain_found = zeros(repeats,size(tcluster,1),8);
rec_found = zeros(repeats,size(tcluster,1),8);
doublerec_found = zeros(repeats,size(tcluster,1),8);

parfor r = 1:repeats
    i_div_found = zeros(size(tcluster,1),8);
    i_con_found = zeros(size(tcluster,1),8);
    i_chain_found = zeros(size(tcluster,1),8);
    i_rec_found = zeros(size(tcluster,1),8);
    i_doublerec_found = zeros(size(tcluster,1),8);
    
    for i = 1:size(tcluster,1)
        
        n1 = tcluster.n_central(i);
        n2 = tcluster.n_partner(i);
        
        cell_central_filter = tcell.Idslice == tcluster.id(i) & tcell.ctype == central_celltype & tcell.layer_sup_deep == layer & ~isnan(tcell.coordinateX);
        cell_partner_filter = tcell.Idslice == tcluster.id(i) & tcell.ctype == partner_celltype & tcell.layer_sup_deep == layer & ~isnan(tcell.coordinateX);
        central_x = tcell.coordinateX(cell_central_filter);
        central_y = tcell.coordinateYorig(cell_central_filter);
        partner_x = tcell.coordinateX(cell_partner_filter);
        partner_y = tcell.coordinateYorig(cell_partner_filter);
        central_xy = [central_x,central_y];
        partner_xy = [partner_x,partner_y];
        all_xy = [central_xy;partner_xy]
        if ~isempty(all_xy)
        [mat] = randneuromat2types_dist(n1,n2,all_xy,bins,dist_out_prob,dist_in_prob)
        
         %outdegree and divergent
        matrix_12filter = zeros(n1+n2);
        matrix_12filter(1:n1,n1+1:n1+n2) = 1; %12 synapses
        [~,divergent,~,~]=secondmotif(mat.*matrix_12filter);
        i_div_found(i,1:n1) = divergent(1:n1);
        
        %indegree and convergent
        matrix_21filter = zeros(n1+n2);
        matrix_21filter(n1+1:n1+n2,1:n1) = 1; %21 synapses
        [convergent,~,~,~]=secondmotif(mat.*matrix_21filter);
        i_con_found(i,1:n1) = convergent(1:n1);
        
        %chain and reciprocal
        [~,~,chain,reciprocal]=secondmotif(mat.*(matrix_12filter + matrix_21filter));
        i_chain_found(i,1:n1) = chain(1:n1);
        
        i_rec_found(i,1:n1) = reciprocal(1:n1);
        
        doublerec = i_rec_found(i,i_rec_found(i,1:n1) >= 2);
        i_doublerec_found(i,1:n1) = sum(doublerec.*(doublerec-1)./2);
        end
    end
    div_found(r,:,:) = i_div_found;
    con_found(r,:,:) = i_con_found;
    chain_found(r,:,:) = i_chain_found;
    rec_found(r,:,:) = i_rec_found;
    doublerec_found(r,:,:) = i_doublerec_found;
end
end