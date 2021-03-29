function [mat] = randneuromat2types_dist(n1,n2,all_xy,bins,dist_out_prob,dist_in_prob)

n = n1+n2;
mat=zeros(n);
mat_er=rand(n);

for i=1:n1 %from central (n1)
    for j=n1+1:n %to partner (n2)
        if i ~= j
            dist = pdist(all_xy([i j],:));
            x=find(dist<bins,1)-1; %in which bin is the x difference
            if isempty(x)
                x=size(dist_out_prob,1);
            end
            mat(i,j)=dist_out_prob(x)-mat_er(i,j) > 0;
        end
    end
end

for i=n1+1:n %from partner
    for j=1:n1 %to central
        if i ~= j
            dist = pdist(all_xy([i j],:));
            x=find(dist<bins,1)-1; %in which bin is the x difference
            if isempty(x)
                x=size(dist_in_prob,1);
            end
            mat(i,j)=dist_in_prob(x)-mat_er(i,j) > 0;
        end
    end
end
end