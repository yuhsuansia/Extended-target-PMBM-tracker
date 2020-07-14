% Code by Yuxuan Xia

function [what,GGIWhat] = mixtureReduction(w,GGIW,model)

% Use a greedy merging method to reduce the number of mixture components

if length(w) == 1
    what = w;
    GGIWhat = GGIW;
    return;
end

%Index set of components
I = 1:length(GGIW);
el = 1;

while ~isempty(I)
    Ij = [];
    %Find the component with the highest weight
    [~,j] = max(w);
    
    for i = I
        %Find other similar GGIW components in the sense of KL divergence
        if i == j
            Ij= [ Ij i ];
        elseif GGIW_KLdiv(GGIW(j),GGIW(i)) <= model.merge
            Ij= [ Ij i ];
        end
    end
    
    %Merge components by moment matching
    [temp,what(el,1)] = normalizeLogWeights(w(Ij));
    [~,GGIWhat(el,1)] = GGIW_merge_wrap(temp,GGIW(Ij));
    
    %Remove indices of merged components from index set
    I = setdiff(I,Ij);
    %Set a negative to make sure this component won't be selected again
    w(Ij,1) = log(eps);
    el = el+1;
end

end