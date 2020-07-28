% Code by Yuxuan Xia

function in_gate = ellipsoidalGating(W,GGIW,model)

if isempty(W)
    in_gate = [];
    return;
end

d = 2;

in_gate = false(size(W,2),1);

H = model.measmodel.H(GGIW.m);
%Take the extent into account when doing ellipsoidal gating
S = GGIW.V/(GGIW.v-2*d-2) + H*GGIW.P*H';
S = (S + S')/2;

nu = W - repmat(model.measmodel.h(GGIW.m),[1,size(W,2)]);
dist= sum((inv(chol(S))'*nu).^2);

%Returns measurement indices inside the gate
in_gate(dist<model.gamma) = true;

end
