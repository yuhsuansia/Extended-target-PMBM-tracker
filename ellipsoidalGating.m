function in_gate = ellipsoidalGating(W,GGIW,model)

d = 2;

in_gate = false(size(W,2),1);

H = model.measmodel.H(GGIW.m);
S = GGIW.V/(GGIW.v-2*d-2) + H*GGIW.P*H' + model.measmodel.R;
S = (S + S')/2;

nu = W - repmat(model.measmodel.h(GGIW.m),[1,size(W,2)]);
dist= sum((inv(chol(S))'*nu).^2);

in_gate(dist<model.gamma) = true;

end
