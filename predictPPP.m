% Code by Yuxuan Xia

function [PPP] = predictPPP(PPP,model)

% Predict existing PPP
PPP.w = PPP.w + log(model.Ps);
PPP.GGIW = arrayfun(@(x) predictGGIW(x,model), PPP.GGIW);

% Incorporate PPP birth
PPP.w = [PPP.w;log(model.birth.w)];
PPP.GGIW = [PPP.GGIW;model.birth.GGIW];

end

