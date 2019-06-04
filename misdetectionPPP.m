function updatePPP = misdetectionPPP(predictPPP,model)

w1 = predictPPP.w + log(model.Qd);
w2 = predictPPP.w + log(model.Pd) + arrayfun(@(x) x.a*log(x.b/(x.b+1)), predictPPP.GGIW);
updatePPP.w = [w1;w2];

updatePPP.GGIW = [predictPPP.GGIW;arrayfun(@(x) bplus(x), predictPPP.GGIW)];

    function GGIW = bplus(GGIW)
        GGIW.b = GGIW.b + 1;
    end

end

