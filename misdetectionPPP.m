% Code by Yuxuan Xia

function updatePPP = misdetectionPPP(predictPPP,model)

%Misdetection update of PPP

% For extended target tracking, misdetection can either due to the fact
% that the target is not detected or the fact that the target is detected
% but generates 0 measurements. We should take both alternatives into
% account.

w1 = predictPPP.w + log(model.Qd);
w2 = predictPPP.w + log(model.Pd) + arrayfun(@(x) x.a*log(x.b/(x.b+1)), predictPPP.GGIW);
updatePPP.w = [w1;w2];

updatePPP.GGIW = [predictPPP.GGIW;arrayfun(@(x) bplus(x), predictPPP.GGIW)];

    function GGIW = bplus(GGIW)
        GGIW.b = GGIW.b + 1;
    end

end

