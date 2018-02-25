function [thetaVec,theta_plusVec]  = AdmRegDiscrete(Theta, thetaMin, thetaMax,delta)
if length(Theta) < 2
    thetaVec = 0;
    theta_plusVec = 0;
    return;
end

thetaVec =      [Theta(1)         Theta(1) min(Theta(1)+delta,thetaMax) Theta(2)       Theta(2) max(Theta(2)-delta,thetaMin)];
theta_plusVec = [min(Theta(1)+delta,thetaMax)  Theta(1) Theta(1)    max(Theta(2)-delta,thetaMin) Theta(2) Theta(2)];
