function [cosA, sinA] = getDirection(pX,pY, strokePoints, U, V, go_normal)

if go_normal
    cosA = -U(pX,pY);
    sinA = V(pX,pY); %normally to gradient
else
    sinA = U(pX,pY);
    cosA = V(pX,pY); %collinearly to gradient
end



if(length(strokePoints.Xs) > 1) %if not a single point, get previous direction
    dX = pX - strokePoints.Xs(end-1);
    dY = pY - strokePoints.Ys(end-1);
    
    %get scalar product
    sc = cosA*dX + sinA*dY;
    if sc < 0 %if in opposite direction
        cosA = -cosA;
        sinA = -sinA;
    end
end


%small perturbation
% dX = 0; dY = 0;
% pf = pi/12; %perturbation factor, rad
% alp = atan2(sinA,cosA);
% if dX ~= 0 && dY ~= 0
%     bet = atan2(dY,dX);
%     pf = 1.2*(bet - alp);
% end
% 
% da = 2*pf*rand - pf;
% cosA = cos(alp + da);
% sinA = sin(alp + da);

end