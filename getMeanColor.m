function [col, meancol] = getMeanColor(img, pX, pY, bs2, bsQuad, m,n)
%% GET MEAN COLOR
xo = pX; yo = pY;
colSum = 0;
ncol = 0;
%в квадрате +/- bs2
for Xl = round(xo - bs2):round(xo + bs2)
    for Yl = round(yo - bs2):round(yo + bs2)
        if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n %если квадрат не вылез за рамки холста
            if (double(Xl) - xo)^2 + (double(Yl) - yo)^2 < bsQuad %если расстояние до точки удовлетворяет уравнению круга
                colSum = colSum + double(img(Xl,Yl,:));
                ncol = ncol + 1;
            end
        end
    end
end
col = colSum./ncol; %цвет средний по области
%meancol = mean(col);
meancol = [col(1,1,1), col(1,1,2), col(1,1,3)];
end