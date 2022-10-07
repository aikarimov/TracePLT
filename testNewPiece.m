%% TEST NEW PIECE
function [err, candidate_out, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgraymix,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2, bsQuad, mixtype, vwhite)
candidate_out = candidate;
accepted = 0;
z3 = zeros(1,1,3);
if candidate.X <= 0 || candidate.X > m || candidate.Y <= 0 || candidate.Y > n %если вылетели за границы кадра
    %ничего не делаем
    err = realmax;
else
    %если точка с малой погрешностью и не цвета холста на
    %исходнике!
    nX = candidate.X;
    nY = candidate.Y;
    doublePix = double(img(nX,nY,:)); %double image pixel
    errPix = mean(abs(double(col) - doublePix));
    meanerrPix = mean(errPix);

    avrcol = 0;
    if (meanerrPix <= pixTol) &&((canvas2(nX,nY,1) == 0)) && ((meanerrPix > canvasColor + canvasEps) || (meanerrPix < canvasColor - canvasEps)) %рассмотрим новую точку
        ncol = 0;
        colSum = 0;
        overlap = 0;
        %line
        N = max(abs(pX - nX), abs(pY - nY));
        if N == 0
            N = 1;
        end
        for t = 0:1/N:1
            xo = round(nX + (pX - nX)*t);
            yo = round(nY + (pY - nY)*t);
            for Xl = round(xo - bs2):round(xo + bs2)
                for Yl = round(yo - bs2):round(yo + bs2)
                    if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n
                        %test for overlap: number of mixtype is >=, and
                        %amount of white is lower
                        if  ((canvasgraymix(Xl,Yl,1) == mixtype) && (canvasgraymix(Xl,Yl,2) < vwhite)) || (canvasgraymix(Xl,Yl,1) == 0) || (canvasgraymix(Xl,Yl,1) < mixtype)  %if canvas is free
                            r2 = (double(Xl) - xo)^2 + (double(Yl) - yo)^2;
                            if r2 < bsQuad
                                colSum = colSum + double(img(Xl,Yl,:));
                                ncol = ncol + 1;
                                if ~isequal(canvas2(Xl,Yl,:),z3)
                                    overlap = overlap + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        if ncol > 0
            avrcol = ((colSum(1) + colSum(2) + colSum(3))./3.0)/ncol;
            overlap = overlap/ncol;
        else
            overlap = ovf + 1; %заведомо слишком большой, чтобы не принять пиксель
        end
        %если усредненный цвет по области мазка  в
        %границах допустимого цвета и перекрытия
        meancol = mean(meancol);
        if abs(double(avrcol) - double(meancol)) <= pixTol2 && overlap <= ovf
            errCand = abs(double(avrcol) - double(meancol)) + overlap/ncol;
            if(errCand < candidate.err)
                accepted = 1;
                candidate_out.X = nX;
                candidate_out.Y = nY;
                candidate_out.err = errCand;
            end
        end
    end
    err = candidate_out.err;
end
end