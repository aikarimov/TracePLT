function [err, canvas_out, canvasgray_out, canvas2_out]= drawPiece(pX,pY,candidate,bs2,bsQuad,canvas, canvasgraymix,canvas2,col,col2,imggray,m,n, mixtype, vwhite)
%% DRAW PIECE
canvas_out = canvas;
canvasgray_out = canvasgraymix;
canvas2_out = canvas2;

N = max(abs(pX - candidate.X), abs(pY - candidate.Y));

for t = 0:1/N:1
    xo = round(candidate.X  + (pX - candidate.X)*t);
    yo = round(candidate.Y + (pY - candidate.Y)*t);
    for Xl = round(xo - bs2):round(xo + bs2)
        for Yl = round(yo - bs2):round(yo + bs2)
            if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n
                if (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2 < bsQuad
                    %test for overlap: number of mixtype is >=, and
                    %amount of white is lower
                    if  ((canvasgraymix(Xl,Yl,1) == mixtype) && (canvasgraymix(Xl,Yl,2) < vwhite)) || (canvasgraymix(Xl,Yl,1) == 0) || (canvasgraymix(Xl,Yl,1) < mixtype)  %if canvas is free
                        canvas_out(Xl,Yl,:) = col;
                        canvasgray_out(Xl,Yl,1) = mixtype;
                        canvasgray_out(Xl,Yl,2) = vwhite;
        
                        canvas2_out(Xl,Yl,:) = col2; %new variant
                    end
                end
%                 if (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2 < bsQuad12
%                     canvas2_out(Xl,Yl,:) = col2;
%                 end
            end
        end
    end
end
meancnv = (canvas_out(:,:,1) + canvas_out(:,:,2) + canvas_out(:,:,3))/3;
err = abs(meancnv - double(imggray));
end