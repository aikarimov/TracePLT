function img2 = map2imgColorCanvas(brushSize,canvas, map, canvasColor)
[m, n] = size(canvas(:,:,1));
img2 = uint8(canvasColor + zeros(m,n,3));
bsQuad = brushSize^2/4;
for i = 1:length(map) %по мазкам
    stroke = map{i,1};
    col = stroke.color;
    pX = stroke.Xs(1); pY = stroke.Ys(1);%начало мазка
    for j = 2:length(stroke.Xs) %по точкам мазка
        candidate = struct('X',stroke.Xs(j),'Y',stroke.Ys(j));
        N = max(abs(pX - candidate.X), abs(pY - candidate.Y));
        for t = 0:1/N:1
            xo = round(candidate.X  + (pX - candidate.X)*t);
            yo = round(candidate.Y + (pY - candidate.Y)*t);
            for Xl = round(xo - brushSize/2):round(xo + brushSize/2)
                for Yl = round(yo - brushSize/2):round(yo + brushSize/2)
                    if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n
                        if (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2 < bsQuad
                            img2(Xl,Yl,:) = col;
                        end
                    end
                end
            end
        end
        pX = candidate.X; pY = candidate.Y;
    end
end
end
