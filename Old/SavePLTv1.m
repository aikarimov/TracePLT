function SavePLTv1(filename,strokes,n,m, canvasW_mm, canvasH_mm, usespp)
fid = fopen(filename, 'w');
fprintf(fid,'IN;');
scl = 40;
W = canvasW_mm*scl; H = canvasH_mm*scl;
sfX = W/n; sfY = H/m;

%масштабируем по меньшему измерению
if sfX < sfY
    sfY = sfX;
end

if sfX > sfY
    sfX = sfY;
end

nStrokes = length(strokes);

col = round(strokes{1,1}.color);

if usespp == 0
    fprintf(fid,'PC%i,%i,%i,%i;',1, col(1), col(2), col(3));
else
    [C,M,Y,B,W] = PCtoPP(col(1), col(2), col(3));
    fprintf(fid,'PP%i,%i,%i,%i,%i;',C,M,Y,B,W);
end

brushSize_mm = 2;
bs = (brushSize_mm*scl);
    
hw = waitbar(0,'Сохранение в файл...');
for i = 1:nStrokes
    waitbar(i/nStrokes,hw, ['Запись в файл/ Число мазков ',num2str(nStrokes)]);
    
    xo = strokes{i,1}.Ys(1)*sfX;
    yo = (m - strokes{i,1}.Xs(1))*sfY;
    
    
    col2 = round(strokes{i,1}.color);
    if ~isequal(col,col2)
        col = col2;
        if usespp == 0
            fprintf(fid,'PC%i,%i,%i,%i;',1, col(1), col(2), col(3));
        else
            [C,M,Y,B,W] = PCtoPP(col(1), col(2), col(3));
            fprintf(fid,'PP%i,%i,%i,%i,%i;',C,M,Y,B,W);
        end
        fprintf(fid,'PU%i,%i;',round(xo),round(yo));
    else
        if(i > 1)
            yprev = (m - strokes{i-1,1}.Xs(end))*sfY;
            xprev = strokes{i-1,1}.Ys(end)*sfX;
            if norm([xo;yo] - [xprev;yprev]) > bs
                fprintf(fid,'PU%i,%i;',round(xo),round(yo)); %PU только если мазки не сливаются в один
            end
        else
            fprintf(fid,'PU%i,%i;',round(xo),round(yo));%или для первого мазка
        end
    end

    for j = 1:length(strokes{i,1}.Xs)
        xo = strokes{i,1}.Ys(j)*sfX;
        yo = (m - strokes{i,1}.Xs(j))*sfY;
        fprintf(fid,'PD%i,%i;',round(xo),round(yo));
    end
end
fclose(fid);
close(hw);