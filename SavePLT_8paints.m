function SavePLT_8paints(filename,strokes,n,m, canvasW_mm, canvasH_mm, brushSize_mm)
%Save file for 8 paint machine

if filename == 0
    return;
end

hw = waitbar(0,'Saving a file...');
fid = fopen(filename, 'w');
fprintf(fid,'IN;');
cmdctr = 1;

scl = 40;
W = canvasW_mm*scl; H = canvasH_mm*scl;
sfX = W/n; sfY = H/m;

% Rescale for lower dimension
if sfX < sfY
    sfY = sfX;
end

if sfX > sfY
    sfX = sfY;
end

nStrokes = length(strokes);
col8paints = strokes{1,1}.col8; %color in 8 paints, in motor ticks

cmdctr = cmdctr + 1;
fprintf(fid,'PP'); 
colstr = [];
for i = 1:7 % record first seven paints
    fprintf(fid,'%i,',col8paints(i));
    colstr = [colstr, num2str(col8paints(i)), ' '];
end
fprintf(fid,'%i;',col8paints(8)); %the last paint with semicolon
disp(['PP ',colstr,' command number ',num2str(cmdctr)]);

bs = (brushSize_mm*scl); %new brush size, scaled
    
for i = 1:nStrokes
    waitbar(i/nStrokes,hw, ['Write to file. Number of strokes...',num2str(nStrokes)]);
    
    xo = strokes{i,1}.Ys(1)*sfX;
    yo = (m - strokes{i,1}.Xs(1))*sfY;
    
    col2 = strokes{i,1}.col8;
    if(~isequal(col2,col8paints))
        col8paints = round(strokes{i,1}.col8);
        
        cmdctr = cmdctr + 1;
        fprintf(fid,'PP');
        colstr = [];
        for k = 1:7
            fprintf(fid,'%i,',col8paints(k));
            colstr = [colstr, num2str(col8paints(k)), ' '];
        end
        
        colstr = [colstr, num2str(col8paints(8))];
        
        disp(['PP ',colstr,' command number ',num2str(cmdctr)]);
        
        fprintf(fid,'%i;',col8paints(8));
        
        cmdctr = cmdctr + 1;
        fprintf(fid,'PU%i,%i;',round(xo),round(yo));
    else
        if(i > 1)
            yprev = (m - strokes{i-1,1}.Xs(end))*sfY;
            xprev = strokes{i-1,1}.Ys(end)*sfX;
            if norm([xo;yo] - [xprev;yprev]) > bs
                cmdctr = cmdctr + 1;
                fprintf(fid,'PU%i,%i;',round(xo),round(yo)); % PU only if strokes are merged in one...
            end
        else
            cmdctr = cmdctr + 1;
            fprintf(fid,'PU%i,%i;',round(xo),round(yo));% ... or for the first stroke
        end
    end

    for j = 1:length(strokes{i,1}.Xs)
        xo = strokes{i,1}.Ys(j)*sfX;
        yo = (m - strokes{i,1}.Xs(j))*sfY;
        cmdctr = cmdctr + 1;
        fprintf(fid,'PD%i,%i;',round(xo),round(yo)); %not entirely correct, we may use only first PD before next commands, but this is simpler
    end
end
fclose(fid);
close(hw);