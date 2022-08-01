function newmap = LloydRelax(map, brushSize, Niter, img)
%Lloyd relaxation method for brushstroke distribution
[m, n] = size(img(:,:,1));
bs2 = brushSize/2;

bsD = 2*brushSize;
bsQuad = brushSize^2/4;

%define the size of cluster array
N = length(map);


for it = 1:Niter
    clusters = cell(N,1);
    for k = 1:N %by brushstrokes
        stroke = map{k};
        clusters{k} = cell(length(stroke.Xs),1);
    end
    %make a hash table
    hw = waitbar(0,'Making a hash table...');
    HT = cell(m,n);
    for k = 1:N %by brushstrokes
        waitbar(k/N,hw);
        stroke = map{k};
        minX = realmax; minY = minX;
        maxX = realmin; maxY = maxX;
        
        for l = 1:length(stroke.Xs) %by stroke points
            pX = stroke.Xs(l); minX = min(minX,pX); maxX = max(maxX,pX);
            pY = stroke.Ys(l); minY = min(minY,pY); maxY = max(maxY,pY);
        end
        
        if (minX - bsD < 1)
            minX = 1;
        else
            minX = minX - bsD;
        end
        if maxX + bsD > m
            maxX = m;
        else
            maxX = maxX + bsD ;
        end
        if (minY - bsD < 1)
            minY = 1;
        else
            minY = minY - bsD;
        end
        if maxY + bsD > n
            maxY = n;
        else
            maxY = maxY + bsD ;
        end
        
        for i = minX:maxX
            for j = minY:maxY
                HT{i,j} = [HT{i,j}; k]; %add k to a hash table
            end
        end
    end
    close(hw);
    hw = waitbar(0,'Clustering...');
    
    for i = 1:m
        waitbar(i/n,hw, ['Processing a point (', num2str(i),',1) on iteration ',num2str(it)]);
        for j = 1:n
            %go along canvas dots
            R2 = realmax; %radius^2 to the closest stroke point
            Pid = 1; %index of a point
            hashCell = HT{i,j};
            if length(hashCell) > 0
                Nhash = length(hashCell);
            else
                Nhash = N;
                hashCell = 1:N;
            end
            
            Sid = hashCell(1); %stroke index
            for k = 1:Nhash %by brushstrokes
                stroke = map{hashCell(k)};
                for l = 1:length(stroke.Xs) %by stroke points
                    r = (i - stroke.Xs(l))^2 + (j - stroke.Ys(l))^2; %R^2
                    if( r < R2)
                        R2 = r;
                        Pid = l;
                        Sid = hashCell(k);
                    end
                end
            end
            %for the point i,j the closest stroke point is indexed as
            %s = map{Sid}; pX = s.Xs(Pid), pY = s.Ys(Pid);
            cls = clusters{Sid};
            
            cls{Pid}  = [cls{Pid}; i, j];
            clusters{Sid} = cls;
        end
    end
    
    close(hw);
    %all clusters are ready
    for k = 1:N %by brushstrokes
        stroke = map{k};
        cls = clusters{k};
        for l = 1:length(stroke.Xs) %by stroke points
            XY  = cls{l};
            newXY = round(mean(XY,1));
            if length(newXY > 0)
                stroke.Xs(l) = newXY(1);
                stroke.Ys(l) = newXY(2);
            else
                r = 5;
            end
        end
        map{k} = stroke;
    end
end

%replace colors with respect to img
for i = 1:N %by strokes
    stroke = map{i,1};
    col(1,1,:) = [0, 0, 0];
    colctr = 0;
    L = length(stroke.Xs);
    pX = stroke.Xs(1); pY = stroke.Ys(1);%начало мазка
    for j = 2:L %по точкам мазка
        candidate = struct('X',stroke.Xs(j),'Y',stroke.Ys(j));
        %draw circles
        Nc = max(abs(pX - candidate.X), abs(pY - candidate.Y));
        xo = candidate.X;
        yo = candidate.Y;
        
        for t = 0:1/Nc:1
            
            xp = xo;
            yp = yo;
            
            xo = round(candidate.X  + (pX - candidate.X)*t);
            yo = round(candidate.Y + (pY - candidate.Y)*t);
            
            dx = xo - xp;
            dy = yo - yp;
            
            
            xl = round(xo - bs2); xu = round(xo + bs2);
            yl = round(yo - bs2); yu = round(yo + bs2);
            
            if dx == 0 && dy == 0
                for Xl = xl:xu
                    for Yl = yl:yu
                        goDraw;
                    end
                end
            else
                for Yl = yl:yu
                    if dx > 0
                        xr = round(xu-bs2):xu;
                        xc = xl:round(xu-bs2-1); %complementary
                    else
                        xr = xl:round(xl+bs2); %range
                        xc = round(xl+bs2+1):xu; %complementary
                    end
                    for Xl = xr
                        goDraw;
                    end
                end
                
                for Xl = xc
                    if dy > 0
                        yr = round(yu-bs2):yu;
                    else
                        yr = yl:round(yl+bs2);
                    end
                    for Yl = yr
                        goDraw;
                    end
                end
            end
        end
        
        pX = candidate.X; pY = candidate.Y;
    end
    stroke.color = col/colctr;
    map{i,1} = stroke;
end
    %вложенная процедура
    function goDraw
        if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n
            r2 = (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2; %r^2
            if r2 < bsQuad %core of the line (bs2 - 1)^2
                col = double(img(Xl,Yl,:)) + col;
                colctr = colctr + 1;
            end
        end
    end

newmap = map;

end



