%“рассировщик PLT 
%% TRACE - MAIN

function TracePLT
close all
canvasColor = 255;
TotalIter = 2;
do_blur = 0;

[filename,pathname] = uigetfile({'*.png;*.jpg;*.bmp','Images (*.png;*.jpg;*.bmp)';'*.*','All files (*.*)'},'Select the file');
img = imread([pathname,filename]); %img contains image
 
%img = imread('C:\\Users\\Sapr1152\\Dropbox\\ARTCYBE\\∆ивописец\\MATLAB\\Lenna256.jpg');

[m, n] = size(img(:,:,1));

[brushWidth_mm, canvasW_mm, canvasH_mm, ~, nColors] = GUI_Trace;
canvasEps = 2; % +\- погрешность задани€ цвета холста


figure(1);
imshow(img);

canvasgray =  (canvasColor + zeros(m,n,1)); %gray version of canvas with tones
canvas = (canvasColor + zeros(m,n,3)); %colored version

figure(2);
imshow(canvas);
err = abs(double(canvas) - double(img)); %белый это цвет фона

strokes = cell(0);
itersMinOverlap = 1; %итераций с малым перекрытием

minOverlap = 0.6; %макс. коэффициент перекрыти€ начальный
% maxOverlap = 0.8;
maxOverlap = 1.0;

pixTol = 6; %возможное отклонение цвета от исходника на конце
pixTol2 = 100; %в среднем
pixTolBest = 4; %погрешность, при которой мазок безоговорочно прин€т

sfX = canvasW_mm/n; %mm->pix
sfY = canvasH_mm/m;

%масштабируем по меньшему измерению
if sfX < sfY
    sfY = sfX;
end
if sfX > sfY
    sfX = sfY;
end

brushSize = round(brushWidth_mm/sfX)%толщина кисти
bs2 = ceil(brushSize/2);
MaxIters = brushSize;%макс. число попыток найти новый отрезок

%размытие в соответствии с кистью - —ƒ≈Ћј“№ ќѕ÷»ќЌјЋ№Ќќ!

img_unblurred = img;
if do_blur
    H = fspecial('disk',bs2);
    img = imfilter(img,H,'replicate');
    figure(19);
    imshow(img); title('–азмытый');
end


%найдем градиенты
imggray = rgb2gray(img);
%cамописна€ функци€
[U, V] = GetGradientByTensor(imggray, brushSize);

%VISUALIZATION

figure(4); 
quiver(U,V);
axis ij
axis equal

Rmax = brushSize;%sqrt(m*n)*0.03; %макс длина 1 отрезка
Rmin = 1;%sqrt(m*n)*0.005;
maxLen = brushSize*10; %макс. длина мазка
nStrokes = 0;
accepted = 0;

bsQuad = (double(brushSize)/2)^2;
hw = waitbar(0,'wait...');
canvas2 = zeros(m,n,3); %for brushstroke images 
%%

for kk = 1:TotalIter

    if kk == 2 %на второй итерации замен€ем изображение на неразмытое
        img = img_unblurred;
    end
    
    i_loop = 1:m;
    j_loop = 1:n;
    if (kk > itersMinOverlap)
        ovf = maxOverlap;
    else
        ovf = minOverlap;
    end
    for i = i_loop
        for j = j_loop
            if (err(i,j) > pixTol) || canvas2(i,j) == 0
                %если в нем высока€ погрешность
                pX = i; pY = j; %prev X,Y
                [col, meancol] = getMeanColor(img, pX, pY, bs2, bsQuad, m, n); %average color on the area
                
                col2 = 0.99*[rand rand rand] + 0.01; %RGB color for stroke map
                
                %stroke points
                strokeStructure = struct('color',col,'Xs',[pX],'Ys',[pY],'length',0,'col8',zeros(1,8)); %added new field col8 for paints
                nPts = 1; %number of points in the stroke
                endedStroke = 0;
                while endedStroke == 0 %until the stroke is not ended
                    %find new direction
                    candidate = struct('X',-1,'Y',-1,'err',realmax);
                    ctr = 1;
                    while  ctr < MaxIters
                        [cosA, sinA] = getDirection(pX,pY, strokeStructure,U,V);
                        r = MaxIters - ctr;

                        candidate.X = pX + round(r*cosA); %new X
                        candidate.Y = pY + round(r*sinA); %new Y
                        
                        %test new fragment of the stroke
                        [~, candidate, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgray,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2,bsQuad, canvas);
                        
                        if accepted %if error is small, accept the stroke immediately
                            ctr = MaxIters;
                        else
                            %also try opposite direction
                            candidate.X = pX - round(r*cosA); %new X
                            candidate.Y = pY - round(r*sinA); %new Y
                            
                            %test new fragment of the stroke
                            [~, candidate, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgray,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2,bsQuad, canvas);
                            if accepted %if error is small, accept the stroke immediately
                                ctr = MaxIters;
                            end
                        end
                        
                        ctr = ctr + 1;
                    end
                    
                    %draw the stroke fragment
                    
                    if (accepted)
                        %new message, update message bar
                        waitbar((i*n+j)/m/n,hw, ['iter = ', num2str(kk), 'i,j = ', num2str(i),',',num2str(j),' pxpy = ',num2str(candidate.X),',',num2str(candidate.Y),' nStrokes = ', num2str(nStrokes)]);
                        
                        [err, canvas, canvasgray, canvas2] = drawPiece(pX,pY,candidate,bs2,bsQuad,canvas,canvasColor,canvasgray,canvas2,meancol,col,col2,imggray, m,n);
                        %[err, canvas, canvasgray, canvas2] = drawPieceAntialiaced(pX,pY,candidate,bs2,bsQuad,canvas,canvasColor,canvasgray,canvas2,meancol,col,col2,imggray, m,n);
                        
                        
                        %determine new length
                        dlen = sqrt((pX - candidate.X)^2 +  (pY - candidate.Y)^2);
                        strokeStructure.length = strokeStructure.length + dlen;
                        
                        %if length is too large
                        if(strokeStructure.length >= maxLen)
                            endedStroke = 1;
                        end
                        
                        vX = candidate.X; vY = candidate.Y;
                        pX = candidate.X; pY = candidate.Y;

                        strokeStructure.Xs = [strokeStructure.Xs, vX];
                        strokeStructure.Ys = [strokeStructure.Ys, vY];
                        nPts = nPts + 1;
                    else
                        endedStroke = 1;
                    end
                end
                if nPts > 1
                    nStrokes = nStrokes + 1;
                    
                    %mutating color - optional
                    %strokeStructure.color = mutateColor(strokeStructure.color);
                    
                    strokes{nStrokes,1} = strokeStructure;
                end
            end
        end
    end
    figure(2);
    imshow(uint8(canvas));title('Canvas');
    
    figure(3);
    imshow(canvas2);title('Brushstroke map');
    
    figure(5);
    [~,ch] =  contourf(err);
    title('Error');
    set(ch,'edgecolor','none');
    ax = gca;
    ax.YDir = 'reverse';
    colorbar
    axis equal
end

colArray = zeros(nStrokes,3);
wArray = zeros(nStrokes,1);

for i = 1:nStrokes
    colArray(i,:) = strokes{i,1}.color; %color array
    wArray(i) = strokes{i,1}.col8(5); %white colors
end

if nColors > nStrokes %in order to avoid error in k-means
    nColors = nStrokes;
end
ids = kmeans(double(colArray),nColors); %clustering the colors, reducing the number of colors

colorsFinal = zeros(nColors,3);
nels = zeros(1,nColors); % number of color array elements
for i = 1:nStrokes %for each stroke, rewrite the color
    dcol = double(strokes{i,1}.color); dcol = [dcol(1,1,1) dcol(1,1,2) dcol(1,1,3)];
    colorsFinal(ids(i),:) = colorsFinal(ids(i),:) + dcol;
    nels(ids(i)) = nels(ids(i)) + 1;
end

for j = 1:nColors
    colorsFinal(j,:) = double(colorsFinal(j,:))./nels(j);
end



%now, transform colors into real proportions
[props, cls, hsvnew] = PredictProportions(rgb2hsv(colorsFinal/255)); %in a format double from 0 to 1

rgbnew = uint8(255*hsv2rgb(hsvnew));

mixvalues = prop2pp(props,cls);

%create mix groups
valuecell = cell(1,3);
for j = 1:nColors
    %in the class cls(j) change 
    valuecell{cls(j)} = mixvalues(j,:);
end

%sort by amount of white color in the overall mix

%[~, ind] = sort(mean(colArray,2),'descend');
[~, ind] = sort(mixvalues(:,5),'ascend');
strokes2 = zeros(nColors,nStrokes); %array of stroke indices
strCount = zeros(nColors,1); %number of strokes of each color

%now, order colorsFinal по количеству белого
%*************************************************************************
tints = wArray;
[~, idsCol] = sort(tints,'ascend');
%*************************************************************************

map2 = cell(nStrokes,1);

for i = 1:nStrokes
    i_cluster = ids(ind(i));
    i_sorted = ind(i);
    strokes{i_sorted,1}.color = uint8(colorsFinal(i_cluster,:));
    strokes{i_sorted,1}.col8 = mixvalues(i_cluster,:);
    strCount(i_cluster) = strCount(i_cluster) + 1;
    strokes2(i_cluster,strCount(i_cluster)) = i_sorted;
    map2{i,1} = strokes{i_sorted,1}; %map with sorted colors
end

%sort strokes by position in each cluster, to reduce the machine tool path
strokes3 = strokes2;
for i = 1:nColors % for each cluster
    waitbar(i/nColors,hw, ['Sorting strokes by position, i = ', num2str(i)]);
    curStrokeId = strokes2(i,1); %index in strokes
    strokes2(i,1) = 0;%0, because the i-th stroke is already taken
    j = 2;
    while j < strCount(i) %for the strokes in a cluster
        %найдем пару дл€ curStroke
        curStroke = strokes{curStrokeId,1};
        curX = curStroke.Xs(end);
        curY = curStroke.Ys(end);
        distMin = realmax;
        nextStrokeId = curStrokeId;
        k = 1;
        nextNumber = 0;
        while k < strCount(i)
            if strokes2(i,k) ~= 0
                candidate = strokes{strokes2(i,k),1};
                dist = sqrt( (curX - candidate.Xs(1))^2 + (curY - candidate.Ys(1))^2 );
                if (dist < distMin)
                    distMin = dist;
                    nextStrokeId = strokes2(i,k);
                    nextNumber = k;
                end
            end
            k = k + 1;
        end
        strokes3(i,j) = nextStrokeId;
        curStrokeId = nextStrokeId;
        strokes2(i,nextNumber) = 0; %processed
        j = j+1;
    end
end

%put them into the array again
map = cell(nStrokes,1);
ctr = 1;
for i = 1:nColors
    iCol = idsCol(i); %correct order of clusters
    for j = 1:strCount(iCol)
        map{ctr,1} = strokes{strokes3(iCol,j),1};
        ctr = ctr + 1;
    end
end
close(hw);

figure(13);
imshow(map2imgColorCanvas(brushSize,canvas,map,canvasColor));
title('Final image');

%% SAVING
filenameforsave = uiputfile([filename(1:end-4),'.plt'],'Save file as...');
SavePLT_8paints(filenameforsave, map,n,m,canvasW_mm,canvasH_mm, brushWidth_mm);

end

%% GET MEAN COLOR
function [col, meancol] = getMeanColor(img, pX, pY, bs2, bsQuad, m,n)
xo = pX; yo = pY;
colSum = 0;
ncol = 0;
%в квадрате +/- bs2
for Xl = round(xo - bs2):round(xo + bs2)
    for Yl = round(yo - bs2):round(yo + bs2)
        if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n %если квадрат не вылез за рамки холста
            if (double(Xl) - xo)^2 + (double(Yl) - yo)^2 < bsQuad %если рассто€ние до точки удовлетвор€ет уравнению круга
                colSum = colSum + double(img(Xl,Yl,:));
                ncol = ncol + 1;
            end
        end
    end
end
col = colSum./ncol; %цвет средний по области
meancol = mean(col);
end

%% GET DIRECTION
function [cosA, sinA] = getDirection(pX,pY, strokePoints, U, V)

cosA = -U(pX,pY);
sinA = V(pX,pY); %normally to gradient

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
end

%% TEST NEW PIECE
function [err, candidate_out, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgray,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2, bsQuad, canvas)
candidate_out = candidate;
accepted = 0;
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
    avrcol = 0;
    %if (errPix <= pixTol) &&(mean(double(canvas(nX,nY,:)) - double(img(nX,nY,:))) > pixTol2 || (canvas2(nX,nY,1) == 0)) && ((mean(doublePix) > canvasColor + canvasEps) || (mean(doublePix) < canvasColor - canvasEps)) %рассмотрим новую точку
    if (errPix <= pixTol) &&((canvas2(nX,nY,1) == 0)) && ((mean(doublePix) > canvasColor + canvasEps) || (mean(doublePix) < canvasColor - canvasEps)) %рассмотрим новую точку
      
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
                        %если цвет пикс. холста темнее (перекрытие светлых)
                        if  (canvasgray(Xl,Yl) < meancol) || (canvasgray(Xl,Yl) == canvasColor)  %если цвет пикс. холста светлее
                            r2 = (double(Xl) - xo)^2 + (double(Yl) - yo)^2;
                            if r2 < bsQuad
                                colSum = colSum + double(img(Xl,Yl,:));
                                ncol = ncol + 1;
                                if(canvas2(Xl,Yl) ~= 0)
                                    overlap = overlap + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        if ncol > 0
            avrcol = (colSum(1) + colSum(2) + colSum(3))./ncol/3;
            overlap = overlap/ncol;
        else
            overlap = ovf + 1; %заведомо слишком большой, чтобы не прин€ть пиксель
        end
        %если усредненный цвет по области мазка  в
        %границах допустимого цвета и перекрыти€
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

%% DRAW PIECE
function [err, canvas_out, canvasgray_out, canvas2_out]= drawPiece(pX,pY,candidate,bs2,bsQuad,canvas,canvasColor, canvasgray,canvas2,meancol,col,col2,imggray,m,n)
canvas_out = canvas;
canvasgray_out = canvasgray;
canvas2_out = canvas2;

N = max(abs(pX - candidate.X), abs(pY - candidate.Y));
bsQuad12 = (double(bs2)/2)^2;
for t = 0:1/N:1
    xo = round(candidate.X  + (pX - candidate.X)*t);
    yo = round(candidate.Y + (pY - candidate.Y)*t);
    for Xl = round(xo - bs2):round(xo + bs2)
        for Yl = round(yo - bs2):round(yo + bs2)
            if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n
                if (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2 < bsQuad
                    if  (canvasgray(Xl,Yl) < meancol) || (canvasgray(Xl,Yl) == canvasColor)
                        canvas_out(Xl,Yl,:) = col;         %ибо темные после светлых
                        canvasgray_out(Xl,Yl) = meancol;
                    end
                end
                if (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2 < bsQuad12
                    canvas2_out(Xl,Yl,:) = col2;
                end
            end
        end
    end
end
err = (abs(double(canvasgray) - double(imggray)));
end


%% DRAW ANTI-ALIACED PIECE
% function [err, canvas_out, canvasgray_out, canvas2_out]= drawPieceAntialiaced(pX,pY,candidate,bs2,bsQuad,canvas,canvasColor, canvasgray,canvas2,meancol,col,col2,imggray,m,n)
% canvas_out = canvas;
% canvasgray_out = canvasgray;
% canvas2_out = canvas2;
% bsQuad12 = (double(bs2)/2)^2;
% 
% bsQuad_1 = (bs2 - 1)^2;
% bsQuadp1 = (bs2 + 1)^2;
% 
% %draw circles
% N = max(abs(pX - candidate.X), abs(pY - candidate.Y));
% xo = candidate.X;
% yo = candidate.Y;
% 
% for t = 0:1/N:1
%     
%     xp = xo;
%     yp = yo;
%     
%     xo = round(candidate.X  + (pX - candidate.X)*t);
%     yo = round(candidate.Y + (pY - candidate.Y)*t);
%     
%     dx = xo - xp;
%     dy = yo - yp;
%     
%    
%     xl = round(xo - bs2); xu = round(xo + bs2);
%     yl = round(yo - bs2); yu = round(yo + bs2);
%     
%     if dx == 0 && dy == 0
%         for Xl = xl:xu
%             for Yl = yl:yu
%                 goDraw;
%             end
%         end
%     else
%         for Yl = yl:yu
%             if dx > 0
%                 xr = xu-bs2:xu;
%                 xc = xl:xu-bs2-1; %complementary
%             else
%                 xr = xl:xl+bs2; %range
%                 xc = xl+bs2+1:xu; %complementary
%             end
%             for Xl = xr
%                 goDraw;
%             end
%         end
%         
%         for Xl = xc
%             if dy > 0
%                 yr = yu-bs2:yu;
%             else
%                 yr = yl:yl+bs2;
%             end
%             for Yl = yr
%                 goDraw;
%             end
%         end
%     end
%     
%     
% end
% err = abs(mean(double(canvas_out),3) - double(imggray));
% 
%     %вложенна€ процедура
%     function goDraw
%         if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n
%             r2 = (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2; %r^2
%             
%             if r2 < bsQuad_1 %core of the line (bs2 - 1)^2
%                 if  (canvasgray(Xl,Yl) < meancol) || (canvas2(Xl,Yl) == 0) %non-painted
%                     canvas_out(Xl,Yl,:) = col;         %ибо темные после светлых
%                     canvasgray_out(Xl,Yl) = meancol;
%                 end
%             else
%                 if r2 < bsQuadp1 % (bs2 + 1)^2
%                     gamma = (bs2 - sqrt(r2))/2; % approximately equals (r - bs2)
%                     if gamma > 1
%                         gamma = 1;
%                     end
%                     if gamma < 0
%                         gamma = 0;%%%%
%                     end
%                     if  (canvasgray(Xl,Yl) < meancol) || (canvas2(Xl,Yl) == 0) %non-painted
%                         canvas_out(Xl,Yl,:) = uint8(canvas_out(Xl,Yl,:)*(1 - gamma)) + uint8(gamma*col);         %ибо темные после светлых
%                         canvasgray_out(Xl,Yl) = uint8(canvasgray_out(Xl,Yl)*(1 - gamma)) + uint8(gamma*meancol); % без антиэлайсинга
%                     end
%                 end
%             end
%             if r2 < bsQuad12 %без анти-элайсинга
%                 canvas2_out(Xl,Yl,:) = col2;
%             end
%         end
%     end
% 
% end

%% MUTATE COLOR
% function color_mutated = mutateColor(col)
% color_mutated = col;
% for i = 1:3
%     newcol = col(1,1,i).*(0.95 + 0.1*rand);
%     if newcol > 255 
%         newcol = 255;
%     end
%     color_mutated(1,1,i) = newcol;
% end
% end


