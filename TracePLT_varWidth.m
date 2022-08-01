%Трассировщик PLT 
%% TRACE - MAIN WITH VARIABLE WIDTH
function TracePLT_varWidth
close all
canvasColor = 255;
TotalIter = 1;
do_blur = 0;

[filename,pathname] = uigetfile({'*.png;*.jpg;*.bmp','Images (*.png;*.jpg;*.bmp)';'*.*','All files (*.*)'},'Select the file');
img = imread([pathname,filename]); %img contains image
 
%img = imread('C:\\Users\\Sapr1152\\Dropbox\\ARTCYBE\\Живописец\\MATLAB\\Lenna256.jpg');

[m, n] = size(img(:,:,1));

[brushWidth_mm, canvasW_mm, canvasH_mm, ~, nColors] = GUI_Trace;
canvasEps = 2; % +\- погрешность задания цвета холста

%for predicting mix color values
tablename = 'ModelTable600.xls';
nmixgroups = 4;
Ycell = cell(1,nmixgroups); % 0 1 2 3 = MY1 MY2 CY CM
Wcell = cell(1,nmixgroups);

for i = 1:nmixgroups
    M1 = readmatrix(tablename,'Sheet',i);
    Ycell{i} = M1(:,1:3);
    Wcell{i} = M1(:,4:6);
end
%read Ycell and Wcell and remember

figure(1);
imshow(img);

canvasgraymix =  zeros(m,n,2); %gray version of canvas with tones
                               %first layer is color class
                               %second is volume of white
canvas = (canvasColor + zeros(m,n,3)); %colored version

figure(2);
imshow(canvas);
err = abs(double(canvas) - double(img)); %белый это цвет фона

strokes = cell(0);


sfX = canvasW_mm/n; %mm->pix
sfY = canvasH_mm/m;

%масштабируем по меньшему измерению
if sfX < sfY
    sfY = sfX;
end
if sfX > sfY
    sfX = sfY;
end

brushSize = round(brushWidth_mm/sfX) %толщина кисти - минимальная
bs2 = ceil(brushSize/2);
MaxIters = brushSize;%макс. число попыток найти новый отрезок

itersMinOverlap = 1; %итераций с малым перекрытием

minOverlap = 0.6; %макс. коэффициент перекрытия начальный
% maxOverlap = 0.8;
maxOverlap = 1.0;

pixTol = 6; % 6; %возможное отклонение цвета от исходника на конце
pixTol2 = 100; %в среднем
pixTolBest = 4; %погрешность, при которой мазок безоговорочно принят

maxLen = brushSize*10; %макс. длина мазка

%размытие в соответствии с кистью - СДЕЛАТЬ ОПЦИОНАЛЬНО!

img_unblurred = img;
if do_blur
    H = fspecial('disk',bs2);
    img = imfilter(img,H,'replicate');
    figure(19);
    imshow(img); title('Размытый');
end


%найдем градиенты
imggray = rgb2gray(img);
%cамописная функция
[U, V] = GetGradientByTensor(imggray, brushSize);

%VISUALIZATION

figure(4); 
quiver(U,V);
axis ij
axis equal

Rmax = brushSize;%sqrt(m*n)*0.03; %макс длина 1 отрезка
Rmin = 1;%sqrt(m*n)*0.005;

nStrokes = 0;
accepted = 0;

bsQuad = (double(brushSize)/2)^2;
hw = waitbar(0,'wait...');
canvas2 = zeros(m,n,3); %for brushstroke images 
%%

for kk = 1:TotalIter

    if kk == 2 %на второй итерации заменяем изображение на неразмытое
        img = img_unblurred;
    end
    
    i_loop = 1:m;
    j_loop = 1:n;
    if (kk > itersMinOverlap)
        ovf = maxOverlap;
    else
        ovf = minOverlap;
    end
    j = 1;
    for i = i_loop
        %new message, update message bar
        waitbar((i*n+j)/m/n,hw, ['iter = ', num2str(kk), 'i,j = ', num2str(i),',',num2str(j),' nStrokes = ', num2str(nStrokes)]);
                       
        for j = j_loop
            if (err(i,j) > pixTol) || canvas2(i,j) == 0
                %если в нем высокая погрешность
                pX = i; pY = j; %prev X,Y
                [col, meancol] = getMeanColor(img, pX, pY, bs2, bsQuad, m, n); %average color on the area

                [props, mixtype, realhsv]  = PredictProportions(rgb2hsv(double(meancol)/255), Ycell, Wcell);
                col8paints = prop2pp(props,mixtype);
                vwhite = col8paints(5);

                %col(1,1,:) = uint8(255*hsv2rgb(realhsv)); %%%%%%%%%%%%%
                
                %col2 = 0.99*[rand rand rand] + 0.01; %RGB color for stroke map

                switch mixtype
                    case 1
                        col2 = [(0.5 + 0.5*rand) 0 0] + 0.5*[0 rand rand];
                    case 2
                        col2 = [(0.5 + 0.5*rand) (0.5 + 0.5*rand) 0] + 0.5*[0 0 rand];
                    case 3
                        col2 = [0 (0.5 + 0.5*rand) 0] + 0.5*[rand 0 rand];
                    case 4
                        col2 = [0 0 (0.5 + 0.5*rand)] + 0.5*[0 0 rand];
                end
                
                %stroke points
                strokeStructure = struct('color',col,'Xs',[pX],'Ys',[pY],'length',0,'col8',col8paints,'props',props,'cls',mixtype); %added new field col8 for paints
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
                        [~, candidate, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgraymix,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2,bsQuad, mixtype, vwhite);
                        
                        if accepted %if error is small, accept the stroke immediately
                            ctr = MaxIters;
                        else
                            %also try opposite direction
                            candidate.X = pX - round(r*cosA); %new X
                            candidate.Y = pY - round(r*sinA); %new Y
                            
                            %test new fragment of the stroke
                            [~, candidate, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgraymix,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2,bsQuad, mixtype, vwhite);
                            if accepted %if error is small, accept the stroke immediately
                                ctr = MaxIters;
                            end
                        end
                        
                        ctr = ctr + 1;
                    end
                    
                    %draw the stroke fragment
                    
                    if (accepted)
                        [err, canvas, canvasgraymix, canvas2] = drawPiece(pX,pY,candidate,bs2,bsQuad,canvas,canvasColor,canvasgraymix,canvas2,meancol,col,col2,imggray, m,n, mixtype, vwhite);
                               
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
                    %make color real
                    %strokeStructure.color = 255*hsv2rgb(realhsv);

                    dcol = strokeStructure.color;
                    strokeStructure.color = [dcol(1,1,1) dcol(1,1,2) dcol(1,1,3)];

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

%%%%%%%%%%%%%%%%%%%

%create array for each mix type
mixCell = cell(1,nmixgroups); % [props, cls]
idcell = cell(1,nmixgroups); 
for i = 1:nStrokes
    mixtype = strokes{i,1}.cls;
    mixCell{mixtype} = [mixCell{mixtype}; strokes{i,1}.props]; % color array
    idcell{mixtype} =  [ idcell{mixtype}; i]; %real order of stroke in cell array
end

%now, rescale K for each mix group and perform internal k-means among
%proportions
nColors2 = 0;

ids = [];
cls = [];
mixvalues = [];
colorsFinal = [];

i2i = [];

ctr = 0;
for i = 1:nmixgroups
    idx = idcell{i};
    i2i = [i2i; idx];
    mixArray = mixCell{i};
    Ncl = size(mixArray,1);

    if Ncl > 0
        Ki = ceil(Ncl / nStrokes * nColors);


        if Ki >= Ncl
            Ki = Ncl;
            idsk = (1:Ncl)'; %indices directly transform into clusters
        else
            idsk = kmeans(mixArray,Ki); %clustering the proportions, reducing the number of colors
        end

        nColors2 = nColors2 + Ki;

        propsFinal = zeros(Ki,3);
        nels = zeros(1,Ki);
        for j = 1:Ncl
            dprops = mixArray(j,:);
            propsFinal(idsk(j),:) = propsFinal(idsk(j),:) + dprops; %sum of colors before averaging
            nels(idsk(j)) = nels(idsk(j)) + 1;
        end

        %averaging proportions in each cluster
        for j = 1:Ki
            propsFinal(j,:) = double(propsFinal(j,:))./nels(j);
        end

        %transform proportions into colors
        hsvcolors = prop2hsv(propsFinal, i*ones(Ncl,1), Wcell,Ycell);
        rgbnew = uint8(255*hsv2rgb(hsvcolors));
        colorsFinali = rgbnew; %rewrite final colors
        mixvaluesi = prop2pp(propsFinal,i*ones(Ncl,1));

        ids = [ids; idsk + ctr]; %global indices
        cls = [cls; i*ones(Ki,1)]; %classes of colors
        mixvalues = [mixvalues; mixvaluesi];
        colorsFinal = [colorsFinal; colorsFinali];

        ctr = ctr + Ki;
    end
end

nColors = nColors2; %if round-off error meakes nColors change, take it into account

%re-order ids
idsnew = ids;
for i = 1:nStrokes
   idsnew(i2i(i)) = ids(i);
end

ids = idsnew;
%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%



% %create color array for colors
% colArray = zeros(nStrokes,3);
% for i = 1:nStrokes
%     colArray(i,:) = strokes{i,1}.color; % color array
% end
% if nColors > nStrokes %in order to avoid error in k-means
%     nColors = nStrokes;
% end
% ids = kmeans(double(colArray),nColors); %clustering the colors, reducing the number of colors
% colorsFinal = zeros(nColors,3);
% nels = zeros(1,nColors); % number of color array elements
% 
% for i = 1:nStrokes %for each stroke, rewrite the color
%     dcol = double(strokes{i,1}.color); 
%     colorsFinal(ids(i),:) = colorsFinal(ids(i),:) + dcol; %sum of colors before averaging
%     nels(ids(i)) = nels(ids(i)) + 1;
% end
% %averaging for obtaining average rgb color among each cluster
% for j = 1:nColors
%     colorsFinal(j,:) = double(colorsFinal(j,:))./nels(j);
% end
% 
% %now, transform colors into real proportions
% [props, cls, hsvnew] = PredictProportions(rgb2hsv(double(colorsFinal)/255),Ycell,Wcell); %in a format double from 0 to 1
% 
% rgbnew = uint8(255*hsv2rgb(hsvnew));
% colorsFinal = rgbnew; %rewrite final colors 'as is possible'
% mixvalues = prop2pp(props,cls);


%%%%%%%%%%%%%%%%%%%





%create mix groups
iclscell = cell(1,nmixgroups);
for j = 1:nmixgroups
    iclscell{j} = [];
end

%re-order clusters
for j = 1:nColors
    %in the class cls(j) change 
    mxarray = iclscell{cls(j)};
    mxarray = [mxarray; j]; %append new value of color mix
    iclscell{cls(j)} = mxarray;
end

icls = []; %indices of values sorted by class

for j = 1:nmixgroups
    icls = [icls; iclscell{j}];
end

%sort by amount of white color - 5th row - in the overall mix
irow = [];
mxall = [];
ctr = 0;
for j = 1:nmixgroups
    if numel(iclscell{j}) > 0 %if the array is not empty
        mxarray = mixvalues(iclscell{j},:);
        mxall = [mxall; mxarray];

        [~,iarray] = sort(mxarray(:,5),'ascend'); %iarray is array of indices how strokes are actually organized
        irow = [irow; icls(ctr + iarray)]; %irow is a sorting order of clusters
        ctr = ctr + length(iarray);
    end
end

%then, obrain correct order of strokes by ind
ind = zeros(nStrokes,1);
ctr = 1;
for i = 1:nColors
    i_corr = irow(i); %ordinal number of i-th cluster
    for j = 1:nStrokes
        if ids(j) == i_corr %for each stroke of i_corr cluster, record its index
            ind(ctr) = j;
            ctr = ctr + 1;
        end
    end
end

strokes2 = zeros(nColors,nStrokes); %array of stroke indices
strCount = zeros(nColors,1); %number of strokes of each color

map2 = cell(nStrokes,1);

for i = 1:nStrokes
    i_cluster = ids(ind(i));  %index of cluster for i-th stroke
    i_sorted = ind(i);  %new index of i-th stroke

    strokes{i_sorted,1}.color = colorsFinal(i_cluster,:);
    %%% for testing %%%%%%%%%%%%
%     switch cls(i_cluster)
%         case 1
%             col2 = [1 0 0];
%         case 2
%             col2 = [1 1 0];
%         case 3
%             col2 = [0 1 0];
%         case 4
%             col2 = [0 0 1];
%     end
%     strokes{i_sorted,1}.color = 255*col2;
    %%%%

    strokes{i_sorted,1}.col8 = mixvalues(i_cluster,:);

    strCount(i_cluster) = strCount(i_cluster) + 1;
    strokes2(i_cluster,strCount(i_cluster)) = i_sorted;
    map2{i,1} = strokes{i_sorted,1}; %map with sorted colors
end

% Uncomment for debug
% figure(16); 
% imshow(map2imgColorCanvas(brushSize,canvas,map2,canvasColor));
% title('Semi-Final image');

%sort strokes by position in each cluster, to reduce the machine tool path
strokes3 = strokes2;
for i = 1:nColors % for each cluster
    waitbar(i/nColors,hw, ['Sorting strokes by position, i = ', num2str(i)]);
    curStrokeId = strokes2(i,1); %index in strokes
    strokes2(i,1) = 0;%0, because the i-th stroke is already taken
    j = 2;
    while j < strCount(i) %for the strokes in a cluster
        %найдем пару для curStroke
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
    iCol = irow(i); %correct order of clusters
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
flnm = [pathname,filename];
[filename, path] = uiputfile([flnm(1:end-4),'.plt'],'Save file as...');
filenameforsave = [path,filename];
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
function [err, candidate_out, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgraymix,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2, bsQuad, mixtype, vwhite)
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
                        %test for overlap: number of mixtype is >=, and
                        %amount of white is lower
                        if  ((canvasgraymix(Xl,Yl,1) == mixtype) && (canvasgraymix(Xl,Yl,2) < vwhite)) || (canvasgraymix(Xl,Yl,1) == 0) || (canvasgraymix(Xl,Yl,1) < mixtype)  %if canvas is free
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

%% DRAW PIECE
function [err, canvas_out, canvasgray_out, canvas2_out]= drawPiece(pX,pY,candidate,bs2,bsQuad,canvas,canvasColor, canvasgraymix,canvas2,meancol,col,col2,imggray,m,n, mixtype, vwhite)
canvas_out = canvas;
canvasgray_out = canvasgraymix;
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
err = abs(mean(canvas_out,3) - double(imggray));
end


