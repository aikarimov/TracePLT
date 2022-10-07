%Трассировщик PLT 
%% TRACE - MAIN
close all
warning off

% ОСНОВНЫЕ НАСТРОЙКИ:
canvasColor = 255; %тон холста (белый)
itersMinOverlap = 1; %1 итераций с малым перекрытием
TotalIter = 3; %3 общее число итераций

minlenfactor = [3 1 0]; % мин. длина мазка в диаметрах кисти
maxlenfactor = 30; % макс. длина мазка в диаметрах кисти

minOverlap = 0.6;%0.6; % макс. коэффициент перекрытия начальный
maxOverlap = 0.8; % 0.8;
pixTol = 9; % 6; %возможное отклонение цвета от исходника на конце, в RGB
pixTol2 = 100; %возможное отклонение цвета в среднем
do_blur = 0; %нужно ли размывать холст
canvasEps = 2; % +\- погрешность задания цвета холста

gonormal = 1; % класть мазки: 0 - по градиенту, 1 - перпендикулярно градиенту

% ЧТЕНИЕ ФАЙЛА РАСТРОВОГО ИЗОБРАЖЕНИЯ

[filename,pathname] = uigetfile({'*.png;*.jpg;*.bmp','Images (*.png;*.jpg;*.bmp)';'*.*','All files (*.*)'},'Select the file');
img = imread([pathname,filename]); %img contains image

[m, n] = size(img(:,:,1));

[brushWidth_mm, canvasW_mm, canvasH_mm, ~, nColors] = GUI_Trace;

% ЧИТАЕМ ТАБЛИЦУ ДАННЫХ С ПРОПОРЦИЯМИ И ЦВЕТАМИ
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

brushSize = round(brushWidth_mm/sfX) %толщина кисти
bs2 = ceil(brushSize/2);

MaxIters = brushSize;%макс. число попыток найти новый отрезок
maxLen = brushSize*maxlenfactor; %макс. длина мазка

%размытие в соответствии с радиусом кисти - опционально
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

nStrokes = 0;
accepted = 0;

bsQuad = (double(brushSize)/2)^2;
hw = waitbar(0,'wait...');
canvas2 = zeros(m,n,3); %холст с искусственными цветами мазков

for kk = 1:TotalIter

    minLen = brushSize*minlenfactor(kk); %мин. длина мазка - параметр применяется только на 1 итерации

    if kk == 2 %на второй итерации заменяем изображение на неразмытое
        img = img_unblurred;
    end
    
    if (kk > itersMinOverlap) %на итерации номер больше itersMinOverlap изменяем коэффициент перекрытия мазков
        ovf = maxOverlap;
    else
        ovf = minOverlap;
    end
    j = 1;
    %randomize stroke seeds
%     ivect = randperm(m);
%     jvect = randperm(n);

    pvect = randperm(m*n);

    for ictr = 1:m            
        for jctr = 1:n
            
            pctr = ((ictr - 1)*n + jctr);
            pij = pvect(pctr) - 1;
            j = mod(pij,n) + 1;
            i = (pij - j + 1)/n + 1;

%             i = ictr;
%             j = jctr;

            if (err(i,j) > pixTol) || canvas2(i,j) == 0
                %если в пикселе высокая погрешность
                pX = i; pY = j; %prev X,Y
                [col, meancol] = getMeanColor(img, pX, pY, bs2, bsQuad, m, n); %найдем средий цвет по круглой области радиусом с радиус кисти

                [props, mixtype, hsvnew]  = PredictProportions(rgb2hsv(double(meancol)/255), Ycell, Wcell);%определим пропорции и тип смеси в этой области
                col8paints = prop2pp(props,mixtype); % пропорции красок в реальном выражении
                vwhite = col8paints(5); % объем белого цвета в смеси - для сортировки мазков по перекрытию

                %нарисуем синтетическую карту мазков с цветом,
                %соответствующим типу смеси

%                 switch mixtype
%                     case 1
%                         col2 = [(0.5 + 0.5*rand) 0 0] + 0.5*[0 rand rand];
%                     case 2
%                         col2 = [(0.5 + 0.5*rand) (0.5 + 0.5*rand) 0] + 0.5*[0 0 rand];
%                     case 3
%                         col2 = [0.01 (0.5 + 0.5*rand) 0] + 0.5*[rand 0 rand];
%                     case 4
%                         col2 = [0.01 0 (0.5 + 0.5*rand)] + 0.5*[0 0 rand];
%                 end
                col2 = hsv2rgb(hsvnew);
                
                %структура мазка
                strokeStructure = struct('color',col,'Xs',[pX],'Ys',[pY],'length',0,'col8',col8paints,'props',props,'cls',mixtype);
                nPts = 1; %number of points in the stroke
                endedStroke = 0; %флаг для определения, завершен ли мазок

                %copy canvases - to backup if the stroke is too short
                canvas_copy =  canvas;
                canvasgraymix_copy = canvasgraymix;
                canvas2_copy = canvas2;

                while endedStroke == 0 %until the stroke is not ended
                    %find new direction
                    candidate = struct('X',-1,'Y',-1,'err',realmax);
                    ctr = 1;
                    while  ctr < MaxIters
                        [cosA, sinA] = getDirection(pX,pY, strokeStructure,U,V,gonormal);

                        %candidate r
                        r = MaxIters - ctr;

                        candidate.X = pX + round(r*cosA); %new X
                        candidate.Y = pY + round(r*sinA); %new Y
                        
                        %test new fragment of the stroke
                        [~, candidate, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgraymix,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2,bsQuad, mixtype, vwhite);
                        
                        if accepted %if error is small, accept the stroke immediately
                            ctr = MaxIters;
                        else
                            %also try opposite direction
%                             candidate.X = pX - round(r*cosA); %new X
%                             candidate.Y = pY - round(r*sinA); %new Y
%                             
%                             %test new fragment of the stroke
%                             [~, candidate, accepted] = testNewPiece(pX,pY,img,canvasColor,canvasEps,canvas2,canvasgraymix,pixTol,pixTol2, col, meancol,m,n,ovf,candidate, bs2,bsQuad, mixtype, vwhite);
%                             if accepted %if error is small, accept the stroke immediately
%                                 ctr = MaxIters;
%                             end
                        end
                        
                        ctr = ctr + 1;
                    end
                    
                    %draw the stroke fragment
                    if (accepted)
                        [err, canvas, canvasgraymix, canvas2] = drawPiece(pX,pY,candidate,bs2,bsQuad,canvas, canvasgraymix,canvas2,col,col2,imggray,m,n, mixtype, vwhite);
                               
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
                    if strokeStructure.length < minLen % && kk <= itersMinOverlap
                        %if the stroke is too short, and the iteration number
                        %is for long non-overlapping strokes
                        %backup canvases
                        canvas =  canvas_copy;
                        canvasgraymix = canvasgraymix_copy;
                        canvas2 = canvas2_copy;
                    else
                        %if the stroke is appropriate
                        nStrokes = nStrokes + 1;
                        dcol = strokeStructure.color;
                        strokeStructure.color = [dcol(1,1,1) dcol(1,1,2) dcol(1,1,3)];
                        strokes{nStrokes,1} = strokeStructure;

                        %new message, update message bar
                        waitbar(((ictr - 1)*n+jctr)/m/n,hw, ['iter = ', num2str(kk), 'i,j = ', num2str(i),',',num2str(j),' nStrokes = ', num2str(nStrokes)]);
                    end
                end
            end
        end
        %after each row iteration, show canvas
        figure(2);
        imshow(uint8(canvas));title('Canvas');

        figure(3);
        imshow(canvas2);title('Brushstroke map');
    end

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
colCell = cell(1,nmixgroups);

idcell = cell(1,nmixgroups); 
for i = 1:nStrokes
    mixtype = strokes{i,1}.cls;
    mixCell{mixtype} = [mixCell{mixtype}; strokes{i,1}.props]; % proportions array
    colCell{mixtype} = [colCell{mixtype}; strokes{i,1}.color]; % color array

    idcell{mixtype} =  [ idcell{mixtype}; i]; %real order of stroke in cell array
end

%now, rescale K for each mix group and perform internal k-means among
%proportions
nColors2 = 0;

ids = [];
cls = [];
mixvalues = []; %arrange col8
colorsFinal = []; %arrange RGB colors
propvalues = []; %arrange proportions

i2i = [];

ctr = 0;
% кластеризуем цвета по каждому типу смеси
for i = 1:nmixgroups
    idx = idcell{i};
    i2i = [i2i; idx];
    mixArray = mixCell{i};
    colArray = colCell{i}; %rgb
    Ncl = size(mixArray,1);

    if Ncl > 0
        Ki = ceil(Ncl / nStrokes * nColors); % число кластеров для данного типа смеси, пропорциональное числу мазков
        if Ki >= Ncl
            Ki = Ncl;
            idsk = (1:Ncl)'; %indices directly transform into clusters
        else
            %idsk = kmeans(mixArray,Ki); %clustering the proportions, reducing the number of colors
            idsk = kmeans(colArray,Ki); %clustering the colors in RGB, reducing the number of colors
        end

        nColors2 = nColors2 + Ki;

        %propsFinal = zeros(Ki,3);
        colFinal = zeros(Ki,3);
        nels = zeros(1,Ki);
        for j = 1:Ncl
            %dprops = mixArray(j,:);
            %propsFinal(idsk(j),:) = propsFinal(idsk(j),:) + dprops; %sum of proportions before averaging
            dcols = colArray(j,:);
            colFinal(idsk(j),:) = colFinal(idsk(j),:) + dcols; %sum of colors before averaging

            nels(idsk(j)) = nels(idsk(j)) + 1;
        end

        %averaging proportions in each cluster
%         for j = 1:Ki
%             propsFinal(j,:) = double(propsFinal(j,:))./nels(j);
%         end

        %averaging colors in each cluster
        for j = 1:Ki
            colFinal(j,:) = colFinal(j,:)./nels(j)/255;
        end

        %transform proportions into colors
        %hsvcolors = prop2hsv(propsFinal, i*ones(Ncl,1), Wcell,Ycell);
        hsvcolors = rgb2hsv(colFinal);
        [propsFinal,mixtyps] = PredictProportions(hsvcolors,Ycell,Wcell);
        colorsFinali = uint8(255*hsv2rgb(prop2hsv(propsFinal,mixtyps,Wcell,Ycell)));

        %rgbnew = uint8(255*colFinal);
        %colorsFinali = rgbnew; %rewrite final colors

        mixvaluesi = prop2pp(propsFinal,i*ones(Ncl,1));

        ids = [ids; idsk + ctr]; % global indices
        %cls = [cls; i*ones(Ki,1)]; % classes of colors %теперь надо учесть, что некоторые типы смесей поменялись!
        cls = [cls;mixtyps];
        mixvalues = [mixvalues; mixvaluesi]; % 8 цветов
        colorsFinal = [colorsFinal; colorsFinali]; % финальные цвета
        propvalues = [propvalues; propsFinal]; % пропорции смесей

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
    strokes{i_sorted,1}.col8 = mixvalues(i_cluster,:);
    strokes{i_sorted,1}.props = propvalues(i_cluster,:);

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

% сохраняем
flnm = [pathname,filename];
%[filename, path] = uiputfile([flnm(1:end-4),'.plt'],'Save file as...');
path = pathname; %automatic save
filenameforsave = [flnm(1:end-4),'.plt'];

%filenameforsave = [path,filename];
SavePLT_8paints(filenameforsave, map,n,m,canvasW_mm,canvasH_mm, brushWidth_mm);











