
%ЗАДАЕМ ИМЯ ФАЙЛА, КОТОРЫЙ ЧИТАЕМ. ОТКРЫВАЕМ ЭТОТ ФАЙЛ
warning off

pixScle = 4; %увеличение масштаба мм -> пиксель
antialiasing = 0; %change to 1 if needed

%table with colors
tablename = 'ModelTable600_.xls';
Ycell = cell(1,4); % 0 1 2 3 = MY1 MY2 CY CM
Wcell = cell(1,4);

for i = 1:4
    M1 = readmatrix(tablename,'Sheet',i);
    Ycell{i} = M1(:,4:6); %proportions
    Wcell{i} = M1(:,1:3); %colors in hsv
end

[filename,pathname] = uigetfile('*.plt','Выберите файл');
fid = fopen([pathname,filename],'r');

c = fscanf(fid,'%s');
cmds = split(c,';'); %РАЗБИВАЕМ С НА МАССИВ СMDS: CMDS[i] - i-я команда
color = zeros(1,1,3); %Начальный цвет
mixtype = 1;
ptX = 0; ptY = 0;

[brushWidth_mm, sideXmm, sideYmm, usespp, nColors] = GUI_Trace;

scalePLT = 1/40;%масштабирующий фактор PLT


n = sideXmm*40*pixScle*scalePLT;
m = sideYmm*40*pixScle*scalePLT;

canvasColor = 255;

canvas = uint8(canvasColor + zeros(m,n,3));
canvas2 = uint8(canvasColor + zeros(m,n,3));
%canvas3 = uint8(canvasColor + zeros(m,n,3)); %canvas for color numeration

canvasprop = zeros(m,n,5); %third dimension: [props,mixtype,colnumber]

scalefactor=scalePLT;% pixScle*

HeigthY = sideYmm*40*scalePLT;


brushSize = brushWidth_mm * pixScle;%толщина кисти
bs2 = ceil(brushSize/2);
bsQuad = bs2^2;
bsQuad2 = 2;

bsQuad_1 = (bs2 - 1)^2;
bsQuadp1 = (bs2 + 1)^2;

PWscl = 300/25.4;

hw = waitbar(0,'Выполнение команды PLT'); %статистика выполнения, вэйтбар
N = length(cmds);
strokectr = 0;
colctr = 0; %color counter
codecolor = [0 0 0];

strokelengths = [];
curlength = 0;

for j=1:N %по списку команд
    s = cmds{j,1}; %s - очередная команда
    if(~isempty(s)) %если она не пустая
        cmd = s(1:2); %берем первые 2 символа
        switch cmd %в зависимости от них делаем действия
            case 'PU'    %поднять кисть, привести в координаты - запомним ptX, ptY как новые координаты            
                i = strfind(s,',');
                sX = s(3: i - 1); 
                sY = s(i:end);
                %преобразум текст в число
                ptX = str2double(sY);
                ptY = str2double(sX);
                %преобразование масштаба
                ptX = (HeigthY*pixScle - round(ptX*scalePLT*pixScle));
                ptY = round(ptY*scalePLT*pixScle);
                
                %color2 = uint8([rand rand rand]*230); %случайный цвет для рисования карты мазков

                switch mixtype
                    case 1 % MY
                        col2 = [(0.5 + 0.5*rand) (0.5 + 0.5*rand) 0] + 0.5*[0 0 rand];
                    case 2 %YC
                        col2 = [0 (0.5 + 0.5*rand) 0] + 0.5*[rand 0 rand];
                    case 3 %CM
                        col2 = [(0.5 + 0.5*rand) 0 (0.5 + 0.5*rand)] + 0.5*[0 rand 0];
                end

                color2 = uint8(col2*230); %случайный цвет для рисования карты мазков
                
                strokectr = strokectr + 1;

                strokelengths = [strokelengths; curlength];
                curlength = 0;
                %j
            case 'PD'
                i = strfind(s,',');
                sX = s(3: i - 1);
                sY = s(i:end); 
                pX = str2double(sY);
                pY = str2double(sX);
                
        %        pX2 = (HeigthY - round(pX*scalePLT))*pixScle;
        %        pY2 = round(pY*scalePLT)*pixScle;
                 
                pX = (HeigthY*pixScle - round(pX*scalePLT*pixScle));
                pY = round(pY*scalePLT*pixScle);

                curlength = curlength + sqrt((ptX - pX)^2 + (ptY - pY)^2)/pixScle; %длина фрагмента
                
                %рисуем линию
                Np = max(abs(ptX - pX), abs(ptY - pY));    %наиб. число пикселов в линии по Х или Y
                for t = 0:1/Np:1 %параметрическое рисование линии, типа: x = at + b; y = ct + d;
                    xo = round(pX  + (ptX - pX)*t);
                    yo = round(pY + (ptY - pY)*t); %хо, уо - координаты точки по центру линии
                    for Xl = round(xo - bs2):round(xo + bs2)
                        for Yl = round(yo - bs2):round(yo + bs2)
                            %закрасим все пикселы вокруг очередной точки (хо, уо) попавшей в круг радиуса кисти (след от кисти)
                            if Xl > 0 && Xl <= m && Yl > 0 && Yl <= n 
                                r2 = (double(Xl) - double(xo))^2 + (double(Yl) - double(yo))^2; %радиус в квадрате
                                
                                if antialiasing
%                                     if r2 < bsQuad_1 %core of the line (bs2 - 1)^2
%                                         canvas(Xl,Yl,:) = color;
%                                         canvas3(Xl,Yl,:) = codecolor;
%                                     else
%                                         if r2 < bsQuadp1 % (bs2 + 1)^2
%                                             gamma = (bs2 - sqrt(r2))/2; % approximately equals (r - bs2)
%                                             if gamma > 1
%                                                 gamma = 1;
%                                             end
%                                             if gamma < 0
%                                                 gamma = 0;%%%%
%                                             end
%                                             canvas(Xl,Yl,:) = uint8(double(canvas(Xl,Yl,:))*(1 - gamma) + gamma*color); %ибо темные после светлых
%                                             %canvas3(Xl,Yl,:) = codecolor;
%                                             canvasprop(Xl,Yl,:) = [props,mixtype,colctr];
% 
%                                         end
%                                     end
                                else
                                    if r2 <= bsQuad %core of the line (bs2)^2
                                        canvas(Xl,Yl,:) = color;
                                        %canvas3(Xl,Yl,:) = codecolor;
                                        canvasprop(Xl,Yl,:) = [props,mixtype,colctr];
                                    end
                                end
                                %paint extra canvases with brushstrokes
                                if r2 < bsQuad2
                                    canvas2(Xl,Yl,:) = color2;
                                end
                               
                            end
                        end
                    end
                end
                ptX = pX; ptY = pY; %запомним текущую точку как предыдущую
                
            case 'PP' %смена цвета
                waitbar(j/N,hw,['Выполнение команды PLT ',num2str(j)]);
                i = strfind(s,',');
                spl = split(s(3:end),',');
                
                j
                s
                
                col8 = zeros(1,8);
                for i= 1:8
                    col8(i) = str2double(spl{i,1});
                end
%                 c_ = str2double(spl{1,1});
%                 m_ = str2double(spl{2,1});
%                 y_ = str2double(spl{3,1});
%                 b_ = str2double(spl{4,1});
%                 w_ = str2double(spl{5,1});
                
%                 [R,G,B] = PPtoPC([c_,m_,y_,b_,w_);
                [hsvcol,props,mixtype] = col82hsv(col8,Ycell,Wcell);
                [R,G,B] = hsv2rgb(hsvcol);
                color(1,1,1) = 255*R;
                color(1,1,2) = 255*G;
                color(1,1,3) = 255*B;

                colctr = colctr + 1

%                 bc = mod(colctr,100);
%                 gc = (mod(colctr,10000) - bc)/100;
%                 rc = (mod(colctr,1000000) - mod(colctr,10000))/10000;
%                 codecolor = [rc gc bc];
                
                figure(7); %выведем на экран холст мазков
                imshow(canvas2);
                figure(6); %выведем на экран холст с рисующимся изображением
                imshow(canvas);

%                 if ( j >= 120000)
%                pause
%                 end

           case 'PC' %смена цвета
                waitbar(j/N,hw,['Выполнение команды PLT ',num2str(j)]);
                i = strfind(s,',');
                spl = split(s(3:end),',');
                color(1,1,1)= str2double(spl{2,1});
                color(1,1,2) = str2double(spl{3,1});
                color(1,1,3) = str2double(spl{4,1});
                
                figure(6); %выведем на экран холст с рисующимся изображением
                imshow(canvas);
          
            case 'PW' %ширина кисти
                i = strfind(s,',');
                brushSize_txt = s(3:end);
                brushSize = str2double(brushSize_txt) * pixScle /PWscl;%толщина кисти
                bs2 = ceil(brushSize/2);
                bsQuad = bs2^2;
                bsQuad2 = 2;
        end
    end
end
close(hw); %закроем вэйтбар
fclose(fid); %закроем файл
figure(7);
imshow(canvas2); %выведем на экран холст мазков
figure(6);
imshow(canvas);%выведем на экран холст с готовым изображением

figure(14);
histogram(strokelengths,100);
title(['Stroke length histogram. Count = ',num2str(strokectr)]);
xlabel('length, mm');
ylabel('count');

imwrite(canvas,'Read.bmp');
save('Read.mat','canvasprop','colctr');