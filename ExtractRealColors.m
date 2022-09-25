%Extract real colors from read file and write into xls file

canvas = imread('Read.bmp'); % canvas contains scanned image
load('Read.mat','canvasprop','colctr'); % load brushstroke map

tablenameout = 'ModelTableExtra.xls';

[m,n,~] = size(canvasprop);

[mc,nc,~] = size(canvas);

%resize image if it does not fit to canvasprop
if (m ~= mc) || (n ~= nc)
    canvas = imresize(canvas,[m n]);
end

%create cell arrays: Y for props, W for hsv colors
Ycell = cell(1,4); % 0 1 2 3 = MY1 MY2 CY CM
Wcell = cell(1,4);

Ccell = cell(1,colctr);
proparray = zeros(colctr,4);

Hdiv = 0.4;

%read canvas data
for i = 1:m
    for j = 1:n
        rgbcol = reshape(canvas(i,j,:),[1 3]);

        colnumber = canvasprop(i,j,5);

        if colnumber > 0
            Ccell{colnumber} = [Ccell{colnumber}; rgbcol]; %add rgb color

            if isequal(proparray(colnumber,:),[0 0 0 0]) %if not recorded yet
                props = reshape(canvasprop(i,j,[1:3]),[1 3]);
                mixtype = reshape(canvasprop(i,j,4), [1 1]);

                hsvcol = rgb2hsv(double(rgbcol)/255);

                newmixtype = mixtype + 1;

                if mixtype == 1 && hsvcol(1) > Hdiv
                    newmixtype = 1;
                end

                proparray(colnumber,:) = [props, newmixtype];
            end
        end
    end
end

%averaging
for i = 1:colctr
    rgbcol = mean(Ccell{i},1);
    mixtype = proparray(i,4);
    props = proparray(i,1:3);
    Ycell{mixtype} = [Ycell{mixtype}; props]; % 0 1 2 3 = MY1 MY2 CY CM
    Wcell{mixtype} = [Wcell{mixtype};rgb2hsv(double(rgbcol)/255)];
end

%draw  figure
figure;
markers = {'o','s','d','^'};
markersizes = {30,20,30,50};
%visualize
for i = 1:4
    W = Wcell{i}; %HSV    
    scatter3(W(:,1),W(:,2),W(:,3),markersizes{i},hsv2rgb(W),markers{i},'filled','MarkerEdgeColor','k'); hold on
end

%save
%write model
for i = 1:4
    writematrix([Wcell{i}, Ycell{i}],tablenameout,'Sheet',i); %HSV first, proportions second
end

