function [hsvcol, props, mixtype] = col82hsv(varargin)
% [hsvcol, props, mixtype] = col82hsv(col8)
%
% [hsvcol, props, mixtype] = col82hsv(col8, Ycell, Wcell)
%
% hsvcol is HSV color, or 1 x 3 vector of hsv colors
% props are proportions of paints
% mixtype = {1 = MY, 2 = YC, 3 = CM} type of primary colors mix
% Ycell, Wcell are cell arrays from ModelTable600.xls
%
% This function predicts HSV color from absolute values of
% 8 color proportions:
% col8 = [C M Y B W 0 0 0]

tablename = 'ModelTable600.xls';

col8 = varargin{1,1};

if nargin == 1
    Ycell = cell(1,4); % 0 1 2 3 = MY1 MY2 CY CM
    Wcell = cell(1,4);

    for i = 1:4
        M1 = readmatrix(tablename,'Sheet',i);
        Ycell{i} = M1(:,4:6); %proportions
        Wcell{i} = M1(:,1:3); %colors in hsv
    end
end

if nargin == 3
    Ycell = varargin{1,2};
    Wcell = varargin{1,3};
end

[props, cls] = col82props(col8);
mixtype = cls;

if cls == 1
    if props(1) < 0.4
         cls = 2;
    end
else
    cls = cls + 1;
end

hsvcol = prop2hsv(props, cls, Ycell,Wcell);

end
