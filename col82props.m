function [props, cls] = col82props(col8)
% [props, cls] = col82hsv(col8)
%
% props are proportions of paints
% cls = {1 = MY, 2 = YC, 3 = CM} type of primary colors mix
%
% Convert col8 into proportions 
% col8 = [C M Y B W 0 0 0]

vTotal = sum(col8);
hue = sum(col8(1:3));

%determine class of the color
cls = 1; %default class
a = 0;
if hue > 0
    if (col8(2) > 0 && col8(3) > 0) || (col8(2) > 0 && col8(1) == 0 && col8(3) == 0)
        cls = 1; %MY
        a = col8(2)/(col8(2) + col8(3));
    else
        if (col8(3) > 0 && col8(1) > 0) || (col8(3) > 0 && col8(1) ==0 && col8(2) == 0)
            cls = 2; %YC
            a = col8(3)/(col8(1) + col8(3));
        else
            cls = 3; %CM
            a = col8(1)/(col8(1) + col8(2));
        end
    end
end

if col8(4) +  col8(5) > 0
    b = col8(4)/sum(col8(4) + col8(5));
else
     b = 0;
end

c = hue/vTotal;

props = [a b c];