function col8paints = prop2pp(props,cls)
%cls = MY1,MY2,YC,CM
% props = [a b c]
% a = color1/(color1 + color2), where color1 and color2 are primary colors
% b = black/(black + white)
% c = hue / all

vTotal = 3000;

[N,~] = size(props);
col8paints = zeros(N,8);
for i = 1:N
    C = 0; M = 0; Y = 0; %initiate colors
    a = props(i,1);
    b = props(i,2);
    c = props(i,3);

    C1 = floor(a*c*vTotal);
    C2 = floor((1 - a)*c*vTotal);
    B  = floor(b*(1-c)*vTotal);
    W  = floor((1 - b)*(1 - c)*vTotal);

    switch cls(i)
        case {1,2}
            M = C1; Y = C2;
        case 3
            Y = C1; C = C2;
        case 4
            C = C1; M = C2;
    end

    col8paints(i,:) = [C, M, Y, B, W, 0, 0, 0]; %Och Or G are not included yet
end
end
