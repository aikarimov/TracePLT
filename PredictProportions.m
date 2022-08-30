function [props, cls, hsvnew] = PredictProportions(varargin)
%
% [props, cls, hsvnew] = PredictProportions(hsvcol)

% [props, cls, hsvnew] = PredictProportions(hsvcol,Ycell, Wcell)
%
% hsvcol is HSV color, or N x 3 matrix of hsv colors
% Ycell, Wcell are cell arrays from ModelTable600.xls
%
% This function predicts proportions of the paint mix
% props = [a b c]
% a = color1/(color1 + color2), where color1 and color2 are primary colors
% b = black/(black + white)
% c = a / (a + b)
%
% cls is a type of mix: MY1, MY2, CY, CM
%
% if HSV is matrix, props and cls are N x 3 and N x 1, respectively
tablename = 'ModelTable600.xls';

K = 20;

tol = 0.2; %tolerance of color error in hsv

sat = @(x) min(max(x, 0), 1); %saturation function

Ycell = cell(1,4); % 0 1 2 3 = MY1 MY2 CY CM
Wcell = cell(1,4);
Tbl = [];
Clss = [];

hsvcol = varargin{1,1};

if nargin == 3
    Ycell = varargin{1,2};
    Wcell = varargin{1,3};
end

for i = 1:4
    if nargin == 1
        M1 = readmatrix(tablename,'Sheet',i);
        Ycell{i} = M1(:,1:3);
        Wcell{i} = M1(:,4:6);
    end

    [Nset,~] = size(Ycell{i});
    Tbl = [Tbl; Ycell{i}];
    Clss = [Clss;i*ones(Nset,1)]; %'Magenta1','Magenta2','Yellow','Cyan'
end

%now, create classification table and class list

[N,~] = size(Tbl);

[Nhsv,~] = size(hsvcol);
hsvnew = hsvcol;
props = zeros(Nhsv,3);
cls = zeros(Nhsv,1);
hsvarray = zeros(Nhsv,3);

for i = 1:Nhsv
%     %first, check for consistency - if the color is outside a possible range
%     tol = 0.1; %tolerance, 1/2 of width of a color stripe
% 
%     H = hsvcol(i,1); %take a slice in H
%     slicepts = [];
%     Npts = 0;
%     for k = 1:N %among all points
%         if abs(Tbl(k,1) - H) <= tol
%             slicepts = [slicepts; Tbl(k,:)];
%             Npts = Npts + 1;
%         end
%     end
%     %inside a slice, divide space into quadrants
%     quads = zeros(1,4);
%     for k = 1:Npts %among points in a slice
%         sk = slicepts(k,2);
%         vk = slicepts(k,3);
%         s = hsvcol(i,2);
%         v = hsvcol(i,3);
% 
%         if  sk > s
%             if vk > v
%                 quads(1) = 1;
%             else
%                 quads(4) = 1;
%             end
%         else
%             if vk > v
%                 quads(2) = 1;
%             else
%                 quads(3) = 1;
%             end
%         end
%     end
%     if sum(quads) < 4 %this means that interpolation is impossible
%         %take Kp closest points and select mean values
%         Kp = 10;
%         dst = sqrt(3)*ones(Npts,1); %distances - taken as maximum range
%         for k = 1:Npts %among points in a slice
%             dst(k) = sqrt( (hsvcol(i,:) - slicepts(k,:))*(hsvcol(i,:) - slicepts(k,:))'); %Euclidean distance in (H,S,V) space
%         end
%         [~,I] = sort(dst,1,'ascend');
%         s = mean(slicepts(I(1:Kp),2));
%         v = mean(slicepts(I(1:Kp),3));
%         hsvnew(i,2:3) = [s, v]; %modify color
%     end

%alternative: no check!
    
    hsvnew(i,:) = hsvcol(i,:);

    %then, replace the current color with the accessible one
    hsvcolcur = hsvnew(i,:);
    
    %make prediction from the closest point
    dst = sqrt(3)*ones(N,1); %distances - taken as maximum range
    for k = 1:N %among all points
        dst(k) = (hsvcolcur - Tbl(k,:))*(hsvcolcur - Tbl(k,:))'; %squared Euclidean distance in (H,S,V) space
        %dst(k) = abs(hsvcolcur(1) - Tbl(k,1)); %Euclidean distance in H only
    end
    [~,I] = sort(dst,1,'ascend');
    clsvect = Clss(I(1));
    cls(i) = clsvect;

    %then, find second possible class
    flag = 1;
    ctr = 2;
    class2 = cls(i);
    while flag
        if Clss(I(ctr)) ~= cls(i) %if 
            flag = 0;
            class2 = Clss(I(ctr)); %this class will be tested if cls(i) is wrong
        else
            ctr = ctr + 1;
        end
    end

    flag = 1; ctr = 1; 
    
    props0   = zeros(1,3); 
    propscur = zeros(1,3);
    
    err0 = 0;

    while flag
    
        %take the corresponding sets
        Y = Ycell{cls(i)};  %Y for HSV
        W = Wcell{cls(i)}; %W for proportions

        [NY,~ ] = size(Y);
        dst = sqrt(3)*ones(NY,1); %distances - taken as maximum range
        for k = 1:NY %take all points
        %    dst(k) = sqrt( (hsvcolcur - Y(k,:))*(hsvcolcur - Y(k,:))'); %Euclidean distance
            dst(k) = (hsvcolcur - Y(k,:))*(hsvcolcur - Y(k,:))'; %squared Euclidean distance
        end
        [~,I] = sort(dst,1,'ascend');
        
        %assign weights
        d = 1./dst(I(1:K));
        D = diag(d);

        %take first K points
%         dmax = 2;
%             [H, T, ~] = PolyRegression(Y(I(1:K),:),W(I(1:K),:),0,dmax,0.01,0.01);
%             propscur = PolyPredict(hsvcolcur, H, T);%predicted proportions
        %
        for j = 1:3
            X = Y(I(1:K),:);
            E = zeros(K,4); %evaluated polynomial
            T = [zeros(1,3); eye(3)]; %monomial orders
            h = eye(4);
            for k = 1:4
                E = E + h(k,:).*prod(X.^repmat(T(k,:),K,1),2);
            end
            
            h = (E'*D*E)\(E'*D*W(I(1:K),j)); %WLS

            %predict proportion
            p = 0;
            X = hsvcolcur;
            for k = 1:4
                p = p + h(k,:).*prod(X.^repmat(T(k,:),1,1),2);
            end

            propscur(j) = p;
        end

        %get into ranges [0,1]
        propscur = sat(propscur);

        %comment out, for test!
        %flag = 0; %go out of the loop

        %test inversion
        hsvinv = prop2hsv(propscur, cls, Wcell,Ycell);
        
        h_alt = [hsvinv(1),hsvinv(1) - 1 ]; %alternative versions of h
        [~,I] = min(abs(h_alt - hsvcolcur(1)));
        hsvinv2 = hsvinv;
        hsvinv2(1) = h_alt(I);
        err = sqrt( (hsvcolcur - hsvinv2)*(hsvcolcur - hsvinv2)');

        if ctr < 2 
            %if the values are calculated for the first time
            if err > tol
                ctr = ctr + 1;
                cls(i) = class2; %try to apply an alternative class
                props0 = propscur;
                err0 = err;
            else
                flag = 0; %go out of the loop
            end
        else
            %if the values are calculated for the second time
            if err > err0 %if error is greater, return to first variant
                propscur = props0;
                cls(i) = clsvect;

                hsvinv = prop2hsv(propscur, cls, Wcell,Ycell);
            end
            %anyway, go out of the loop
            flag = 0; %go out of the loop
        end
    end
    props(i,:) = propscur;
    hsvarray(i,:) = hsvinv;%hsvnew(i,:);
end
hsvnew = hsvarray;
end
