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

%K = 20;
K = 25;

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

%matrix wor weighting errors
M = diag([1, 1/16, 1]);

for i = 1:Nhsv   
    hsvnew(i,:) = hsvcol(i,:);
    hsvcolcur = hsvnew(i,:);

    %make prediction of class from the Ks closest point
    Ks = 10;
    dst = sqrt(3)*ones(N,1); %distances - taken as maximum range
    for k = 1:N %among all points
        dst(k) = (hsvcolcur - Tbl(k,:))*(hsvcolcur - Tbl(k,:))'; %squared Euclidean distance in (H,S,V) space
        %dst(k) = abs(hsvcolcur(1) - Tbl(k,1)); %Euclidean distance in H only
    end
    [~,I] = sort(dst,1,'ascend'); %find
    clvec = Clss(I(1:Ks));
    %clsvect = mode(clvec);
    %my faster code for mode
    classes = zeros(1,4);
    for j = 1:Ks
        classes(clvec(j)) = classes(clvec(j)) + 1;
    end
    [~,clsvect] = max(classes);




    cls(i) = clsvect;

    %then, find second possible class: mode of other K closest
    %points, excluding class y
    ctr = 1; ctri = 1;
    class2vect = zeros(K,1);
    while (ctr <= length(Clss)) && (ctri <= K)
        if Clss(I(ctr)) ~= cls(i) %if
            class2vect(ctri) = Clss(I(ctr)); %this class will be tested if cls(i) is wrong
            ctri = ctri + 1;
        end
        ctr = ctr + 1;
    end
    %class2 = mode(class2vect(1:(ctri-1)));
    %my faster code for mode
    classes = zeros(1,4);
    for j = 1:(ctri - 1)
        classes(class2vect(j)) = classes(class2vect(j)) + 1;
    end
    [~,class2] = max(classes);


    flag = 1; 
    
    ctr = 1; 
    
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
            dst(k) = (hsvcolcur - Y(k,:))*M*(hsvcolcur - Y(k,:))'; %squared Euclidean distance
        end
        [~,I] = sort(dst,1,'ascend');
        
        %assign weights - optional
        d = 1./dst(I(1:K));
        D = diag(d);

        %for non-weighted regression
        %D = eye(K);

        %for determining possible color - weighted linear interpolation
        ds = sqrt(d);
        hsvpossible = ds'*Y(I(1:K),:)/sum(ds); %linear interpolation
        al = 0.7;
        hsvnew(i,:) = (1 - al)*hsvpossible + al*hsvcolcur; %set 0.3 parameter for real color, for expolring the space
        hsvcolcur = hsvnew(i,:);

        %

        %take first K points, make linear regression
        for j = 1:3
            X = Y(I(1:K),:);
            E = zeros(K,4); %evaluated polynomial
            T = [zeros(1,3); eye(3)]; %monomial orders
            h = eye(4);
            for k = 1:4
                E = E + h(k,:).*prod(X.^repmat(T(k,:),K,1),2);
            end
            
            delt = 1e-14; %for Tikhonov regularization

            h = (E'*D*E + delt*eye(4))\(E'*D*W(I(1:K),j)); %WLS

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         flag = 0;
%         hsvinv = hsvcolcur;
%         continue; %go out of the loop

        %test inversion
        hsvinv = prop2hsv(propscur, cls, Wcell,Ycell);
        
        h_alt = [hsvinv(1),hsvinv(1) - 1 ]; %alternative versions of h
        [~,I] = min(abs(h_alt - hsvcolcur(1)));
        hsvinv2 = hsvinv;
        hsvinv2(1) = h_alt(I);

        hsvdemanded = hsvcol(i,:);
        err = sqrt( (hsvdemanded - hsvinv2)*(hsvdemanded - hsvinv2)');

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
                %hsvinv = prop2hsv(propscur, cls, Wcell,Ycell);
            end
            %anyway, go out of the loop
            flag = 0; %go out of the loop
        end
    end
    props(i,:) = propscur;
    %hsvarray(i,:) = hsvinv;
    hsvarray(i,:) = hsvnew(i,:);
end
hsvnew = hsvarray;
end
