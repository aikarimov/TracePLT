function [hsvcol] = prop2hsv(varargin)
% hsvcol = prop2hsv(props, cls)
%
% hsvcol = prop2hsv(props, cls, Ycell,Wcell)
%
% hsvcol is HSV color, or 1 x 3 vector of hsv colors
% Ycell is a cell array of proportions, Wcell is a cell of HSV colors from
% ModelTable600.xls
%
% This function predicts HSV color from proportions and class
% cls = {1 2 3 4} <-> {MY1 MY2 CY CM}

tablename = 'ModelTable600.xls';

%K = 150;
K = 250;

sat = @(x) min(max(x, 0), 1); %saturation function
props = varargin{1,1};
cls =  varargin{1,2};

if nargin == 2
    Ycell = cell(1,4); % 0 1 2 3 = MY1 MY2 CY CM
    Wcell = cell(1,4);

    for i = 1:4
        M1 = readmatrix(tablename,'Sheet',i);
        Ycell{i} = M1(:,4:6); %proportions
        Wcell{i} = M1(:,1:3); %colors in hsv
    end
end

if nargin == 4
    Ycell = varargin{1,3};
    Wcell = varargin{1,4};
end

Npts = size(props,1);
hsvcol = zeros(Npts,3);

for i = 1:Npts

    %take the corresponding sets
    if cls(i) <= 2
        Y = [Ycell{1};Ycell{2}];  %Y for proportions
        W = Wcell{1}; %W for HSV
        W(:,1) = W(:,1) - 1; %for MY with H > 0.4 make H negative
        W = [W; Wcell{2}];
    else
        Y = Ycell{cls(i)};  %Y for proportions
        W = Wcell{cls(i)}; %W for HSV
    end
    [NY,~ ] = size(Y);
    dst = sqrt(3)*ones(NY,1); %distances - taken as maximum range

    for k = 1:NY %take all points
        %dst(k) = sqrt( (props(i,:) - Y(k,:))*(props(i,:) - Y(k,:))'); %Euclidean distance
        dst(k) = (props(i,:) - Y(k,:))*(props(i,:) - Y(k,:))'; %squared Euclidean distance
    end
    [~,I] = sort(dst,1,'ascend');

    %assign weights
    d = 1./(dst(I(1:K)) + 1e-4);
    D = diag(d);

    hsvcolcur = zeros(1,3);

    N = 4;
    T = [0     0     0;
         1     0     0;
         0     1     0;
         0     0     1]; %deglexord(0,1,3);
    
    %take first K points
    K = min(NY,K); %decrease K if needed
    X = Y(I(1:K),:);
    
    E = zeros(K,N); %evaluated polynomial
    h = eye(N);
    for k = 1:N
        pwrs = repmat(T(k,:),K,1);
        E = E + h(k,:).*prod(X.^pwrs,2);
    end

    for j = 1:3
        V = W(I(1:K),j);
        %h = (E'*E)\(E'*V); %OLS
        h = (E'*D*E + 1e-16)\(E'*D*V); %WLS
        %predict proportion
        p = 0;
        x = props(i,:);
        for k = 1:N
            p = p + h(k,:).*prod(x.^T(k,:),2);
        end
        hsvcolcur(j) = p;
    end

    if hsvcolcur(1) < 0
        hsvcolcur(1) = hsvcolcur(1) + 1; %use this because models predict with shift
    end
    %get into ranges [0,1]
    hsvcol(i,:) = sat(hsvcolcur);
end
end
