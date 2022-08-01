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

K = 20;

sat = @(x) min(max(x, 0), 1); %saturation function
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

%take the corresponding sets
if cls == 1
    Y = [Ycell{1};Ycell{2}];  %Y for proportions
    W = Wcell{1}; %W for HSV
    W(:,1) = W(:,1) - 1; %for MY with H > 0.4 make H negative
    W = [W; Wcell{2}];
else
    Y = Ycell{cls + 1};  %Y for proportions
    W = Wcell{cls + 1}; %W for HSV
end
[NY,~ ] = size(Y);
dst = sqrt(3)*ones(NY,1); %distances - taken as maximum range
for k = 1:NY %take all points
    %dst(k) = sqrt( (props - Y(k,:))*(props - Y(k,:))'); %Euclidean distance 
    dst(k) = (props - Y(k,:))*(props - Y(k,:))'; %squared Euclidean distance 
end
[~,I] = sort(dst,1,'ascend');

%assign weights
d = 1./dst(I(1:K));
D = diag(d);

% [H, T, ~] = PolyRegression(Y(I(1:K),:),W(I(1:K),:),0,dmax,eta,eps);
% hsvcol = PolyPredict(props, H, T);%predicted proportions

hsvcol = zeros(1,3);
%take first K points
for j = 1:3
    X = Y(I(1:K),:);
    E = zeros(K,4); %evaluated polynomial
    T = [zeros(1,3); eye(3)]; %monomial orders
    h = eye(4);
    for k = 1:4
        E = E + h(k,:).*prod(X.^repmat(T(k,:),K,1),2);
    end
    %h = (E'*E)\(E'*W(I(1:K),j)); %OLS
    h = (E'*D*E)\(E'*D*W(I(1:K),j)); %WLS

    %predict proportion
    p = 0;
    X = props;
    for k = 1:4
        p = p + h(k,:).*prod(X.^repmat(T(k,:),1,1),2);
    end
    hsvcol(j) = p;
end

if hsvcol(1) < 0
    hsvcol(1) = hsvcol(1) + 1; %use this because models predict with shift
end
%get into ranges [0,1]
hsvcol = sat(hsvcol); 
mixtype = cls;

end
