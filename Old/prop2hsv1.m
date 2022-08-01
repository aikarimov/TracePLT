function [hsvcol] = prop2hsv(varargin)
%hsvcol = prop2hsv(props, cls)
%
%hsvcol = prop2hsv(props, cls, Ycell,Wcell)
%
% hsvcol is HSV color, or 1 x 3 vector of hsv colors
% Wcell is a cell array of proportions, Ycell is a cell of HSV colors from
% ModelTable600.xls
%
% This function predicts HSV color from proportions and class
% cls = {1 2 3 4} <-> {MY1 MY2 CY CM}

tablename = 'ModelTable600.xls';

K = 150;

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

%take the corresponding sets
if cls <= 2
    Y = [Ycell{1};Ycell{2}];  %Y for proportions
    W = Wcell{1}; %W for HSV
    W(:,1) = W(:,1) - 1; %for MY with H > 0.4 make H negative
    W = [W; Wcell{2}];
else
    Y = Ycell{cls};  %Y for proportions
    W = Wcell{cls}; %W for HSV
end
[NY,~ ] = size(Y);
dst = sqrt(3)*ones(NY,1); %distances - taken as maximum range
for k = 1:NY %take all points
    dst(k) = sqrt( (props - Y(k,:))*(props - Y(k,:))'); %Euclidean distance 
end
[~,I] = sort(dst,1,'ascend');

hsvcol = zeros(1,3);
% N = 10;
% T = [0     0     0;
%      1     0     0;
%      0     1     0;
%      0     0     1;
%      2     0     0;
%      1     1     0;
%      1     0     1;
%      0     2     0;
%      0     1     1;
%      0     0     2]; %deglexord(0,2,3);

N = 4;
T = [0     0     0;
     1     0     0;
     0     1     0;
     0     0     1]; %deglexord(0,1,3);

%take first K points
for j = 1:3
    X = Y(I(1:K),:);
    E = zeros(K,N); %evaluated polynomial
    h = eye(N);
    for k = 1:N
        E = E + h(k,:).*prod(X.^repmat(T(k,:),K,1),2);
    end
    h = (E'*E)\(E'*W(I(1:K),j)); %OLS

    %predict proportion
    p = 0;
    X = props;
    for k = 1:N
        p = p + h(k,:).*prod(X.^repmat(T(k,:),1,1),2);
    end
    hsvcol(j) = p;
end

if hsvcol(1) < 0
    hsvcol(1) = hsvcol(1) + 1; %use this because models predict with shift
end
%get into ranges [0,1]
hsvcol = sat(hsvcol); 

end
