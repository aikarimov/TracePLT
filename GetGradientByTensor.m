% MODIFIED GRADIENT WITH A TENSOR APPROACH
function [U, V] = GetGradientByTensor(imggray, brushSize)
s = 6;
g = s*brushSize; % filter parameter

% p1 = 0.183;
% Dx = 0.5*[p1, 0, -p1;
%     1 - 2*p1, 0, 2*p1 - 1;
%     p1, 0, -p1];
% 
% Dy = transpose(Dx);
% fx = conv2(imggray, -Dx); % filter2(Dx,imggray);
% fy = conv2(imggray, -Dy);% filter2(Dy,imggray);

% by Farid and Simoncelli
% k  = [0.030320  0.249724  0.439911  0.249724  0.030320];
% d  = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
k  = [ 0.004711  0.069321  0.245410  0.361117  0.245410  0.069321  0.004711];
d  = [ 0.018708  0.125376  0.193091  0.000000 -0.193091 -0.125376 -0.018708];

fx = conv2(k, d, imggray, 'same');  % derivative horizontally (wrt X)
fy = conv2(d, k, imggray, 'same');  % derivative vertically (wrt Y)


%smoothing
% h =     [0         0         0         0         0         0;
%          0         0         0         0         0         0;
%          0         0    0.2500    0.2500         0         0;
%          0         0    0.2500    0.2500         0         0;
%          0         0         0         0         0         0;
%          0         0         0         0         0         0]; %fspecial('gaussian',6,0.1);

h = fspecial('gaussian',s*g,s);
%h = ones(g)/(g^2);

U = conv2(fx,h,'same');
V = conv2(fy,h,'same');

%tensor code
E = U.*U;
F = U.*V;
G = V.*V;

D = sqrt((E - G).^2 + 4*F.^2);
lam1 = 0.5*(E + G + D);

U = F;
V = lam1 - E;
 
%additional filtering U and V

%H = ones(g)/(g^2); % instead of fspecial('average',g);
%H = fspecial('gaussian',g,s);

%add to components for direction modulation
% epsx = 0.25*mean2(abs(U)); %tend to be horizontal
% epsy = 0; %tend to be vertical
% 
% U = U + epsx;
% V = V + epsy;
% 
% U  = conv2(U,H,'same'); 
% V  = conv2(V,H,'same');

% U  = imfilter(U,H,'replicate');
% V  = imfilter(V,H,'replicate');

scl = sqrt(U.^2 + V.^2); %scaling factor
scl = (scl == 0) + scl;
U = U./scl;
V = V./scl;

end