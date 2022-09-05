%process data from files
%backward approximation
%color to proportions
% using closest K points
close all

Hdiv =  0.4; %H value divising magenta-yellow pair into 2 sets

tablename = 'ModelTable600_very_old.xls';
Ycell = cell(1,4); % 0 1 2 3 = MY1 MY2 CY CM
Wcell = cell(1,4);

for i = 1:4
    M1 = readmatrix(tablename,'Sheet',i);
    Ycell{i} = M1(:,1:3);
    Wcell{i} = M1(:,4:6);
end

textcell = {'MY1', 'MY2', 'CY', 'CM'};

%now, in each group perform a regression using polynomial fits

%save model
%save('ClassificationModel',Mdl);

%try to distinguish visually
figure(100);

markers = {'o','o','d','^'};
markersizes = {20,20,30,50};

for i = [2,3,4,1]
    Y = Ycell{i}; %HSV
    scatter3(Y(:,1),Y(:,2),Y(:,3),markersizes{i},hsv2rgb(Y),markers{i},'filled','MarkerEdgeColor','k'); hold on
    axis square
end

legend('MY','YC','CM');

xlabel('H');
ylabel('S');
zlabel('V');


%visualize in prop system

spn = [1 1 2 3]; %subplot number
figure(200);

titlecell = {'MY', 'CY', 'CM'};

for i = 1:4
    subplot(1,3,spn(i));
    Y = Ycell{i}; %HSV 
    W = Wcell{i}; %proportions real

    scatter3(W(:,1),W(:,3),W(:,2),[],hsv2rgb(Y),markers{i},'filled','MarkerEdgeColor','k'); hold on
    set(gca, 'Zdir', 'reverse','Xdir', 'reverse');

%     xlabel('Col1 / (Col1 + Col2)');
%     ylabel('B / (B + W)');
%     zlabel('Hue / All');
    axis square
    xlabel('a');
    ylabel('c');
    zlabel('b');

    title(titlecell(spn(i)));
end

%visualize in a cylinder
figure(101);
for i = [2,3,4,1]
    Y = Ycell{i}; %HSV
    alph = 2*pi*Y(:,1); %hue is alpha
    x = cos(alph).*Y(:,2);
    y = sin(alph).*Y(:,2);
    scatter3(x,y,Y(:,3),markersizes{i},hsv2rgb(Y),markers{i},'filled','MarkerEdgeColor','k'); hold on
    axis square
end

%plot a circle
phi = 0:0.01:2*pi;
plot3(cos(phi),sin(phi),0.*phi,'-k');

legend('MY','YC','CM');

xlabel('x');
ylabel('y');
zlabel('V');

return;

%approximation params

eta = 1e-6; %ABM parameter
eps = 1e-6; %elimination parameter

sat = @(x) min(max(x, 0), 1); %saturation function

for dmax = 1:4
    Nexp = 21;
    wp = waitbar(0);
    Ks = linspace(10,50,Nexp);
    errs = zeros(1,Nexp);
    for l = 1:Nexp %number of closest points
        K = Ks(l);
        waitbar(l/Nexp,wp,['Wait... dmax = ',num2str(dmax),', K = ',num2str(K)]);
        errcell = {[],[],[],[]};
        for ctr = 1:4
            [Ncols,~] = size(Ycell{ctr});
            err = zeros(Ncols,1);
            for j = 1:Ncols
                Yactual = Ycell{ctr}; %current color
                Wactual = Wcell{ctr};
                hsvcol = Yactual(j,:);

                i = ctr; %suppose, classification is correct

                %take the corresponding sets
                Y = Ycell{i};  %Y for HSV
                W = Wcell{i}; %W for proportions

                [NY,~ ]= size(Y);

                dst = sqrt(3)*ones(NY,1); %distances - taken as maximum range

                for k = 1:NY
                    if k ~= j && ctr == i % do not include current point if classification is correct
                        %dst(k) = sqrt( (hsvcol - Y(k,:))*(hsvcol - Y(k,:))' ); %Euclidean distance
                        dst(k) = (hsvcol - Y(k,:))*(hsvcol - Y(k,:))'; %squared Euclidean distance
                    end
                end

                [~,I] = sort(dst,1,'ascend');

                %assign weights
                d = 1./dst(I(1:K));
                D = diag(d);

                %take first K points
                %[H, T, ~] = PolyRegression(Ycur,Wcur,0,dmax,eta,eps);
                %props = PolyPredict(hsvcol, H, T);%predicted proportions
                
                O = deglexord(dmax,3); %create order ideal
                [L,~] = size(O);
                
                Ycur = Y(I(1:K),:);

                props = [0 0 0];
                E = EvalPoly(eye(L),Ycur,O);

                for k = 1:3
                    Wcur = W(I(1:K),k);
                    h = (E'*D*E)\(E'*D*Wcur); %WLS
                    props(k) = EvalPoly(h,hsvcol,O); %predict proportions for hsvcol
                end


                props = sat(props);
                err(j) = norm(props - Wactual(j,:));
            end
            errcell{ctr} = err;
        end
        errtotal = [];
        for i = 1:4
            errtotal = [errtotal; errcell{i}];
        end
        errs(l) = mean(errtotal);
    end

    close(wp);

    figure(1); hold on
    plot(Ks,errs);
end

figure(1);
xlabel('$K$','interpreter','latex');
ylabel('err','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
legend('$d_{max} = 1$','$d_{max} = 2$','$d_{max} = 3$','$d_{max} = 4$','interpreter','latex');

