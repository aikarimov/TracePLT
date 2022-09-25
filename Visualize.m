%process data from files
%backward approximation
%color to proportions
% using closest K points
close all

Hdiv =  0.4; %H value divising magenta-yellow pair into 2 sets

%tablename = 'ModelTable600_very_old.xls';
tablename = 'ModelTable600.xls';
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

legend('red, yellow','yellow, blue','blue, red');

xlabel('$H$','interpreter','latex');
ylabel('$S$','interpreter','latex');
zlabel('$V$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

%visualize in prop system

spn = [1 1 2 3]; %subplot number

titlecell = {'red, yellow','yellow, blue','blue, red'};

for i = 1:4
    
    Y = Ycell{i}; %HSV 
    W = Wcell{i}; %proportions real

    figure;
    scatter3(W(:,1),W(:,2),W(:,3),[],hsv2rgb(Y),markers{i},'filled','MarkerEdgeColor','k'); hold on
    set(gca, 'Zdir', 'reverse','Xdir', 'reverse');
    xlabel('$a$','interpreter','latex');
    ylabel('$b$','interpreter','latex');
    zlabel('$c$','interpreter','latex');
    
    figure(200);
    subplot(1,3,spn(i));
    scatter3(W(:,1),W(:,3),W(:,2),[],hsv2rgb(Y),markers{i},'filled','MarkerEdgeColor','k'); hold on
    set(gca, 'Zdir', 'reverse','Xdir', 'reverse');

    axis square
    xlabel('$a$','interpreter','latex');
    ylabel('$c$','interpreter','latex');
    zlabel('$b$','interpreter','latex');

    title(titlecell(spn(i)));
    set(gca,'TickLabelInterpreter','latex');
end

%visualize in a cylinder
figure(101);
for i = [2,3,4,1]
    Y = Ycell{i}; %HSV
    alph = 2*pi*Y(:,1); %hue is alpha
    x = cos(alph).*Y(:,2);
    y = sin(alph).*Y(:,2);
    scatter3(x,y,Y(:,3),markersizes{i},hsv2rgb(Y),markers{i},'filled','MarkerEdgeColor','k'); hold on
    scatter3(x,y,0*x,7,[0.8 0.8 0.8],'o','filled','MarkerEdgeColor','none'); hold on
    axis square
end

%plot a circle
phi = 0:0.01:2*pi;
plot3(cos(phi),sin(phi),0.*phi,'-k','LineWidth',1);

legend('red, yellow','yellow, blue','blue, red');

xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$V$','interpreter','latex');

xtickformat('$%g$');
ytickformat('$%g$');

set(gca,'TickLabelInterpreter','latex');


