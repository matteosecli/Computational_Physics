x = dlmread('X.txt','\t');
u = dlmread('U.txt','\t');
e = dlmread('E.txt','\t');
save('data.mat');
axes('FontSize', 20);
h = plot(x,u,'.');
title('Solution u(x)', 'FontSize', 20);
xlabel('x','FontSize', 20);
ylabel('y=u(x)','FontSize', 20);
saveas(h, 'solution','eps');
exit
