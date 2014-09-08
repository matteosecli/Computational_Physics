errors = dlmread('errors.txt','\t');
save('errors.mat');
axes('FontSize', 20);
h = loglog(errors(:,1),errors(:,2),'-s');
title('Percentage error of the solution', 'FontSize', 24);
grid on;
xlabel('N','FontSize', 20);
ylabel('Error [%]','FontSize', 20);
saveas(h, 'errors','eps');
exit
