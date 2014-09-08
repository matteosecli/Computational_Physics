errors = dlmread('errors.txt','\t');
save('errors.mat');
axes('FontSize', 20);
h = loglog(errors(:,1),errors(:,2),'-s');
title('Percentual error of the solution', 'FontSize', 24);
xlabel('N','FontSize', 20);
ylabel('Error [%]','FontSize', 20);
saveas(h, 'errors','eps');
exit
