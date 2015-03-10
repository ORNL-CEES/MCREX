N=[100 200 300 500 700];
P=[3 6 12 24 48];

time=[0.402869 4.30778 16.028 78.5828 218.702;
    0.241363 2.88537 8.67008 41.3041 139.926;
    0.496757 2.64127 7.78725 43.5681 126.481;
    0.464905 3.70239 7.43161 32.1722 92.141;
    0.209484 1.84922 5.0329 28.51 76.9692];

loglog(P, time(:,1), '-ro', 'LineWidth', 1, 'MarkerSize', 7);  
hold on
loglog(P, time(:,2), '-bd', 'LineWidth', 1, 'MarkerSize', 7); 
loglog(P, time(:,3), '-g*', 'LineWidth', 1, 'MarkerSize', 7); 
loglog(P, time(:,4), '-k+', 'LineWidth', 1, 'MarkerSize', 7);
loglog(P, time(:,5), '-c*', 'LineWidth', 1, 'MarkerSize', 7);
legend('N=100', 'N=200', 'N=300', 'N=500', 'N=700');
title('MCSA')
xlabel('Number of processes - logarithmic scale');
ylabel('Time (seconds) - logarithmic scale');  


gran=[1/2; 1/4; 1/8; 1/16];
gran_time=[%0.496757 3.70239 7.43161 32.1722 92.141;
    0.476507 2.94749 7.42754 22.6319 73.7607;
    0.542164 2.26947 6.00373 12.1798 45.8213;
    0.889152 1.70438 4.21667 12.0108 29.9483;
    0.627307 1.79536 2.19337 10.3864 19.3599];

hold off
figure()
loglog(gran, gran_time(:,1), '-ro', 'LineWidth', 1, 'MarkerSize', 7);  
hold on
loglog(gran, gran_time(:,2), '-bd', 'LineWidth', 1, 'MarkerSize', 7); 
loglog(gran, gran_time(:,3), '-g*', 'LineWidth', 1, 'MarkerSize', 7); 
loglog(gran, gran_time(:,4), '-k+', 'LineWidth', 1, 'MarkerSize', 7); 
loglog(gran, gran_time(:,5), '-c*', 'LineWidth', 1, 'MarkerSize', 7); 


legend('N=100', 'N=200', 'N=300', 'N=500', 'N=700');
title('MCSA')
xlabel('proportion of the # of histories');
ylabel('Time (seconds) - logarithmic scale');  
