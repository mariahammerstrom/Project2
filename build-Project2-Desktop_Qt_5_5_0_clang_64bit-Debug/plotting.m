% Plot data from C++

n = 200;

% Omega 5.00:
omega = 5.00;
filename = sprintf('Eigenvectors_%0.6f.txt',omega);
[x,y] = textread(filename,'%f %f',n);

% Omega 1.00:
omega = 1.00;
filename = sprintf('Eigenvectors_%0.6f.txt',omega);
[x2,y2] = textread(filename,'%f %f',n);

% Omega 0.50:
omega = 0.50;
filename = sprintf('Eigenvectors_%0.6f.txt',omega);
[x3,y3] = textread(filename,'%f %f',n);

% Omega 0.01:
omega = 0.01;
filename = sprintf('Eigenvectors_%0.6f.txt',omega);
[x4,y4] = textread(filename,'%f %f',n);

% PLOTTING:
figure;
plot(x,y,x2,y2,x3,y3,x4,y4);
xlabel('\rho','FontSize', 18);
ylabel('\psi','FontSize', 18);
title(sprintf('Wave function, n = %d', n));
legend('$\omega =5.00$','$\omega=1.00$','$\omega=0.5$','$\omega=0.01$','Location','southeast');
set(legend,'FontSize',14,'Interpreter','latex');
print(sprintf('wave_function.png'),'-dpng');



%%%%%%%%%%%%%%%%

% Omega 5.00:
omega = 5.00;
filename = sprintf('Eigenvectors2_%0.6f.txt',omega);
[x,y] = textread(filename,'%f %f',n);

% Omega 1.00:
omega = 1.00;
filename = sprintf('Eigenvectors2_%0.6f.txt',omega);
[x2,y2] = textread(filename,'%f %f',n);

% Omega 0.50:
omega = 0.50;
filename = sprintf('Eigenvectors2_%0.6f.txt',omega);
[x3,y3] = textread(filename,'%f %f',n);

% Omega 0.01:
omega = 0.01;
filename = sprintf('Eigenvectors2_%0.6f.txt',omega);
[x4,y4] = textread(filename,'%f %f',n);

% PLOTTING:
figure;
plot(x,normc(y),x2,normc(y2),x3,normc(y3),x4,normc(y4)); % normalized
%plot(x,y,x2,y2,x3,y3,x4,y4);
xlabel('\rho','FontSize', 18);
ylabel('$|\psi|^2$','FontSize', 18,'Interpreter','latex');
title(sprintf('Probability distribution, n = %d', n));
legend('$\omega =5.00$','$\omega=1.00$','$\omega=0.5$','$\omega=0.01$','Location','northeast');
set(legend,'FontSize',14,'Interpreter','latex');
print(sprintf('probability_distribution.png'),'-dpng');

figure;
n = [5,10,20,50,100,200];
%i = [28,42,676,4378,17789,72044];
i = [11,41,176,1727,11315,61308];
plot(n,i,n,n.^(2.1));
xlabel('n');
ylabel('no. of iterations');
legend('numerical','approximation');