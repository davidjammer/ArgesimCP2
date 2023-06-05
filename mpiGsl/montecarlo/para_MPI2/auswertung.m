close all;
clear all;

num_sim1 = [3370840, 5865795, 10854329, 14122004];
N1 = [1, 2, 4, 6];

num_sim2 = [2419100, 4843744, 9775779, 19582790, 38779100, 77444826, 146807592, 277403730, 541733102, 592138641];
N2 = [1, 2, 4, 8, 16, 32, 64, 128, 256, 288];

subplot(2,2,1);
plot(N2, num_sim2./num_sim2(1), 'o--'); hold on;
plot(N2,N2)
xlabel('Anzahl der Prozesse');
ylabel('N_p/N_1');
grid on;
xlim([1 288]);
ylim([1 288]);
legend('Seneca', 'ideal');

subplot(2,2,2);
plot(N2, num_sim2, 'o--'); hold on;
plot(N2, num_sim2(1)*N2);
xlabel('Anzahl der Prozesse');
ylabel('Anzahl der Simulationen');
grid on;
xlim([1 288]);
legend('Seneca', 'ideal');


subplot(2,2,3);
plot(N1, num_sim1./num_sim1(1), 'o--'); hold on;
plot(N1,N1)
xlabel('Anzahl der Prozesse');
ylabel('N_p/N_1');
grid on;
xlim([1 8]);
ylim([1 8]);
legend('Davids Laptop', 'ideal');

subplot(2,2,4);
plot(N1, num_sim1, 'o--'); hold on;
plot(N1, num_sim1(1)*N1);
xlabel('Anzahl der Prozesse');
ylabel('Anzahl der Simulationen');
grid on;
xlim([1 8]);
legend('Seneca', 'ideal');
