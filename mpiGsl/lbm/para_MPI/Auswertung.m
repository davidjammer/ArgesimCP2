clear all;
close all;

%Laufzeiten / Speedup Sceneca
%mpicc -O3 main.c -lm

%parallele Impl.:

%Tasks: 256
%real	1m39,968s
%user	122m35,182s
%sys	12m30,788s


%Tasks: 128
%real	3m9,188s
%user	250m44,834s
%sys	12m0,312s

%Tasks: 64
%
%real	1m27.786s
%user	83m6.947s
%sys	7m35.228s
%Tasks: 32
%
%real	2m0.132s
%user	60m53.493s
%sys	2m3.453s
%Tasks: 16
%
%real	3m4.466s
%user	48m12.823s
%sys	0m28.468s
%Tasks: 8
%
%real	29m23.821s
%user	234m45.401s
%sys	0m10.850s
%Tasks: 4
%
%real	74m20.943s
%user	297m11.468s
%sys	0m3.790s
%Tasks: 2
%
%real	149m27.370s
%user	298m47.527s
%sys	0m1.512s
%Tasks: 1
%
%real	295m16.595s
%user	295m11.247s
%sys	0m1.079s

%Laufzeiten / Speedup David Laptop
%mpicc -O3 main.c -lm

%parallele Impl.:
%Tasks: 1
%real	39m36,514s
%user	39m36,007s
%sys	0m0,081s
%Tasks: 2
%real	12m25,258s
%user	24m49,156s
%sys	0m0,543s
%Tasks: 4
%real	7m39,513s
%user	30m35,906s
%sys	0m0,510s


tp1=[295*60+16, 149*60+27, 74*60+20, 29*60+23, 3*60+4, 2*60, 1*60+27, 3*60+9, 1*60+39];
N1=[1, 2, 4, 8, 16, 32, 64, 128, 256];
s1= tp1(1)./tp1(:);

tp2=[39*60+36, 12*60+25, 7*60+39];
N2=[1, 2, 4];
s2= tp2(1)./tp2(:);

subplot(2,2,1)
plot(N1, N1); hold on;
plot(N1, s1, 'o--');
plot(N2, s2, 'o--');
xlim([1 256]);
ylim([0 250]);
legend('Ideal', 'Seneca');
title('Speed up Lattice Bolzmann');
xlabel('Anzahl  der Prozesse');
ylabel('tp_1/tp_n');
grid on;

subplot(2,2,2)
plot(N1,tp1);
set(gca, 'YScale', 'log')
title('Laufzeit Lattice Bolzmann');
legend('Seneca');
xlim([1 256]);
ylim([10 10^5]);
xlabel('Anzahl  der Prozesse');
ylabel('Laufzeit in s');
grid on;

subplot(2,2,3)
plot(N2, N2); hold on;
plot(N2, s2, 'o--');
xlim([1 8]);
ylim([0 8]);
legend('Ideal', 'Davids Laptop');
title('Speed up Lattice Bolzmann');
xlabel('Anzahl  der Prozesse');
ylabel('tp_1/tp_n');
grid on;

subplot(2,2,4)
plot(N2,tp2);
xlim([1 8]);
ylim([10 10^5]);
set(gca, 'YScale', 'log')
title('Laufzeit Lattice Bolzmann');
legend('Davids Laptop');
xlabel('Anzahl  der Prozesse');
ylabel('Laufzeit in s');
grid on;
