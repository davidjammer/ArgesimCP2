close all;

%%%256x256 350000 Iter.
t_seq = (60.739072) * 350000/10000;
t_para = [  (61.956733)  * 350000/10000,...     %N=1
            (50.598090) * 350000/20000,...      %N=2
            (60.298954) * 350000/40000,...      %N=4
            (76.889463) * 350000/80000,...      %N=8
            (124.292929) * 350000/160000,...    %N=16
            (313.240020),...                    %N=32            
            (446.716659),...                    %N=64
            (2032.421425),...                   %N=128
            (215.326372) * 350000/10000         %N=256
            ];

speedup = t_seq./t_para;
N = [1,2,4,8,16,32,64,128,256];

figure("Name","256x256")
subplot(2,1,1)
plot(N,t_para,'o--');
set(gca, 'YScale', 'log');
grid on;
xlabel("Number of Cores")
ylabel("Runtime in s")
legend("Matlab PCT","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "Matlab PCT","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")

%%%512x512 350000 Iter.
t_seq = (24.517492) * 350000/1000;
t_para = [  (30.094460)  * 350000/1000,...     %N=1
            (138.203450) * 350000/10000,...      %N=2
            (122.084783) * 350000/20000,...      %N=4
            (191.561294) * 350000/40000,...      %N=8
            (144.207286) * 350000/80000,...    %N=16
            (266.774283) * 350000/160000,...                    %N=32            
            (650.103636),...                    %N=64
            (2372.990044),...                 %N=128
            (5935.349036)                   %256
            ];

speedup = t_seq./t_para;
N = [1,2,4,8,16,32,64,128,256];

figure("Name","512x512")
subplot(2,1,1)
plot(N,t_para,'o--');
set(gca, 'YScale', 'log');
grid on;
xlabel("Number of Cores")
ylabel("Runtime in s")
legend("Matlab PCT","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "Matlab PCT","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")

%%%1024x1024 500000 Iter.
t_seq = (100.973118) * 500000/1000;
t_para = [  (140.767552)  * 500000/1000,...     %N=1
            (112.627005) * 500000/2000,...      %N=2
            (146.416393) * 500000/4000,...      %N=4
            (134.497528) * 500000/8000,...      %N=8
            (232.883495) * 500000/16000,...     %N=16
            (407.914700) * 500000/40000,...     %N=32            
            (971.241087) * 500000/100000,...    %N=64
            (1015.682545) * 500000/100000,...   %N=128
            (1996.835212) * 500000/100000       %N=256
            ];

speedup = t_seq./t_para;
N = [1,2,4,8,16,32,64,128,256];

figure("Name","1024x1024")
subplot(2,1,1)
plot(N,t_para,'o--');
set(gca, 'YScale', 'log');
grid on;
xlabel("Number of Cores")
ylabel("Runtime in s")
legend("Matlab PCT","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "Matlab PCT","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")


%%%2048x2048 500000 Iter.
t_seq = (445.377815) * 500000/1000;
t_para = [  (521.713920)  * 500000/1000,...     %N=1
            (577.097167) * 500000/2000,...      %N=2
            (589.238048) * 500000/4000,...      %N=4
            (790.592359) * 500000/8000,...      %N=8
            (1263.949360) * 500000/16000,...     %N=16
            (2241.615313) * 500000/32000,...     %N=32            
            (4202.342240) * 500000/64000,...     %N=64
            (2928.864765) * 500000/80000,...      %N=128
            (3528.173159) * 500000/80000        %N=256
            ];

speedup = t_seq./t_para;
N = [1,2,4,8,16,32,64,128,256];

figure("Name","2048x2048")
subplot(2,1,1)
plot(N,t_para,'o--');
set(gca, 'YScale', 'log');
grid on;
xlabel("Number of Cores")
ylabel("Runtime in s")
legend("Matlab PCT","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "Matlab PCT","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")