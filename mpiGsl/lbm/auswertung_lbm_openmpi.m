close all;

%%%2048x2048 500000 Iter.
t_seq = (29*60+10.023) * 500000/1000;
t_para = [  (29*60+30.869) * 500000/1000,...    %N=1
            (16*60+19.939)  * 500000/1000,...    %N=2
            (17*60+48.529)  * 500000/2000,...    %N=4
            (16*60+13.824)  * 500000/4000,...    %N=8
            (14*60+56.617) * 500000/8000,...    %N=16
            (40*60+11.280) * 500000/20000,...   %N=32
            (42*60+28.629) * 500000/40000,...   %N=64
            (43*60+20.450)  * 500000/80000,...   %N=128
            (125*60+43.013)                     %N=256
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
legend("C OpenMPI","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "C OpenMPI","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")

%%%1024x1024 500000 Iter.
t_seq = (8*60+42.847) * 500000/1000;
t_para = [  (8*60+35.681) * 500000/1000,...    %N=1
            (18*60+57.738) * 500000/4000,...   %N=2
            (21*60+25.034) * 500000/10000,...  %N=4
            (36*60+58.343) * 500000/20000,...  %N=8
            (42*60+31.684) * 500000/50000,...  %N=16
            (40*60+55.567) * 500000/100000,...  %N=32
            (98*60+6.094),...                 %N=64
            (57*60+15.139),...                 %N=128
            (12*60+58.099)                     %N=256
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
legend("C OpenMPI","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "C OpenMPI","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")

%%%512x512 350000 Iter.
t_seq = (10*60+56.713) * 350000/5000;
t_para = [  (12*60+56.810)  * 350000/5000,...   %N=1
            (17*60+52.682) * 350000/10000,...   %N=2
            (15*60+29.504) * 350000/20000,...   %N=4
            (35*60+55.926) * 350000/100000,...  %N=8
            (32*60+9.376) * 350000/100000,...   %N=16
            (14*60+48.730) * 350000/200000,...   %N=32            
            (2*60+7.196),...                   %N=64
            (3*60+11.471),...                   %N=128
            (2*60+12.915)                       %N=256
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
legend("C OpenMPI","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "C OpenMPI","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")


%%%256x256 350000 Iter.
t_seq = (3*60+51.618) * 350000/10000;
t_para = [  (7*60+50.318)  * 350000/10000,...   %N=1
            (3*60+20.782) * 350000/20000,...    %N=2
            (5*60+3.643) * 350000/40000,...    %N=4
            (1*60+32.630) * 350000/100000,...   %N=8
            (1*60+42.204) * 350000/200000,...   %N=16
            (1*60+55.194),...                   %N=32            
            (1*60+18.945),...                   %N=64
            (2*60+6.776),...                   %N=128
            (1*60+22.603)                       %N=256
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
legend("C OpenMPI","Location","best");
title("Runtime")

subplot(2,1,2)
plot(N,N); hold on;
plot(N,speedup,'o--'); hold off;
legend("perfect", "C OpenMPI","Location","best");
grid on;
xlabel("Number of Cores")
ylabel("ts/tp_n")
title("Speedup")