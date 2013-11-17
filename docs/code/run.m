time=[0:.1:30];
beta= 2.5;
gamma = 0.85;
N=10000;

YInit=[N 10 0];

options=[];
[t,PredY]=ode15s(@Iode, time, YInit, options, beta, gamma, N);

PredS=PredY(:,1);
PredI=PredY(:,2);
PredR=PredY(:,3);

figure
plot(t,PredS, t, PredI, t, PredR)
xlabel('Time (in Days)')
ylabel('Number Infected')
axis([0 30 0 N])
legend('S', 'I', 'R')

