function [Ydot]=SIRode(t,Y, beta,gamma,N)
% [Ydot]=SIRode(t,Y, beta,gamma,N)
% Input:
% t: time
% YInit: Initial number in susceptible, infected, and recovered
S=Y(1);
I=Y(2);
R=Y(3);

Ydot(1)=-beta*S*I/N;
Ydot(2)=beta*S*I/N-gamma*I;
Ydot(3)=gamma*I;

Ydot=Ydot';
