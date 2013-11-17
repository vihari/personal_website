function [Ydot]=SIRode(t,Y, beta,gamma,N)
% [Ydot]=SIRode(t,Y, beta,gamma,N)
% Input:
% t: time
% YInit: Initial number in susceptible, infected, and recovered
S=Y(1);
I=Y(2);

Ydot(1)=-beta*S*I/N+gamma*I;
Ydot(2)=beta*S*I/N-gamma*I;


Ydot=Ydot';
