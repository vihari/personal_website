function [Ydot]=SEIode(t,Y, beta,gamma,p,r,N)
% [Ydot]=SIRode(t,Y, beta,gamma,N)
% Input:
% t: time
% YInit: Initial number in susceptible, infected, and recovered
S=Y(1);
L=Y(2);
I=Y(3);

Ydot(1)=-beta*S*I/N+gamma*I+gamma*L;
Ydot(2)=(1-p)*beta*S*I/N-r*L-gamma*L;
Ydot(3)=p*beta*S*I/N+r*L-gamma*I;

Ydot=Ydot';
