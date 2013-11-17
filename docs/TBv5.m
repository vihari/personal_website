clear ALL
clc
num_members=100;
numstates=1;
num_iterations=24;
diff_num_states = 35;
betalast = [];
plast = [];
Lcountry = [];
Icountry = [];
sumrate = [1:num_iterations]-[1:num_iterations];
suminfection = [1:num_iterations]-[1:num_iterations];
sumI = [1:num_iterations]-[1:num_iterations];
sumL = [1:num_iterations]-[1:num_iterations];
WHOcmp = 0;
under_estimation = 1.;
store_beta_state = [];
store_p_state = [];
max_beta = [];
max_p = [];
max_state = 0;
max_state_p = 0;
min_beta = [];
min_p = [];
min_state = 0;
min_state_p = 0;
product_p_beta = [];
pstore = [];
betastore = [];
store_beta_country = [];
store_p_country = [];
RMSE = [];
beta_2012 = [];
p_2012 = [];
L_2012 = [];
I_2012 = [];
for state_iterator=1:diff_num_states
Lstate = [];
Istate = [];
Lstateweighted = [];
Istateweighted = [];
sp= under_estimation*Sp(state_iterator,:);
spn= under_estimation*Ln(state_iterator,:);
spe= under_estimation*Le(state_iterator,:);

time = 0:1:9;
beta = 5;

gamma = 0.85;
N = P(state_iterator,:)*100000;
x=sp(5);
y=spn(5)+spe(5);
p=0.15;
r=0.02;


for i=1:1
    w(i,1)=0.5;
    w(i,2)=0.02;
end

a_estbar(:,1)=5;
p_estbar(:,1)=p;

for j=1:num_members
for i=1:4
   ahist(i,j)=5;
   phist(i,j)=p;
end
end

for j=1:num_members
     for k=1:numstates
      a_est(k,j)=a_estbar(k,1);
      p_est(k,j)=p_estbar(k,1);
     end
end
    
Zcov=eye(numstates); %create measurement noise covariance matrix

YInit=[N-x-y y x];

options=[];
tt=[0:1:2];

  
  PredY=ones(3,3);
  S=ones(num_iterations,1);
  L=ones(num_iterations,1);
  I=ones(num_iterations,1);
  Spred=ones(num_iterations,1);
  Lpred=ones(num_iterations,1);
  Ipred=ones(num_iterations,1);
  Lassim=ones(num_iterations,1);
  Iassim=ones(num_iterations,1);
  
  S(1,1)=N-x-y;
  L(1,1)=y;
  I(1,1)=x;
 means(1) = sp(1);
 plotbeta = [];
 plotp = [];
 index=0;
 for i=1:num_iterations
   
     index=index+1;
     if (index>4)
         index=1; 
     end 
     

     if(i==1)
         store_beta_state = [];
         store_p_state = [];
     end
     % mesurement noise is 5%
    for m=1:numstates
     z(m,1)=0.05*L(i,1);
     z(m,2)=0.05*I(i,1);
     Zcov(m,m)=z(m,1)^2;
     Zcov(m+1,m+1)=z(m,1)^2;
    end
    
    for j=1:num_members
     W(:,j)=w(1,1).*randn(numstates,1);      
     WW(:,j)=w(1,2).*randn(numstates,1);      
     
     Z(:,j)=z(1,1).*randn(numstates,1);
     ZZ(:,j)=z(1,2).*randn(numstates,1);
       
     for m=1:numstates
         a_est(m,j)=ahist(index,j)+W(m,j);             
         p_est(m,j)=phist(index,j)+WW(m,j);          
%     a_est(m,j)=a_est(m,j)+W(m,j);             
%     p_est(m,j)=p_est(m,j)+WW(m,j); 
            
         
         y(m,j)= spn(i)+spe(i)+Z(m,j);    
         y1(m,j)= sp(i)+ZZ(m,j) ;       

         y_f1(m,j)=a_est(m,j)*(1-p_est(m,j))*S(i,1)*I(i,1)/N-r*L(i,1);
         y_f2(m,j)=a_est(m,j)*p_est(m,j)*S(i,1)*I(i,1)/N+r*L(i,1);
     end 
    end
    
     a_estbar(:,i) = mean(a_est,2) ;
     p_estbar(:,i) = mean(p_est,2) ;
     
     beta=a_estbar(:,i);
     p=p_estbar(:,i);
     Lpred(i) = (1-p)*(beta*S(i,1)*I(i,1)/N-r*L(i,1));
     Ipred(i) = p*(beta*S(i,1)*I(i,1)/N+r*L(i,1));
     
     y_f1bar=mean(y_f1,2);
     y_f2bar=mean(y_f2,2);
     ybar=mean(y,2);
     y1bar=mean(y1,2);
 
   
for j=1:numstates
    for k=1:num_members
     Ex(j,k)=[a_est(j,k)-a_estbar(j,i)];
     Ex(numstates+j,k)=[p_est(j,k)-p_estbar(j,i)];
   end
end

for j=1:numstates
    for k=1:num_members
     Ey(j,k)=[y_f1(j,k)-y_f1bar(j)];
     Ey(numstates+j,k)=[y_f2(j,k)-y_f2bar(j)];
    end
end

   Pxy=Ex*Ey'/(num_members-1);
   Pyy=Ey*Ey'/(num_members-1)+Zcov;
   K=Pxy*inv(Pyy);
   
 %plot(L)
 means(i+1) = y_f1bar+y_f2bar;

   for j=1:numstates
    for k=1:num_members
     inov(j,k)=y(j,k)-y_f1(j,k);
     inov(numstates+j,k)=y1(j,k)-y_f2(j,k);
    end
   end

  gain=K*inov;
 
   for j=1:numstates
    for k=1:num_members
     a_est(j,k)=a_est(j,k)+gain(j,k);
     p_est(j,k)=p_est(j,k)+gain(numstates+j,k);
     ahist(index,k)=a_est(j,k);
     phist(index,k)=p_est(j,k);
    end
   end
     
   a_estbar(:,i)= mean(a_est,2);
   p_estbar(:,i)= mean(p_est,2) ;
   p=p_estbar(:,i);
   beta=a_estbar(:,i);
   store = beta;
  
%    plotbeta = [plotbeta;beta];
%    plotp = [plotp;p];

%    if(i>4)
%     Lpred(i+1) = (1-store_p_state(i-4))*(store_beta_state(i-4)*S(i,1)*I(i,1)/N-r*L(i,1));
%     Ipred(i+1) = store_p_state(i-4)*(store_beta_state(i-4)*S(i,1)*I(i,1)/N+r*L(i,1));
%    else


Lassim(i) = (1-p)*(beta*S(i,1)*I(i,1)/N-r*L(i,1));
Iassim(i) = p*(beta*S(i,1)*I(i,1)/N+r*L(i,1));

%    end
   
[t,PredY]=ode15s(@SEIode, tt, YInit, options, beta, gamma,p,r, N);

S(i+1,1)=PredY(2,1);
L(i+1,1)=PredY(2,2);
I(i+1,1)=PredY(2,3);

if(i==num_iterations-1)
    I_2012 = [I_2012;Iassim(i)];
    L_2012 = [L_2012;Lassim(i)];
end

% Lpred(i+1) = L(i+1,1);
% Ipred(i+1) = I(i+1,1);

YInit=[S(i+1,1) L(i+1,1) I(i+1,1)];

Lstate = [Lstate,L(i+1,1)];
Istate = [Istate,I(i+1,1)];

if(det(Pyy)>1)
    Lstateweighted = [Lstateweighted,store*L(i+1,1)];
    Istateweighted = [Istateweighted,store*I(i+1,1)];
else
    Lstateweighted = [Lstateweighted,1*L(i+1,1)];
    Istateweighted = [Istateweighted,1*I(i+1,1)];
end

store;
store_beta_state = [store_beta_state,store];
store_p_state = [store_p_state,p];
if(i == num_iterations)
    beta_2012 = [beta_2012;store];
    p_2012 = [p_2012;p];
end
if(i == num_iterations)
    store_beta_state;
    if(det(Pyy)>1)%Or else the values are not reliable.
        store_beta_country = [store_beta_country;store_beta_state];
        store_p_country = [store_p_country;store_p_state];
        product_p_beta = [product_p_beta,mean(store_beta_state.*store_p_state)];
        betastore = [betastore,mean(store_beta_state)];
        pstore = [pstore,mean(store_p_state)];
    else
        betastore = [betastore,1.5];%Some mean values
        pstore = [pstore,0.5];
    end
    if state_iterator == 1
        max_beta = store_beta_state;
        max_state = 1;
        min_beta = store_beta_state;
        min_state = 1;
        max_p = store_p_state;
        max_state_p = 1;
        min_p = store_p_state;
        min_state_p = 1;
    else
        %mean from only after 5 steps is considered, so that erroneous
        %initial values of beta and p does not play a role 
        if mean(max_beta(5:num_iterations-1))<mean(store_beta_state(5:num_iterations-1))
            if(det(Pyy)>1)
                max_beta = store_beta_state;
                max_state = state_iterator;
            end
        
        elseif mean(min_beta(5:num_iterations-1))>mean(store_beta_state(5:num_iterations-1))
            if(det(Pyy)>1)
                min_beta = store_beta_state;
                min_state = state_iterator;
            end
        end
        if mean(max_p(5:num_iterations-1))<mean(store_p_state(5:num_iterations-1))
            if(det(Pyy)>1)
                max_p = store_p_state;
                max_state_p = state_iterator;
            end
        
        elseif mean(min_p(5:num_iterations-1))>mean(store_p_state(5:num_iterations-1))
            if(det(Pyy)>1)
                min_p = store_p_state;
                min_state_p = state_iterator;
            end
        end
    end
end

 end

% plot(Lpred,'^');
% hold on
% plot(spn+spe,'*')
% return;
 
 suminfection  = suminfection + Lstate+Istate;
 sumrate = sumrate + Lstateweighted+Istateweighted;
 sumI = sumI+Istate;
 sumL = sumL+Lstate;
Lcountry = [Lcountry;Lstate];
Icountry = [Icountry;Istate];
% plot(plotbeta,'--rs','LineWidth',2,...
%                 'MarkerEdgeColor','k',...
%                  'MarkerFaceColor','g',...
%                     'MarkerSize',10)
% legend(sprintf('Infection Rate when Under estimation = 1 \n when no under estimation is considered'))
% xlabel('Time (in quarter)')
% ylabel(sprintf('Number of people to which the \ndisease is spread\n by a single infected individual \nper quarter'))
% while(1)
% w = waitforbuttonpress;
% if w~=0
% break;
% end
% end

%        hold on;
%plot(sp+spn+spe,'--rs','LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor','b',...
%             'MarkerSize',10);
%     m= sprintf('Ratio of number of\nsmear Positive Patients');
%     legend(m);
%       %      return
%  %figure
% %plot(time, PredI)
% %plot(t,PredS, t, PredI)
% %plot(t, PredI, t, sp)
% xlabel('Time (in quarter)')
% ylabel('Ratio')
% while(1)
% w = waitforbuttonpress;
% if w~=0
% break;
% end
% end
%axis([0 40 0 PredI])
%legend('I')
rmse = 0;
for i=0:23
    if((det(Pyy)>1)&&(Le(state_iterator,i+1) + Ln(state_iterator,i+1)>0)&&(Sp(state_iterator,i+1)>0))
        rm = L(i+1,1) - Le(state_iterator,i+1) - Ln(state_iterator,i+1);
        rm = rm/(Le(state_iterator,i+1) + Ln(state_iterator,i+1));
        rm = rm*rm;
        se = I(i+1,1) - Sp(state_iterator,i+1);
        se = se/(Sp(state_iterator,i+1));
        se = se*se;
        rmse = rmse+se+rm;
    end
end
RMSE = [RMSE;sqrt(rmse/24)];
if(strcmp(states(state_iterator),'Manipur'))
    LstateManipur = Lpred;
    IstateManipur = Ipred;
    LAssimMani = Lassim
    IAssimMani = Iassim
    Mani = state_iterator;
    ActualLMani = spn(1:24)+spe(1:24)
    ActualIMani = sp(1:24)
end
if(strcmp(states(state_iterator),'Pondicherry'))
    LstatePondicherry = Lpred;
    IstatePondicherry = Ipred;
    LAssimPond = Lassim
    IAssimPond = Iassim
    Pond = state_iterator;
    ActualLPond = spn+spe
    ActualIPond = sp
end
% WHOcmp = WHOcmp+sum(Lpred(21:24));
 WHOcmp = WHOcmp+sum(Ipred(21:24));
end



betacountry = sumrate./suminfection;
plot(max_beta,'^:','LineWidth',2,...
                  'MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                   'MarkerSize',10)
   hold on;
  plot(betacountry,'o-','LineWidth',2,...
                  'MarkerEdgeColor','b',...
                  'MarkerFaceColor','b',...
                   'MarkerSize',10)
               hold on;
plot(min_beta,'v--','LineWidth',2,...
                  'MarkerEdgeColor','g',...
                  'MarkerFaceColor','g',...
                   'MarkerSize',10)
               hold on;
 max_state_name = char(states(max_state));
 min_state_name = char(states(min_state));
 legend(max_state_name,'India',min_state_name);

 text(12.5,0.25,...
'Time(Yr)',...
'EdgeColor','red');

 w = waitforbuttonpress;
 pcountry = Istate./(Lstate+Istate);
 
 plot(Sp(Mani,:),'o:','LineWidth',2,...
                  'MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                   'MarkerSize',10)
               hold on;
 plot(pcountry,'o-','LineWidth',2,...
                  'MarkerEdgeColor','b',...
                  'MarkerFaceColor','b',...
                   'MarkerSize',10)
               hold on;
 plot(min_p,'v--','LineWidth',2,...
                  'MarkerEdgeColor','g',...
                  'MarkerFaceColor','g',...
                   'MarkerSize',10)
               hold on;
 max_state_name = char(states(max_state_p));
 min_state_name = char(states(min_state_p));
 legend(max_state_name,'India',min_state_name);

 text(12.5,0.2,...
    'Time(Yr)',...
    'EdgeColor','red');

%plot(Le(Mani,:)+Ln(Mani,:),'r')
%plot(LstateManipur,'b')

%% Graphing part for L for state:Manipur
plot(ActualLMani,'^--','LineWidth',2,...
                  'MarkerEdgeColor','b',...
                  'MarkerFaceColor','b',...
                   'MarkerSize',10)
hold on;
plot(LstateManipur,'o--','LineWidth',2,...
                  'MarkerEdgeColor','k',...
                  'MarkerFaceColor','k',...
                   'MarkerSize',10)
hold on;
plot(LAssimMani,'v--','LineWidth',2,...
                  'MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                   'MarkerSize',10)
hold on;
% plot(LAssimManiRaw,'x-','LineWidth',2,...
%                   'MarkerEdgeColor','g',...
%                   'MarkerFaceColor','g',...
%                    'MarkerSize',10)
% hold on;
% plot(LstateManipurRaw,'d-','LineWidth',2,...
%                   'MarkerEdgeColor','m',...
%                   'MarkerFaceColor','m',...
%                    'MarkerSize',10)

% legend('O','P','A','A\prime','P\prime');
legend('Observed','Predicted','Assimilated');

%% Graphing part for I for state:Manipur
plot(ActualIMani,'^--','LineWidth',2,...
                  'MarkerEdgeColor','b',...
                  'MarkerFaceColor','b',...
                   'MarkerSize',10)
hold on;
plot(IstateManipur,'o--','LineWidth',2,...
                  'MarkerEdgeColor','k',...
                  'MarkerFaceColor','k',...
                   'MarkerSize',10)
hold on;
plot(IAssimMani,'v--','LineWidth',2,...
                  'MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                   'MarkerSize',10)
hold on;
% plot(IAssimManiRaw,'x-','LineWidth',2,...
%                   'MarkerEdgeColor','g',...
%                   'MarkerFaceColor','g',...
%                    'MarkerSize',10)
% hold on;
% plot(IstateManipurRaw,'d-','LineWidth',2,...
%                   'MarkerEdgeColor','m',...
%                   'MarkerFaceColor','m',...
%                    'MarkerSize',10)

               

legend('Observed','Predicted','Assimilated');

%%
plot(Le(Pond,:)+Ln(Pond,:),'^:','LineWidth',2,...
                  'MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                   'MarkerSize',10)
               
   hold on;
plot(LstatePondicherry,'v--','LineWidth',2,...
                  'MarkerEdgeColor','b',...
                  'MarkerFaceColor','b',...
                   'MarkerSize',10)


               
               
