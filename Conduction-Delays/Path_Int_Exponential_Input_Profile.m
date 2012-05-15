%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 1-d Continuous Attractor Neural Network for path integration
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %clear;  %hold on; view(0,90);
   nn = 100; dx=100/nn; % number of node an resolution in deg
   tau_inv = 1./1;      % inverse of membrane time constant
   beta =.07; alpha=.0; % transfer function is 1/(1+exp(-beta*(u-theta))
 %weight matrices
   sig = 5/dx;
   [ws,wa]=hebb_trace_path_or(nn,sig);
   w_inh=7*(sqrt(2*pi)*sig)^2/nn;
   
   %%%% addended to code
   I_in_start=15;
   thickness=9;
   %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%   Experiment    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% external input to initiate bubble, no idiothetic cue, w sym.
   u0 = zeros(nn,1)-10; tspan=[0,10]; param=0;
   w=ws-w_inh;
   I_ext=zeros(nn,1); for i=I_in_start-thickness:I_in_start+thickness; I_ext(i)=50; end
   [t1,u1]=ode45('rnn_ode_or',tspan,u0,[],nn,tau_inv,dx,beta,alpha,w,I_ext);
 %%%% no external input to equilibrate, no idiothetic cue.
   u0 = u1(size(t1,1),:);  tspan=[10,20];
   I_ext=zeros(nn,1);
   [t2,u2]=ode45('rnn_ode_or',tspan,u0,[],nn,tau_inv,dx,beta,alpha,w,I_ext);
%  %%%% clockwise rotation node on --> weight matrix asymetric
%    u0 = u2(size(t2,1),:); tspan=[20,40];
    w=(ws-w_inh).*(1+wa(:,:,2));
%    [t3,u3]=ode45('rnn_ode_or',tspan,u0,[],nn,tau_inv,dx,beta,alpha,w,I_ext);
%  %%%% no idiothetic cue --> weight matrix symetric
%    u0 = u3(size(t3,1),:); tspan=[40,50];
%    w=ws-w_inh; I_ext=zeros(nn,1);
%    [t4,u4]=ode45('rnn_ode_or',tspan,u0,[],nn,tau_inv,dx,beta,alpha,w,I_ext);
%  %%%% counter-clockwise rotation node on --> weight matrix asymetric
%    u0 = u4(size(t4,1),:); tspan=[50,70];
%    w=(ws-w_inh).*(1+2*wa(:,:,1)); I_ext=zeros(nn,1);
%    [t5,u5]=ode45('rnn_ode_or',tspan,u0,[],nn,tau_inv,dx,beta,alpha,w,I_ext);
%  %%%% no idiothetic cue --> weight matrix symetric
%    u0 = u5(size(t5,1),:); tspan=[70,80];
%    w=ws-w_inh; I_ext=zeros(nn,1);
%    [t6,u6]=ode45('rnn_ode_or',tspan,u0,[],nn,tau_inv,dx,beta,alpha,w,I_ext);
%  %%%% plot results
%    u=[u1; u2; u3; u4; u5; u6]; t=[t1; t2; t3; t4; t5; t6];
%    r=1./(1+exp(-beta.*(u-alpha)));
%    %surf(t',1:nn,r','linestyle','none')
%    h=surf(t',1:nn,r');
%    set(h,'linestyle','none');



t1=200;
t2=2000;
% t3=6000*t_f;
% t4=8000*t_f;
% t5=10000*t_f;
w1=ws-w_inh;
% w2=(ws-w_inh).*(1+wa(:,:,2));
% w3=ws-w_inh;
% w4=(ws-w_inh).*(1+2*wa(:,:,1));
% w5=ws-w_inh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clockwise rotation, Euler method
u0 = u2(size(t2,1),:); %% don't mess with this (it's not a trivially distributed, so don't 'manually' construct it.
dt=.01; %%% 'infintesimal' time step
tau_m=1;
t_span=t2; %%% number of time steps
I_ext_time=1000; %%% time duration for which to apply I_ext
delay_r_mu=1; %%% ditto from above except now...
delay_r_sig=0; %%% ditto ...
delay_matrix=round(triu(normrnd(delay_r_mu*ones(nn,nn),delay_r_sig),1));
delay_matrix=delay_matrix+delay_matrix'+eye(nn);
r=max(max(delay_matrix)); %%% radius of the head direction cell system, must be less than or equal to I_ext_time and greater than or equal to 1
if min(min(delay_matrix))<0;
    display('delays out of bounds!!!')
    return
end

%% Extenal Input
I_ext_max=50;
amp=I_in_start;
lower_by=10; %% adjusts (lowers) the maximum rate of change of the input, at the last time step
B=log((nn-thickness-lower_by)/amp)/(t_span-t1);
I_ext=zeros(nn,t_span);
for a=t1+1:t_span;
    for b=1:nn;
        I_ext(b,a)=I_ext_max*(abs(amp*exp(B*(a-t1))-b)<thickness);
    end
end



u_output=zeros(nn,t_span);  %%% matrix holding activity of the nodes
for i=1:r; u_output(:,i)=u0; end;  %%% the first column in the above matrix is u0
%%t1:
for i=r+1:t1  %%% for loop that determines the activity at the rest of the time steps
    for j=1:nn;
        u_out_j=zeros(nn,1);
        for k=1:nn;
            u_out_j(k)=u_output(k,i-delay_matrix(j,k)); %% vector holding the activity from all neurons with delays that neuron j needs to integrate
        end
        w_sum=w1(j,:)*f1_or(u_out_j);
        u_output(j,i)=u_output(j,i-1)*(1-dt*tau_m)+dt*tau_m*(w_sum*dx+I_ext(j,i));
    end
end
%%t2:
for i=t1+1:t2  %%% for loop that determines the activity at the rest of the time steps
    for j=1:nn;
        u_out_j=zeros(nn,1);
        for k=1:nn;
            u_out_j(k)=u_output(k,i-delay_matrix(j,k)); %% vector holding the activity from all neurons with delays that neuron j needs to integrate
        end
        w_sum=w1(j,:)*f1_or(u_out_j);
        u_output(j,i)=u_output(j,i-1)*(1-dt*tau_m)+dt*tau_m*(w_sum*dx+I_ext(j,i));
    end
end

%surf(f1(u_output),'linestyle','none'); view(0,90); %%Hill of Activity, top
%view

active_head_n=f1(u_output)>.9; %%% Gauge the lower limit for active_head_n
%% to a smooth head direction curve
%% Decoder
head_dir=zeros(1,t_span);
for i=1:t_span; 
    n_sum=0;
    for j=1:nn;
        if active_head_n(j,i)==1;
            n_sum=n_sum+j;
        end
    end
    head_dir(i)=n_sum/sum(active_head_n(:,i));
end

%% Decodes input so it can be plotted alongside input
input_dir=zeros(1,t_span);
for i=1:t_span; 
    n_sum=0;
    for j=1:nn;
        if I_ext(j,i)==I_ext_max;
            n_sum=n_sum+j;
        end
    end
    input_dir(i)=n_sum/(sum(I_ext(:,i))/I_ext_max);
end

plot(1:t_span,head_dir,1:t_span,input_dir)

