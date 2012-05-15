 function udot=rnn(t,u,flag,nn,tau_inv,dx,beta,alpha,w,I_ext)
 % odefile for recurrent network
   r=1./(1+exp(-beta.*(u-alpha)));
   sum=w*r*dx;
   udot=tau_inv*(-u+sum+I_ext);
 return