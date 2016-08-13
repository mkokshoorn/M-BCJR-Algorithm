# M-BCJR Algorithm

Matlab implementation of the M-BCJR Algorithm from the paper:

     Anderson, J. B., & Prlja, A. (2010, October). Turbo equalization and an 
     M-BCJR algorithm for strongly narrowband intersymbol interference. In 
     Information Theory and its Applications (ISITA), 2010 International 
     Symposium on (pp. 261-266). IEEE
     
## Code Example

Define an intersymbol interference (ISI) channel, v and create M_BCJR object.
 
 >\>\>v=[1,0.8,0.3,0.15,0.07];
 
 >\>\>BCJR_dec=M_BCJR_decoder(v); 
 
 >\>\>M_T=length(v); 
 
 
Create test data and define noise power.
 
 >\>\>data=[0,1,0,1];
 
 >\>\>N_0=0.01;
 
 >\>\>data_len=length(data);
 
BPSK modulate data and padd signal for M-BCJR ternmination.
 
 >\>\>tx_sig=[-1 * ones(1,M_T),2*data-1,-1 * ones(1,M_T)];  % Pad with -1 for BCJR algorithm.
 
 >\>\>rx_sig=conv(tx_sig,v) % Add noise here.

 >\>\>y=rx_sig(M_T+1:2*M_T+data_len-1);
 
Decode the recieved signal.

 >\>\>BCJR_dec.step(y,zeros(data_len,1) ,N_0*ones(data_len,1),M)
