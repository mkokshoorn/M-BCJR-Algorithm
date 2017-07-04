# M-BCJR Algorithm

A soft input, soft-output MATLAB implementation of the M-BCJR Algorithm from the paper:

     Anderson, J. B., & Prlja, A. (2010, October). Turbo equalization and an 
     M-BCJR algorithm for strongly narrowband intersymbol interference. In 
     Information Theory and its Applications (ISITA), 2010 International 
     Symposium on (pp. 261-266). IEEE
     
## Class Inputs/Outputs     

CONSTRUCTOR:

function obj = M_BCJR_decoder(v)

@input 'v' ISI channel of length M_T taps.

STEP: 

[a_APP_LLR] = step(obj,y,A_ext_LLR,N_0,M,SO)

@input 'y' recieved symbols.

@input 'a_ext_LLR' Extrinsic LLR information each symbol.

@input 'N_0' Noise information for each symbol.

@input 'M' number of survivors at each trellis step.


@output 'a_APP_LLR' APP LLR of each symbol (i.e. soft output)
     
## Code Example

Define an inter-symbol interference (ISI) channel, v and create M_BCJR object.
 
 >\>\>v=[1,0.8,0.3,0.15,0.07];
 
 >\>\>BCJR_dec=M_BCJR_decoder(v); 
 
 >\>\>M_T=length(v); 
 
 
Create test data and define noise power.
 
 >\>\>d=[0,1,0,1]; % Short binary  example sequence.
 
 >\>\>N_0=0.01;
 
 >\>\>data_len=length(d);
 
BPSK modulate data and pad signal for M-BCJR termination.
 
 >\>\>x =[ -1 * ones(1,M_T),2*data-1,-1 * ones(1,M_T)];  % Pad with -1 for BCJR algorithm.
 
 >\>\>rx_sig = conv(x,v) + randn(size(conv(x,v))).*sqrt(N_0)% Add noise here.

 >\>\>y=rx_sig(M_T+1:2*M_T+data_len-1);
 
Decode the received signal.

 >\>\> x_prior_LLR = zeros(data_len,1); % Assume 1/-1 with equal prob.  
 >\>\> x_LLR = BCJR_dec.step(y, ,N_0*ones(data_len,1),M)
 
x_LLR =

   -10    10   -10    10
 
 
  ## Turbo Decoder 
  
  In a turbo decoder configuration (coupled with an outer code), x_prior_LLR can be 
  updated on each iteration by feeding back the extrinsic information of the codeword. 
  
 
 
 
