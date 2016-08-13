classdef M_BCJR_decoder

% -------------------------------------------------------------------------
% Author:        Matthew L Kokshoorn 
% Date:          10/10/2015   
% Organization:  The University of Sydney
% Position:      PhD Candidate
% -------------------------------------------------------------------------
% CLASS: Implementation of the M-BCJR Algorithm from the paper [1]:
% [1] Anderson, J. B., & Prlja, A. (2010, October). Turbo equalization and an 
%     M-BCJR algorithm for strongly narrowband intersymbol interference. In 
%     Information Theory and its Applications (ISITA), 2010 International 
%     Symposium on (pp. 261-266). IEEE
% ------------------------------------------------------------------------%
% CONSTRUCTOR:
% @input 'v' ISI channel of length M_T taps.
%
% STEP:
% @input 'y' recieved .
% @input 'A_ext_LLR' APP information for each symbol.
% @input 'N_0' Noise information for each symbol.
% @input 'M' number of survivors at each trellis step.
%
% @output 'SO' Extrinsic LLR of each symbol.
%
% EXAMPLE:
% 
% 
% v=[1,0.8,0.3,0.15,0.07];
% BCJR_dec=M_BCJR_decoder(v);
% M_T=length(v)
% 
% %Create Data
% data=[0,1,0,1];
% data_len=length(data);
% 
% %%BPSK tranmission.
% tx_sig=[-1*ones(1,M_T),2*data-1,-1*ones(1,M_T)];  %Pad with -1 for BCJR algorithm.
% rx_sig=conv(tx_sig,v)  % + Noise.
% 
% 
% y=rx_sig(M_T+1:2*M_T+data_len-1)
% 
% N_0=0.01;
% BCJR_dec.step(y,zeros(data_len,1) ,N_0*ones(data_len,1),M)

% -------------------------------------------------------------------------
% Notes: Assumes that intiial state is the all '-1'
% -------------------------------------------------------------------------    
    
    properties
        
        M_T;
        R;
        C;
        
        l;
        S_in;
        
        pos;
        neg;
               
        S_alpha_set;
        S_beta_set;
        l_beta_out;
    end
        
        
    methods
        
        function obj = M_BCJR_decoder(v)
            
            %Perpare Variables
            obj.M_T=length(v);
            A_n=npermutek([-1,1],obj.M_T);
            [obj.R,obj.C]=size(A_n);
           
            %Generate Trellis -----------------------------------------------------
            obj.l=zeros(obj.R,obj.R);
            obj.S_in=zeros(obj.R,obj.R);
            for state_num=1:obj.R
                for input_sym=[-1,1]

                    %Calculate New State given input symbol.
                    new_state=[A_n(state_num,2:end),input_sym];
                    new_state_index= find(sum(A_n==ones(obj.R,1)*new_state,2)==(obj.C));

                    %Calculate output given new state.
                    output_sym=0;
                    for n=1:length(new_state)
                        output_sym =output_sym+ A_n(new_state_index,n)*v(end-n+1);
                    end

                    %Place in trellis.
                    obj.l(state_num,new_state_index)=output_sym;
                    obj.S_in(state_num,new_state_index)=input_sym;

                end
            end
            
            obj.pos=(A_n(:,end)==(1) );
            obj.neg=(A_n(:,end)==(-1));
            
            
            %Create Set - All states that can be reached from state I.
            obj.S_alpha_set=zeros(obj.R,2);              
            for I=1:obj.R
                obj.S_alpha_set(I,:)=find(obj.S_in(I,:));
            end
            
            
            %Create Set - All states that can reach state J.
            obj.S_beta_set=zeros(obj.R,2);
            obj.l_beta_out=zeros(obj.R,2);
            for J=1:obj.R
                obj.S_beta_set(J,:)=find(obj.S_in(:,J));
                obj.l_beta_out(J,:)=obj.l(abs(obj.l(:,J))>0,J);
            end        
            

                           
        end

        function [SO] = step(obj,y, A_ext_LLR ,N_0,M)
            
            %Prepare Extrinsic Info. --------------------------------------
            A_ext_LLR_appended=[A_ext_LLR(:)',0*ones(1,obj.M_T-1)];
            N_0=[N_0;mean(N_0).*ones(obj.M_T,1)];
            pi_noise=1./sqrt(pi*N_0);

            %M can not be larger than the number of states
            M=min(obj.R,M);
            %Get Word Length.
            N=length(y);
                        
            %Compute APP Probabilities. -----------------------------------
            P_APP_n_1=zeros(1,N);
            for n=1:N  
                P_APP_n_1(n)=(tanh(A_ext_LLR_appended(n)/2)/2+0.5); 
            end            
            
            %Compute Alpha Probabilities ----------------------------------
            [alpha,alpha_survivors]=obj.compute_alphas(y,M,P_APP_n_1,N_0,pi_noise);
            %Calculate Betas ----------------------------------------------
            [beta]=obj.compute_betas(y,M,P_APP_n_1,N_0,pi_noise,alpha_survivors);

            
            %Calculate Final Probabilities --------------------------------
            lambda=alpha.*beta;
            LLL_a_n=zeros(1,N);

            for n=1:N-obj.M_T+1
                pos_set=lambda(obj.pos,n);
                neg_set=lambda(obj.neg,n);
                
                neg_sum=sum(neg_set);
                pos_sum=sum(pos_set);
                
                if(pos_sum==0 || neg_sum==0)
                     LLL_a_n(n)=10*sign(log((pos_sum/neg_sum)));
                else
                    LLL_a_n(n)=  log(pos_sum/neg_sum);
                end

            end
            % -------------------------------------------------------------

            %Subtract APP info.
            LLL_a_n=LLL_a_n-A_ext_LLR_appended;            
            
            %Clip
            max_LLR=38.1400;
            LLL_a_n(LLL_a_n> max_LLR) = max_LLR;
            LLL_a_n(LLL_a_n<-max_LLR) =-max_LLR;      
           
            
            %Calculate Soft Outputs.
            SO=LLL_a_n(1:end-obj.M_T+1);
            SO=SO(:);
            
            %Remove any Nans
            if(sum(isnan(SO))>0)
                SO(isnan(SO))=0;
            end

        end
        

        
        function [alpha_trim,alpha_survivors] = compute_alphas(obj,y,M,P_APP_n_1,N_0,pi_noise)

                %Get Word Length.
                N=length(y);            

                %Calculate Alphas -----------------------------------------
                alpha=zeros(obj.R,N+1);
                alpha(1,1)=1; %i.e., starts in the all -1 state.

                alpha_survivors=zeros(M,N);
                alpha_survivors(1,1)=1; %i.e., starts in the all -1 state. 

                J_0=1;
                J_1=2;                  


                for n=1:N  
                    for I=alpha_survivors(:,n)'    
                        % Only compute remaining if greater than 0.
                        if(not(I==0))
                            % Only compute if APP greater than 0.
                            if(P_APP_n_1(n)<1)
                                %Calculate Gamma prob for -1
                                gamma_sing_0=pi_noise(n)*exp(-(y(n)-obj.l(I,obj.S_alpha_set(I,J_0)))^2/N_0(n)); 
                                %Calulate Alpha prob
                                alpha(obj.S_alpha_set(I,J_0),n+1)=alpha(obj.S_alpha_set(I,J_0),n+1) + alpha(I,n)*gamma_sing_0*(1-P_APP_n_1(n)); 
                            end

                            %Only compute remaining if greater than 0.
                            if(P_APP_n_1(n)>0)
                                %Calculate Gamma prob
                                gamma_sing_1=pi_noise(n)*exp(-(y(n)-obj.l(I,obj.S_alpha_set(I,J_1)))^2/N_0(n)); 
                                %Calulate Alpha prob
                                alpha(obj.S_alpha_set(I,J_1),n+1)=alpha(obj.S_alpha_set(I,J_1),n+1) + alpha(I,n)*gamma_sing_1*P_APP_n_1(n); 
                            end  
                        end
                    end

                    %Normalize
                    alpha_sum=sum(alpha(:,n+1));
                    if(alpha_sum==0)
                        J_1_set=obj.S_alpha_set((alpha(:,n)>0),J_1);
                        J_0_set=obj.S_alpha_set((alpha(:,n)>0),J_0);

                        prob=1/(length(J_1_set)+length(J_0_set));
                        alpha(J_1_set,n+1)=prob;
                        alpha(J_0_set,n+1)=prob;
                    else
                        alpha(:,n+1)=alpha(:,n+1)/alpha_sum; 
                    end

                    % Compute M - Survivors
                    sorted_alpha=sort(alpha(:,n+1));
                    if(sorted_alpha(end)==1)
                        alpha_survivors(1,n+1)=find(alpha(:,n+1)==sorted_alpha(end));
                    else
                        val_thresh=sorted_alpha(end-M+1);
                        if(val_thresh==0)
                            temp=sorted_alpha;
                            temp(temp==0)=[];
                            val_thresh=temp(1);
                        end
                        surviv_indices=find((alpha(:,n+1)>=val_thresh));
                        alpha_survivors(1:length(surviv_indices),n+1)=surviv_indices;
                        all_indices=1:obj.R;
                        all_indices(alpha_survivors(1:length(surviv_indices),n+1))=[];
                        alpha(all_indices,n+1)=0;            
                    end

                    %Error Checking
        %                 %Check for errors.
        %                 if(alpha_sum==0)
        %                     error
        %                 end
        %                 if(sum(alpha(:,n+1))==0)
        %                     error
        %                 end                
        %                 if(length(alpha_survivors(:,n+1))>M)
        %                    error
        %                 end

                end
                
                alpha_trim=alpha(:,2:N-1);
        end

    
    

        function [beta_trim] = compute_betas(obj,y,M,P_APP_n_1,N_0,pi_noise,alpha_survivors)
            
            
            %Get Word Length.
            N=length(y);                

            beta=zeros(obj.R,N+1);
            beta(1,N+1)=1; %i.e., ends int he all -1 state
            P_APP_n_1=[P_APP_n_1,-38]; % update APP info.


            beta_survivors=zeros(M,1);
            beta_survivors(1,1)=1; %i.e., starts in the all -1 state. 

  
            for n=N:-1:1
                for J=beta_survivors(:)'
                    if(not(J==0))
                        
                        %Branch metrics for +1.
                        origin_pos=obj.S_beta_set(J,1);
                        comp_val= obj.l_beta_out(origin_pos,1);

                        %Determine Sign.
                        if( obj.pos(origin_pos) )
                            P_APP= P_APP_n_1(n);
                        else
                            P_APP = 1-P_APP_n_1(n);
                        end

                        %Calulate Beta prob
                        gamma_sing=pi_noise(n)*exp(-(y(n)-comp_val)^2/N_0(n)); 
                        beta(origin_pos,n)=  beta(origin_pos,n)+ gamma_sing*P_APP*beta(J,n+1); 

                         %Branch metrics for -1.
                        origin_pos=obj.S_beta_set(J,2);
                        comp_val= obj.l_beta_out(origin_pos,2);

                        %Determine Sign.
                        if( obj.pos(origin_pos) )
                            P_APP= P_APP_n_1(n);
                        else
                            P_APP = 1-P_APP_n_1(n);
                        end

                        %Calulate Beta prob
                        gamma_sing=pi_noise(n)*exp(-(y(n)-comp_val)^2/N_0(n)); 
                        beta(origin_pos,n)= beta(origin_pos,n)+ gamma_sing*P_APP*beta(J,n+1);        
             
                    end
                end

                 %Normalize
                 beta_sum=sum(beta(:,n));
                 if(beta_sum==0)
                    num_alphas=sum(alpha_survivors(:,n+1)>0);
                    alpha_survivors_list=alpha_survivors(1:num_alphas,n+1);
                    alpha_survivors_list(alpha_survivors_list==0)=[];
                    beta(alpha_survivors_list,n)=1/length(alpha_survivors_list);
                 else
                     beta(:,n)=beta(:,n)./beta_sum;
                 end

                % Use alpha survivor list.
                num_alphas=sum(alpha_survivors(:,n+1)>0);
                alpha_survivors_list=alpha_survivors(1:num_alphas,n+1);

                % Compute M - Survivors
                sorted_beta=sort(beta(:,n));
                if(sorted_beta(end)==1)
                    loser_indices=(beta(:,n)<sorted_beta(end));
                    beta_survivors=find(beta(:,n)==sorted_beta(end));
                    beta(loser_indices,n)=0;   
                else
                    if(num_alphas>=M)
                        all_indices=1:obj.R;
                        all_indices(alpha_survivors_list)=[];               
                        beta(all_indices,n)=0;
                        beta_survivors=alpha_survivors_list;
                    else
                        m_ex=M-num_alphas;

                        beta_temp=beta(:,n);                 
                        beta_temp(alpha_survivors_list)=0;

                        sorted_beta=sort(beta_temp);
                        beta_surv=find( (beta(:,n)>sorted_beta(end-m_ex)) );
                        if(length(beta_surv)>m_ex)
                            beta_surv=beta_surv(1:m_ex);
                        end

                        beta_survivors=[alpha_survivors_list ; beta_surv ] ;

                        all_indices=1:obj.R;
                        all_indices(beta_survivors)=[];               
                        beta(all_indices,n)=0;                
                    end
                end

                %Error Checking
%                if(beta_sum==0)
%                       error
%                end
%                if(sum(beta(:,n))==0)
%                        error
%                end                
%                if(length(beta_survivors)>M)
%                        error
%                end
            end
            beta_trim=beta(:,1:N-2);
        end
    
    end
    
    
    
  
    
    
end



