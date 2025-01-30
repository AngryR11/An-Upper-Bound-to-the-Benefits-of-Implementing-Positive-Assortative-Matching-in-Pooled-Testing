aa=.2#alpha parameter for the equity measure
equity_f<-function(x){((x)^(1-aa))/(1-aa)}#equity function

#library("demography")#used for the spline







Ordered_ET_FN_FP_EQT<-function(q,k,lambda,dilution){
  Se=dilution(1,1)
  Sp=1-dilution(0,1)
  utilities=rep(0,length(q))#expected utility from each subject
  PFN_I_M=matrix(0, nrow = 1, ncol = (k+1))#probability the test accuses false negative of an infected patient when there are I infected in the group
  PFP_I_M=matrix(0, nrow = 1, ncol = (k+1))
  for(i in 1:k)
  {
    PFP_I_M[1,(i)]=(1-Sp)*dilution((i-1),k)	
    PFN_I_M[1,(i+1)]=1-Se*dilution(i,k)
  }	
  
  E_FN_i=0#expected numb of false negatives
  E_tests_i=0#expected numb of tests
  E_FP_i=0#expected numb of false positives
  E_equity=0#expected equitability
  q=sort(q)#ording patients on prob of infection
  # prob_flip=true_proportion_in_population[2]* prob_misreport/(true_proportion_in_population[2]* prob_misreport+true_proportion_in_population[1])
  # high_q=max(prevalence_per_group)
  # low_q=min(prevalence_per_group)
  # for(l in 1:length(q)){
  # if(q[l]<high_q){
  # do_i_flip=rbinom(n=1,size=1,prob=prob_flip)
  # q[l]=do_i_flip*high_q+(1-do_i_flip)*q[l]
  # }
  # } 
  if(k==1){#if everyone is tested individually
    E_tests_i=length(q)
    E_FN_i=sum(q)*(1-Se)
    E_FP_i=sum(1-q)*(1-Sp)
    for(l in 1:length(q)){
      E_equity=(equity_f(1-lambda*q[l]*(1-Se)-(1-lambda)*(1-q[l])*(1-Sp)))+E_equity
      utilities[l]=-lambda*q[l]*(1-Se)-(1-lambda)*(1-q[l])*(1-Sp)
    }
  }else{						
    initial_l=1
    Numb_Groups=length(q)%/%k
    for(l in 1:Numb_Groups){
      for(b in 1:k){
        coordinate=initial_l+b-1
        qs=q[initial_l:(k+initial_l-1)]
        qi=qs[b]#prob infection of current subject
        qs[b]=0
        q_i=qs[!qs==0]
        FN_m=0
        FP_m=0
        for(numb_inf in 0:(k-1)){
          FN_m=PFN_I_M[(numb_inf+2)]*qi*dpoisbinom(numb_inf, q_i, log_d = FALSE)+FN_m	
          FP_m=PFP_I_M[(numb_inf+1)]*(1-qi)*dpoisbinom(numb_inf, q_i, log_d=FALSE)+FP_m	
        }	
        E_equity=equity_f(1-lambda*FN_m-(1-lambda)*FP_m)+E_equity
        utilities[coordinate]=-lambda*FN_m-(1-lambda)*FP_m
      }		
      
      
      for(numb_inf in 0:k){
        E_FN_i = numb_inf*(1-Se*dilution(numb_inf,k))*dpoisbinom(numb_inf,q[initial_l:(k+initial_l-1)], log_d = FALSE)+E_FN_i
        
        E_FP_i = (k-numb_inf)*(1-Sp)*dilution(numb_inf,k)*dpoisbinom(numb_inf,q[initial_l:(k+initial_l-1)], log_d = FALSE)+E_FP_i
        
        E_tests_i=(k*dilution(numb_inf,k))*dpoisbinom(numb_inf,q[initial_l:(k+initial_l-1)], log_d = FALSE)+E_tests_i
        
        #E_equity_i= PFN_I((numb_inf),k)*(numb_inf/k)*dpoisbinom(numb_inf,q[initial_l:(k+initial_l-1)], log_d = FALSE)+E_equity_i	
      }
      E_tests_i=E_tests_i+1
      #E_equity=(k*equity_f(1-E_equity_i))+E_equity
      initial_l=initial_l+k
    }
    if(length(q)%%k==1){#if the remainder between n/k is 1, the agent with highest prob. of infection is tested alone.
      E_tests_i= E_tests_i+1	
      E_FP_i=E_FP_i+(1-Sp)*(1-q[length(q)])
      E_FN_i=E_FN_i+(1-Se)*q[length(q)]
      utilities[length(q)]=-lambda*q[length(q)]*(1-Se)-(1-lambda)*(1-q[length(q)])*(1-Sp)	
      E_equity=(equity_f(1-lambda*q[length(q)]*(1-Se)-(1-lambda)*(1-q[length(q)])*(1-Sp)))+E_equity
    }
    
    if(length(q)%%k>1){
      #if the remainder between n/k is greater than 1, the agents with highest prob. of infection are tested together.
      rem=length(q)%%k	
      #E_equity_i=0
      for(b in 1:rem){
        coordinate=initial_l+b-1
        qs=q[initial_l:(initial_l+rem-1)]
        qi=qs[b]#prob infection of current subject
        qs[b]=0
        q_i=qs[!qs==0]
        FN_m=0
        FP_m=0
        for(numb_inf in 0:(rem-1)){
          FN_m= PFN_I_M[(numb_inf+2)]*qi*dpoisbinom(numb_inf, q_i, log_d = FALSE)+FN_m	
          FP_m=PFP_I_M[(numb_inf+1)]*(1-qi)*dpoisbinom(numb_inf, q_i, log_d=FALSE)+FP_m	
        }	
        utilities[coordinate]=-lambda*FN_m-(1-lambda)*FP_m
        E_equity=equity_f(1-lambda*FN_m-(1-lambda)*FP_m)+E_equity
      }	
      
      
      
      for(numb_inf in 0:rem){
        E_tests_i=(rem*dilution(numb_inf,rem))*dpoisbinom(numb_inf,q[(1+length(q)-rem):length(q)], log_d = FALSE)+E_tests_i	
        E_FP_i=(rem-numb_inf)*(1-Sp)*dilution(numb_inf, rem)*dpoisbinom(numb_inf,q[(1+length(q)-rem):length(q)], log_d = FALSE)+E_FP_i
        E_FN_i = numb_inf*(1-Se*dilution(numb_inf,rem))*dpoisbinom(numb_inf,q[(1+length(q)-rem):length(q)], log_d = FALSE)+E_FN_i	
        #E_equity_i= PFN_I((numb_inf),rem)*(numb_inf/rem)*dpoisbinom(numb_inf,q[(1+length(q)-rem):length(q)], log_d = FALSE)+E_equity_i
      }
      E_tests_i=E_tests_i+1
      #E_equity=(rem*equity_f(1-E_equity_i))+E_equity
    }
  }
  if(length(q)%%k==0){#if all groups have the same size
    if(lambda==0){
      min_util_last_group=utilities[1]#min(utilities[1:k])
    }
    else{
      min_util_last_group=utilities[length(q)]#min(utilities[(length(q)-k+1):length(q)])	
    }
  }
  else{
    min_util_last_group=0#NA
  }
  return(c(E_tests_i/length(q),E_FP_i/length(q),E_FN_i/length(q)))#c(E_tests_i,E_FN_i,E_FP_i,E_equity,min(utilities), min_util_last_group))
}



#==========Random Matching=========
#To compute the expected numb of tests, false neg., etc., from random pooling, one only needs information on prevalence: p=probability a generic agent is infected, the pool sizes k, and the number of people tested, n.


prob_y_random<-function(x,p,k){(factorial(k)/(factorial(k-x)*factorial(x)))*(p^x)*((1-p)^(k-x))}
prob_y_random_i<-function(x,p,k){(factorial(k-1)/(factorial(k-1-x)*factorial(x)))*(p^x)*((1-p)^(k-1-x))}

Random_ET_FN_FP_EQT<-function(p,k,n,dilution){
  Se=dilution(1,1)
  Sp=1-dilution(0,1)
  Prob_infected_random=0
  E_false_negatives_random=0
  E_false_positives_random=0
  FN=0
  Numb_Groups=n%/%k
  if(k==1){#if everyone is tested individually
    E_test_random=n
    E_false_negatives_random=n*p*(1-Se)
    E_false_positives_random =n*(1-p)*(1-Sp)
    E_equity_random =n*(equity_f(1-p*(1-Se)))
  }
  else{	
    for(y in 0:k){
      if(y<k){
        FN=FN+prob_y_random_i(y,p,k)*(1-Se*dilution(y+1,k))*p
      }
      Prob_infected_random=Prob_infected_random+prob_y_random(y,p,k)*dilution(y,k)
      E_false_negatives_random=E_false_negatives_random+ Numb_Groups*y*prob_y_random(y,p,k)*(1-Se*dilution(y,k))
      E_false_positives_random=E_false_positives_random+Numb_Groups*(k-y)*prob_y_random(y,p,k)*dilution(y,k)*(1-Sp)
    }	
    E_equity_random=n*equity_f(1-FN)
    E_test_random=n%/%k+k*n%/%k*Prob_infected_random
    
    if(n%%k==1){#if the remainder between n/k is 1, the agent with highest prob. of infection is tested alone.
      E_test_random = E_test_random +1	
      E_false_positives_random = E_false_positives_random +(1-Sp)*(1-p)
      E_false_negatives_random = E_false_negatives_random +(1-Se)*p
      E_equity_random = E_equity_random+(equity_f(1-p*(1-Se)))
    }
    
    if(n%%k>1){#if the remainder between n/k is 1, the agent with highest prob. of infection is tested alone.
      rem=n%%k	
      Prob_infected_random=0
      FN=0
      for(y in 0: rem){
        if(y<rem){
          FN=FN +prob_y_random_i(y,p,rem)*(1-Se*dilution(y+1,rem))*p
        }
        Prob_infected_random=Prob_infected_random+prob_y_random(y,p,rem)*dilution(y,rem)
        E_false_negatives_random=E_false_negatives_random+y*prob_y_random(y,p,rem)*(1-Se*dilution(y,rem))
        E_false_positives_random=E_false_positives_random+(rem-y)*prob_y_random(y,p,rem)*dilution(y,rem)*(1-Sp)
      }	
      E_equity_random=E_equity_random+rem*equity_f(1-FN)
      E_test_random=E_test_random+1+rem*Prob_infected_random
    }
  }	
  return(c(E_test_random/n, E_false_positives_random/n,E_false_negatives_random/n))#c(E_test_random, E_false_negatives_random, E_false_positives_random, E_equity_random))
}



max_var_q_2=function(a,b,n,p,k){
	n=floor(n/k)*k#to force n to be a multiple of K
    q<-rep(b,n)
       i=1
		while(mean(q)>p){
			q[i:(i+k-1)]=rep(a,k)
			i=i+k
		}
		i=i-k#adjusting the probability of infection of the intermediate group
		if(i<n-k & i>1){
		q[i:(i+k-1)]=rep((n*p-sum(q[1:(i-1)])-sum(q[(i+k):n]))/k,k)
		}else if(i==1){
		  q[i:(i+k-1)]=rep((n*p-sum(q[(i+k):n]))/k,k)
		}
		else{q[(n-k+1):n]=rep((n*p-sum(q[1:(i-1)]))/k,k)}
		#if there is only one pool, everyone must have probability of infection equal to p
		if(n==k){
			q<-rep(p,n)
		}
	return(q)
}


#==ET==
asymptotic_UB_ET=function(a,b,qbar,k,dilution){
  alpha=(b-qbar)/(b-a)
  ET=0
  for(I in 0:k){
    ET= ET+dilution(I,k)*(dbinom(I, size=k, prob=qbar)-alpha*dbinom(I, size=k, prob=a)-(1-alpha)*dbinom(I, size=k, prob=b))
  }	
  return(ET)
}


#==FP==
asymptotic_UB_FP=function(a,b,qbar,k,dilution){
  alpha=(b-qbar)/(b-a)
  ET=0
  for(I in 0:(k-1)){
    ET= ET+dilution(0,1)*dilution(I,k)*(dbinom(I, size=(k-1), prob=qbar)*(1-qbar)-alpha*dbinom(I, size=(k-1), prob=a)*(1-a)-(1-alpha)*dbinom(I, size=(k-1), prob=b)*(1-b))
  }	
  return(ET)
}
#==FN==
asymptotic_UB_FN=function(a,b,qbar,k,dilution){
  alpha=(b-qbar)/(b-a)
  ET=0
  for(I in 0:(k-1)){
    ET= ET+(1-dilution(1,1)*dilution(I+1,k))*(dbinom(I, size=(k-1), prob=qbar)*qbar-alpha*dbinom(I, size=(k-1), prob=a)*a-(1-alpha)*dbinom(I, size=(k-1), prob=b)*b)
  }	
  return(ET)
}



#Delete later for the love of God:
# qa=0.01
# qb=.3
# qbar=0.05
# k=10
# n_grid=seq(1,200)
# upper_bound_diff_ET_k=n_grid
# random_ET_k=n_grid
# asymptotic_UB=asymptotic_UB_ET(qa,qb,qbar,k,dilution)
# for(i in 1:length(n_grid)){
#   print(n_grid[i])
#   n=n_grid[i]%*%k
#   extreme_qs=max_var_q_2(qa,qb,n,qbar,k)
#   
#   #random pooling:
#   Variables=Random_ET_FN_FP_EQT(qbar,k,n,dilution)
#   ET= Variables[1]
#   random_ET_k[i]=ET
#   upper_bound_diff_ET_k[i]=ET-Ordered_ET_FN_FP_EQT(extreme_qs,k,lambda,dilution)[1]
# }
# 
# 
# 
# #pdf(file="/Users/gustavosaraiva/Documents/Chile_Ano_2/Coronavirus_Pool_Testing/Latex/asymptotic_UB.pdf", width=7, height=5)  # Adjust 'width' and 'height' as needed
# old_par <- par()  # Save current settings to restore later
# par(mar=c(5, 5, 4, 2) + 0.1)  # Default is c(5, 4, 4, 2) + 0.1, increase the left margin
# plot(n_grid, upper_bound_diff_ET_k, type='l',
#      col="black", # Changes the line color to blue
#      lwd=2, # Increases the line width
#      xlab="n", # X-axis label
#      ylab=expression(UB[o]^T * "(" * mu * "," * n * "," * K * ")"), # Y-axis label
#      cex.lab=1.2, # Increases the size of the axis labels
#      cex.main=1.4, # Increases the size of the main title
#      cex.axis=1.1, # Increases the size of the axis annotation (numbers)
#      ylim=c(min(upper_bound_diff_ET_k), max(upper_bound_diff_ET_k) * 1.1) # Adjusts Y-axis limits for better visibility
# )
# # Add a horizontal line at the specified 'asymptotic_UB' point
# abline(h=asymptotic_UB, col="blue", lwd=2, lty=2)  # 'lty=2' makes the line dashed
# # Label the horizontal line
# # Find a suitable x-position for the text: e.g., the right-most x-value or a specific value that fits well within the plot
# x_pos_label <- max(n_grid)-45
# y_pos_label <- asymptotic_UB+.0039
# text(x=x_pos_label, y=y_pos_label, labels=expression(lim(UB[o]^T* "(" * mu * "," * n * "," * K * ")",n %->% infinity)), pos=4, col="blue", cex=0.8)
# 
# # Restore previous par settings after plotting
# par(old_par)
# #dev.off()
# 
# 
# 
# 
