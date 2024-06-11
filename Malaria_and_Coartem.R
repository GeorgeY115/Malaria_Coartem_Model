rm=list(ls())

install.packages('deSolve')
install.packages('ggplot')
library(ggplot)
library(deSolve)

Treatment_Resistance_Model_w_ITNs <- function(a,b,m,w,e,s,o,l,d,z,t,i){
  require(deSolve)
  init <- c(S=1-1e-6, I=1e-6, Z=0, R=0, X=1-1e-6, Y=1e-6)
  parameters <- c(alpha=a,beta=b,mu=m,delta=d,lambda=l,epsilon=e,sigma=s,omega=o,tau=w,zeta=z,omicron=i)
  time <- seq(0,t,by=t/(2*length(1:t)))
  eqn <- function(time, state, parameters) {
  
    S <- state[1]
    I <- state[2]
    Z <- state[3]
    R <- state[4]
    X <- state[5]
    Y <- state[6]

    
    with(as.list(c(parameters)), { 
      dS <- mu*(S+I+Z+R) - ((omicron*beta*S*Y)/(X+Y)) - mu*S + tau*R
      dI <- ((omicron*beta*S*Y)/(X+Y)) + epsilon*Z - lambda*I - mu*I - sigma*I
      dZ <- sigma*I - epsilon*Z - mu*Z - alpha*Z
      dR <- lambda*I - mu*R - tau*R + alpha*Z
      dX <- delta*(X+Y) - ((omega*X*I)/(S+I+Z+R)) - zeta*X
      dY <- ((omega*X*I)/(S+I+Z+R)) - zeta*Y

      
      return(list(c(dS, dI, dZ, dR, dX, dY)))
    })
  }
  out<-ode(y=init,times=time,eqn,parms=parameters) 
  out.df<-as.data.frame(out) #create a data frame 
  print(tail(out.df))
  require(ggplot2) 
  title <- bquote("Coartem as Malaria Treatment in Kenya") 
  subtit<-bquote(list(alpha==.(a),beta==.(b),delta==.(d),lambda==.(l),tau==.(w),zeta==.(z),epsilon==.(e),omega==.(o),mu==.(m),sigma==.(s),omicron==.(i)))
  res<-ggplot(out.df,aes(x=time))+ 
    ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit))))))+ 
    geom_line(aes(y=S,colour="Susceptible"))+
    geom_line(aes(y=I,colour="Infected"))+ 
    geom_line(aes(y=Z,colour="Treated"))+
    geom_line(aes(y=R,colour="Recovered"))+
    geom_line(aes(y=X,colour="Susceptible Mosquitoes"))+
    geom_line(aes(y=Y,colour="Infected Mosquitoes"))+ # hash out the X and Y lines if you only want the human compartments on the graph.

    ylim(0,1)+
    xlim(0,t)+
    ylab(label="Proportion")+ #y-axis label
    xlab(label="Time (days)")+ #x-axis label 
    labs(colour='Group') +
    theme_bw()

print(res) #print  plot
}

Treatment_Resistance_Model_w_ITNs(t=250, a=0.225, b=0.6, m=0.0005, l=0.009, w=0.002, e=0.009, s=0.128, d=0.05, o=0.6, z=0.05, i=10)

#PARAMETER SOURCES:
  # a (treatment success) = 0.225 -> effectiveness for coartem estimated to be around 99%,
  # over a period of 3-6 days (lumefantrine half life).
  # https://pubmed.ncbi.nlm.nih.gov/33397548/#:~:text=In%20Kenya%2C%20artemether%2Dlumefantrine%20(,line%20treatment%20of%20uncomplicated%20malaria.

  # b (effective bite rate) = 0.6 -> the average infected Anopheles mosquito has ~500 sporozoites,
  # a level which is correlated with a roughly 60% chance of transmission per bite.
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5230737/

  # m (human birth/death rate) = 0.0005 -> in Kenya, both rates are ~2 per 100 people per year. (0.027/365).
  # https://www.macrotrends.net/global-metrics/countries/KEN/kenya/birth-rate

  # l (recovery without treatment) = 0.009 -> 50% recovery rate without treatment, 8 week av infection.
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2607295/

  # w (waning immunity) = 0.002 -> immunity shown to wane after 14 months.
  # https://pubmed.ncbi.nlm.nih.gov/22144890/

  # e (treatment failure) = 0.009 -> baseline drug failure = 1%, along with recent Coartem resistance observed at
  # average of 3% in Kenya, taken from a weighted average of several studies. Divided by 4.5 (av days until treatment).
  # Using data from https://www.macrotrends.net/global-metrics/countries/KEN/kenya/birth-rate

  # s (proportion of cases treated) = 0.128 -> Of 81.4% uncomplicated cases, 51.8% were treated. Of 18.6% which were severe, 83.9% were treated.
  # Weighted, this means ~57.8% of cases in Kenya recieve treatment. Divided by 4.5.

  # d,z (birth and death rate of female mosquitoes) = 0.5 -> arbitrary, assumed to be equal.

  # i (scaling factor) = 10 -> arbitrary, assumed 



