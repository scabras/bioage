## ----setup, include=FALSE---------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = FALSE,warning = FALSE,cache=TRUE,
                      message=FALSE,comment=FALSE)
library(DiagrammeR)
#library(DiagrammeRsvg)
library(rsvg)

library(latex2exp)
library(kableExtra)
library(keras)
library(tensorflow)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(LaplacesDemon)
library(Anthropometry)
library(fastDummies)


## ----fig.cap="BA point estimation for the three subjects from different models in the toy example along with the observed values (black bullets) and the true BA values \\label{fig-toyexlines}"----
# Define the data
theta <- c(5, -10, 0)
cage <- c(20, 30, 40)  # Chronological Age
phenotype <- cage + theta # Phenotype values
zcage=(cage-0)/sqrt(mean((cage-mean(cage))^2))
olsmod=lm(phenotype~cage)

data <- data.frame(
  phenotype = phenotype,
  cage = cage,
  true.bioage = cage + theta,
  beta1m0=1,
  beta1m1=2,
  beta1m2=0.5,
  beta1mols=coef(olsmod)[2],
  beta0mols=coef(olsmod)[1],
  biom1 = phenotype / 2, 
  biom2 = phenotype / 0.5) %>%
  mutate(biomols= (phenotype-beta0mols)/beta1mols,
         sigmam0=sqrt(mean((phenotype-beta1m0*cage)^2)),
         sigmam1=sqrt(mean((phenotype-beta1m1*cage)^2)),
         sigmam2=sqrt(mean((phenotype-beta1m2*cage)^2)),
         sigmay=sqrt(mean((phenotype-mean(phenotype))^2)),
         sigmamols=sqrt(mean((phenotype-beta0mols-beta1mols*cage)^2)))

datainv=data
datainv$theta=theta

modnames=c("Correct (True)", "Underestimation", "Overestimation","OLS")
# Create a long format data frame for correct mapping
plot_data <- data.frame(
  phenotype = rep(data$phenotype, 4),
  age = c(data$true.bioage, data$biom1, data$biom2,data$biomols),
  category = rep(modnames, each = nrow(data))
)

# Define colors for legend
colors <- c("Correct (True)" = "grey", "Underestimation" = "red", 
            "Overestimation" = "orange","OLS"="black")

# Plot
p1 <- ggplot() +
  # Add main scatter points
  geom_point(data = data, aes(x = cage, y = phenotype), size = 2) + xlim(c(0,80))+
  
  # Add legend-associated points
  geom_point(data = plot_data, aes(x = age, y = phenotype, 
                                   color = category, shape = category), 
             size = 4) +
  
  # Add identity reference lines with legend colors
  geom_abline(aes(intercept = 0, slope = 1, color = "Correct (True)"), size = 1) +
  geom_abline(aes(intercept = 0, slope = 2, color = "Underestimation"), size = 1) +
  geom_abline(aes(intercept = 0, slope = 0.5, color = "Overestimation"), size = 1) +
  geom_abline(aes(intercept = data$beta0mols, slope = data$beta1mols, color = "OLS"), size = 2) +
  
  # Set colors and shapes for legend
  scale_color_manual(values = colors, name = "Conditional Estimations") +
  scale_shape_manual(values = c(2, 3, 3,3), name = "Conditional Estimations") +

  # Labels and theme
  xlab("Chronological and Biological Age") +
  ylab("Phenotype") +
  ggtitle("Observed Phenotype vs Chronological Age") +
  theme_minimal()

# Print the plot
print(p1)


## ----fig.cap="Normalized profile likelihoods for $\\boldsymbol \\theta$ from different models in the toy example along with the true BA values  \\label{fig-toyexlik}"----
# Define the likelihood function
llprofile <- function(theta, y,x, beta0=0, beta1,sigma=1) 
  dnorm(y, mean = beta0 + beta1 * (x + theta), sd = sigma, log = TRUE)

# Calculate the second derivative of the log likelihood
llprofile2 <- function(theta, y,x, beta0=0, beta1,sigma) -beta1^2 / sigma^2



# Lets draw for each one of the three subject the corresponding 
# log likelihood profile for each of the three models
ss=seq(-50,50,length=10000)
llres=data.frame(NULL)
for(i in 1:3) for(j in 1:4){
  nll=llprofile(theta=ss,y = data$phenotype[i],x = data$cage[i],
                beta1 = data[i,grep("beta1",colnames(data))[j]],
                beta0 = data$beta0mols[1]*(i==4),
                sigma = data[i,grep("sigmam",colnames(data))[j]])
  nll=nll-max(nll)
  llres=rbind(llres,data.frame(theta=ss,loglik=nll,model=paste0("Model",j-1),
                               subject=paste0("Subject",i)))
}  
                                        
llres %>% 
  mutate(model=factor(model,labels=modnames),subject=factor(subject)) -> llres

# Create a separate dataset for vlines with unique bioage per subject
vline_data <- data.frame(
  subject = factor(paste0("Subject", 1:3)), # Ensure factors match facet levels
  xintercept = theta
)

# Plot with vertical lines appearing only once per subject
ggplot(data = llres, aes(x = theta, y = exp(loglik), color = factor(model))) +
  geom_line() +
  facet_grid(subject ~ ., scales = "free_y") +  # Allow different y-scales per subject
  theme_minimal() +

  # Add vlines using separate dataset
  geom_vline(data = vline_data, aes(xintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +

  # Define color legend
  scale_color_manual(values = colors, name = "Conditional Estimation") +

  # Labels and aesthetics
  xlab(TeX("$\\theta$")) + 
  ylab("Normalized profile likelihood") + 
  ggtitle("Profile likelihoods for the three Subjects") +
  
  # Remove duplicate y-axis labels and improve facet formatting
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1), # Align facet labels
    axis.title.y = element_blank(),                   # Remove global y-axis title
    strip.placement = "outside",                      # Move facet labels outside
    panel.spacing = unit(1, "lines"),                 # Space panels properly
    axis.text.y = element_blank(),                    # Remove individual y-axis labels
    axis.ticks.y = element_blank()                    # Remove y-axis ticks
  )


## ---------------------------------------------------------------
infoexact=data %>% select(starts_with("beta1m"),starts_with("sigmam"))
infoexact=infoexact[1,1:4]^2/infoexact[1,5:8]^2
infonum=NULL
for(i in 1:4) infonum=c(infonum,
                        -numDeriv::hessian(llprofile,
                                           theta=(data[c("true.bioage","biom1",
                                                         "biom2","biomols")[i]]-
                                                    data$cage)[1,],
               y = data$phenotype[1],x = data$cage[1],
               beta1 = data[c("beta1m0", "beta1m1", "beta1m2","beta1mols")[i]][1,],
               sigma = data[c("sigmam0", "sigmam1", "sigmam2",
                              "sigmamols")[i]][1,],
               beta0=data$beta0mols[1]*(i==4)))

sigmacage=sd(data$cage)*sqrt(2/3)

R2=1-(data$sigmamols[1]^2/data$sigmay[1]^2)
R2=(data$beta1mols[1]^2)/data$sigmay[1]^2*sigmacage^2
# Signal to noise ratio increase the uncertainty using different age people
infoR2=(R2/(1-R2))/sigmacage^2

names(infonum)=c("Overestimation","Underestimation","Correct","OLS")
infonum=round(infonum,3)


## ----results='asis'---------------------------------------------
kable(data.frame(t(infonum)), format = "latex", escape = FALSE,
      caption = "\\label{tab:obsfi} Observed Fisher Information of $\\theta_i$ for the considered models",
      row.names = 0,align = "c",booktabs = TRUE) %>%
  kable_styling(latex_options = c("HOLD_position"))


## ----fig.cap="Full conditional posterior for $\\theta_1$, $\\theta_2$ and $\\theta_3$ in the toy example along with the true values.  \\label{fig-toyexfctheta}"----

# Define the prior st for theta such that with high probability the 
# BA is between -15 and 15 years from the chronological age.
sigma.theta=7.5

# Define the exact full conditional posterior
lpostfullcond <- function(theta, y,x, beta0, beta1,sigmaeps){
  mutheta=(beta1/sigmaeps^2)*(y-beta0-beta1*x)/(1/sigma.theta^2+beta1^2/sigmaeps^2)
  sigmatheta=sqrt(1/(1/sigma.theta^2+beta1^2/sigmaeps^2))
  dnorm(theta,mean=mutheta,sd=sigmatheta,log=TRUE)}

rpostfullcondtheta <- function(y,x, beta0, beta1,sigmaeps,sigma.theta=7.5,n){
  mutheta=(beta1/sigmaeps^2)*(y-beta0-beta1*x)/(1/sigma.theta^2+beta1^2/sigmaeps^2)
  sigmatheta=sqrt(1/(1/sigma.theta^2+beta1^2/sigmaeps^2))
  rnorm(n,mean=mutheta,sd=sigmatheta)}

# lets draw the conditional posterior for the three subjects
ss=seq(-50,50,length=10000)
postres=data.frame(NULL)
for(i in 1:3) for(j in 1:3){
  yi=data$phenotype[i]
  ci=data$cage[i]
  beta1=data[i,grep("beta1",colnames(data))[j]]
  sigma=data[i,grep("sigmam",colnames(data))[j]]
  exactpost=exp(lpostfullcond(theta=ss,y = yi,x =ci,
                              beta1 = beta1,beta0 = 0,sigma = sigma))
  postres=rbind(postres,data.frame(theta=ss,
                                   cpost=exactpost,model=paste0("Model",j-1),
                                   subject=paste0("Subject",i)))
}

postres %>% mutate(model=factor(model,labels=modnames[1:3]),
                   subject=factor(subject)) -> llres 

# Create a separate dataset for vlines with unique bioage per subject
vline_data <- data.frame(
  subject = factor(paste0("Subject", 1:3)), # Ensure factors match facet levels
  xintercept = theta
)

# Plot with vertical lines appearing only once per subject
ggplot(data = llres, aes(x = theta, y = cpost, color = factor(model))) +
  geom_line() +
  facet_grid(subject ~ ., scales = "free_y") +  # Allow different y-scales per subject
  theme_minimal() +

  # Add vlines using separate dataset
  geom_vline(data = vline_data, aes(xintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +

  # Define color legend
  scale_color_manual(values = colors, name = "Full conditional posterior") +
  # Labels and aesthetics
  xlab(TeX("$\\theta$")) + 
  ylab("Posterior Density") + 
  ggtitle("Conditional full posterior for the subjects") +
  
  # Remove duplicate y-axis labels and improve facet formatting
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1), # Align facet labels
    axis.title.y = element_blank(),                   # Remove global y-axis title
    strip.placement = "outside",                      # Move facet labels outside
    panel.spacing = unit(1, "lines"),                 # Space panels properly
    axis.text.y = element_blank(),                    # Remove individual y-axis labels
    axis.ticks.y = element_blank()                    # Remove y-axis ticks
  )


## ---------------------------------------------------------------
# Define the full conditional posterior for beta0 and beta1
rpostfullcondbeta0 <- function(beta1, y,x, sigmaeps,theta,n){
  mu0=(1/(n*sigmaeps^2))*sum(y-beta1*(x+theta))
  sig0=sigmaeps^2/n
  rnorm(1,mean=mu0,sd=sqrt(sig0))}

rpostfullcondbeta1 <- function(beta0, y,x, sigmaeps,theta,n){
  mu1=sum((y-beta0)*(x+theta))/sum((x+theta)^2)
  sig1=sigmaeps^2/sum((x+theta)^2)
  rnorm(1,mean=mu1,sd=sqrt(sig1))}

# Define the full conditional posterior for sigma.epsilon^2
rpostfullcondsigmaeps2 <- function(y,x, beta0,beta1,theta,sigmaeps.a,
                                   sigmaeps.b,n){
  n=length(y)
  an=(n/2+sigmaeps.a)
  bn=0.5*(sum((y-beta0-beta1*(x+theta))^2)+sigmaeps.b)
  rinvgamma(1,shape = an,scale = bn)
  }


## ---------------------------------------------------------------
gibbsapproxuni=function(R=10000,burnin=1000,data,sigmaeps.a=0.0001,
                        sigmaeps.b=0.0001,thinning=10){
  n=nrow(data)
  y=data$phenotype
  x=data$cage
  # Initialize the chain
  theta.init=rep(0,n)
  beta0.init=0
  beta1.init=1
  sigmaeps2.init=1
  chain=matrix(NA,nrow=R+burnin,ncol=n+3)
  colnames(chain)=c(paste0("theta",1:n),c("beta0","beta1","sigmaeps2"))
  chain[1,]=c(theta.init,beta0.init,beta1.init,sigmaeps2.init)
  for(r in 2:(R+burnin)){
    # Sample of theta
    chain[r,1:n]=rpostfullcondtheta(y=y,x=x,
                                  beta0=chain[r-1,n+1],beta1=chain[r-1,n+2],
                                  sigmaeps=sqrt(chain[r-1,n+3]),n=n)
    # Sample beta0 and beta1
    for(j in 1:3){
      chain[r,n+1]=rpostfullcondbeta0(beta1=chain[r-1,n+2],y=y,x=x,
                                        sigmaeps=sqrt(chain[r-1,n+3]),
                                      theta=chain[r,1:n],n=n)
      chain[r,n+2]=rpostfullcondbeta1(beta0=chain[r,n+1],y=y,x=x,
                                        sigmaeps=sqrt(chain[r-1,n+3]),
                                      theta=chain[r,1:n],n=n)
    }
    # Sample sigmaeps
    chain[r,n+3]=rpostfullcondsigmaeps2(y=y,x=x,beta0=chain[r,n+1],
                                    beta1=chain[r,n+2],theta=chain[r,1:n],
                                    sigmaeps.a=sigmaeps.a,
                                    sigmaeps.b=sigmaeps.b,n=n)
  }
  idx_keep <- seq(burnin + 1, R + burnin, by = thinning)
  return(chain[idx_keep,])
}



## ---------------------------------------------------------------
if (file.exists("res.gibbsapproxuni.RData")) {
  load("res.gibbsapproxuni.RData")
} else {
  res.gibbsapproxuni=data.frame(gibbsapproxuni(data=data))
  save(res.gibbsapproxuni, file = "res.gibbsapproxuni.RData")
}


## ----fig.cap="Posterior approximation of all parameters in the toy example (one phenotype). \\label{fig-toyexpostuni}"----
rr=res.gibbsapproxuni %>% mutate(sigmaeps=sqrt(sigmaeps2)) %>% 
  select(-sigmaeps2)
trueval=c(theta,beta0=0,beta1=1,sigmaeps=1)
names(trueval)=colnames(rr)

rr=rr %>% mutate(time=1:nrow(res.gibbsapproxuni)) %>% 
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")


rr$parameter= factor(rr$parameter, 
                             levels=c("theta1", "theta2", "theta3",
                                      "beta0", "beta1", "sigmaeps"),
                      labels = c("theta[1]", "theta[2]", "theta[3]", 
                                 "beta[0]", "beta[1]", "sigma[epsilon]"))


# Create a separate dataset for vlines with true values of parameters
vline_data <- data.frame(
  parameter = names(trueval), # Ensure factors match facet levels
  xintercept = trueval
)

vline_data$parameter = factor(vline_data$parameter,
                              levels=c("theta1", "theta2", "theta3",
                                      "beta0", "beta1", "sigmaeps"),
                      labels = c("theta[1]", "theta[2]", "theta[3]", 
                                 "beta[0]", "beta[1]", "sigma[epsilon]"))


# Time series with true values
p.tracechainuni=ggplot(data=rr,aes(x=time,y=value))+
  geom_line()+
      # Add vlines using separate dataset
  geom_hline(data = vline_data, aes(yintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +

  facet_wrap(~parameter,scales="free_y")+
  theme_minimal()+
  ggtitle("Gibbs sampler chains for toy example (one phenotype)")

p.histuni=ggplot(data=rr,aes(x=value))+
  geom_histogram()+
      # Add vlines using separate dataset
  geom_vline(data = vline_data, aes(xintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +
  facet_wrap(~parameter,scales="free",label=label_parsed)+ylab("Counts")+xlab("Parameter value")+
  theme_minimal()+
  ggtitle("Posterior approximation for toy example parameters (one phenotype)")

print(p.histuni)



## ----include=FALSE----------------------------------------------
res.gibbsapproxuni$theta1 %>% quantile(c(0.025,0.975)) %>% diff()
res.gibbsapproxuni$theta2 %>% quantile(c(0.025,0.975)) %>% diff()
res.gibbsapproxuni$theta3 %>% quantile(c(0.025,0.975)) %>% diff()


## ---------------------------------------------------------------
pattern.on.phenotypes=c(1,-1,1)
dat2=data.frame(phenotype1=data$phenotype) %>%
  mutate(phenotype2=phenotype1+pattern.on.phenotypes,
         phenotype3=phenotype1-pattern.on.phenotypes,
         phenotype4=c(20,25,19))
dat2$cage=data$cage


## ----results='asis'---------------------------------------------
# Create a table using kable with escape = FALSE so that LaTeX code is not escaped
kable(dat2, format = "latex", escape = FALSE,
      caption = "\\label{tab:datatoyexmulti} Toy example with multiple phenotypes (columns) for the three subjects (rows) along with their Chronological Age (last column)",col.names = c("$\\boldsymbol{y}_{j=1}$","$\\boldsymbol{y}_{j=2}$","$\\boldsymbol{y}_{j=3}$","$\\boldsymbol{y}_{j=4}$","Chronological Age ($\\boldsymbol c$)"),aling="c",booktabs = TRUE,row.names = TRUE) %>%
  kable_styling(latex_options = c("HOLD_position"))


## ---------------------------------------------------------------
rmpostfullcondtheta <- function(J,n,y,x,beta0, beta1,sigmaeps,sigma.theta=7.5){
  mutheta=(beta1/sigmaeps^2)*(t(y)-beta0-beta1%*%t(x))
  sigma2theta=(1/sigma.theta^2+beta1^2/sigmaeps^2)
  ss=sum(sigma2theta)
  mus=colSums(mutheta)/ss
  rnorm(n,mean=mus,sd=sqrt(ss))
}


# Define the full conditional posterior for beta0 and beta1
rmpostfullcondbeta0 <- function(J,n,beta1, y,x, sigmaeps,theta){
  aa=(1/(n*sigmaeps^2))
  mu0=aa*rowSums(t(y)-beta1%*%t(x+theta))
  sig0=sigmaeps^2/n
  beta0j=rnorm(J,mean=mu0,sd=sqrt(sig0))
  beta0j
}

rmpostfullcondbeta1 <- function(J,n,beta0, y,x, sigmaeps,theta){
  sst=sum((x+theta)^2)
  mu1=colSums(t((t(y)-beta0))*(x+theta))/sst
  sig1=sigmaeps^2/sst
  beta1j=rnorm(J,mean=mu1,sd=sqrt(sig1))
  beta1j
}

# Define the full conditional posterior for sigma.epsilon^2
rmpostfullcondsigmaeps2 <- function(J,n,y,x, beta0,beta1,theta,
                                    sigmaeps.a=0.0001,sigmaeps.b=0.0001){
  an=(n/2+sigmaeps.a)
  bn=0.5*(rowSums((t(y)-beta0-beta1%*%t(x+theta))^2)+sigmaeps.b)
  sigmaeps2j=rinvgamma(J,shape = an,scale = bn)
  sigmaeps2j
}




## ---------------------------------------------------------------

gibbsapproxmulti.ind=function(R=10000,burnin=1000,data,sigmaeps.a=0.0001,
                              sigmaeps.b=0.0001,thinning=10){
  J=ncol(data)-1
  y=data[,1:J]
  x=data$cage
  n=nrow(y)
  # Initialize the chain
  theta.init=rep(0,n)
  beta0.init=rep(0,J)
  beta1.init=rep(1,J)
  sigmaeps2.init=rep(1,J)
  mchain=matrix(NA,nrow=R+burnin,ncol=n+3*J)
  colnames(mchain)=c(paste0("theta",1:n),
                    paste0("beta0",1:J),paste0("beta1",1:J),
                    paste0("sigmaeps2",1:J))
  # Define where to find parameters in chain
  idtheta=1:n
  idbeta0=(n+1):(n+J)
  idbeta1=(n+J+1):(n+2*J)
  idsigmaeps2=(n+2*J+1):(n+3*J)
  mchain[1,]=c(theta.init,beta0.init,beta1.init,sigmaeps2.init)
  for(r in 2:(R+burnin)){
    # Sample Theta 
    mchain[r,idtheta]=rmpostfullcondtheta(J,n,y,x,
                                         beta0=mchain[r-1,idbeta0],
                                         beta1=mchain[r-1,idbeta1],
                                         sigmaeps=sqrt(mchain[r-1,idsigmaeps2]))
    # Sample beta0 and beta1
    mchain[r,idbeta0]=rmpostfullcondbeta0(J,n,beta1=mchain[r-1,idbeta1],y=y,x=x,
                                         sigmaeps=sqrt(mchain[r-1,idsigmaeps2]),
                                         theta=mchain[r,idtheta])
    
    mchain[r,idbeta1]=rmpostfullcondbeta1(J,n,beta0=mchain[r,idbeta0],y=y,x=x,
                                         sigmaeps=sqrt(mchain[r-1,idsigmaeps2]),
                                         theta=mchain[r,idtheta])
    # Sample sigmaeps
    mchain[r,idsigmaeps2]=rmpostfullcondsigmaeps2(J,n,y=y,x=x,
                                                  beta0=mchain[r,idbeta0],
                                                  beta1=mchain[r,idbeta1],
                                                  theta=mchain[r,idtheta],
                                                  sigmaeps.a=sigmaeps.a,
                                                  sigmaeps.b=sigmaeps.b)
  }
  idx_keep <- seq(burnin + 1, R + burnin, by = thinning)
  return(mchain[idx_keep,])
}


if (file.exists("res_gibbsapproxmulti_ind_2.RData")) {
  load("res_gibbsapproxmulti_ind_2.RData")
} else {
  res.gibbsapproxmulti.ind.2 <- gibbsapproxmulti.ind(data = dat2 %>% 
    select(-phenotype3,-phenotype4))
  save(res.gibbsapproxmulti.ind.2, file = "res_gibbsapproxmulti_ind_2.RData")
}

if (file.exists("res_gibbsapproxmulti_ind_3.RData")) {
  load("res_gibbsapproxmulti_ind_3.RData")
} else {
  res.gibbsapproxmulti.ind.3 <- gibbsapproxmulti.ind(data = dat2 %>% 
    select(-phenotype4))
  save(res.gibbsapproxmulti.ind.3, file = "res_gibbsapproxmulti_ind_3.RData")
}

if (file.exists("res_gibbsapproxmulti_ind_4.RData")) {
  load("res_gibbsapproxmulti_ind_4.RData")
} else {
  res.gibbsapproxmulti.ind.4 <- gibbsapproxmulti.ind(data = dat2)
  save(res.gibbsapproxmulti.ind.4, file = "res_gibbsapproxmulti_ind_4.RData")
}



## ---------------------------------------------------------------
J=4
trueval=c(theta,beta0=rep(0,J),beta1=rep(1,J),sigmaeps=rep(1,J))

rr=res.gibbsapproxmulti.ind.4 %>% data.frame() %>% 
  mutate(time=1:nrow(res.gibbsapproxmulti.ind.4)) %>% 
  mutate(across(starts_with("sigmaeps2"), ~ sqrt(.)))%>% 
  rename_with(~ gsub("^sigmaeps2", "sigmaeps", .), starts_with("sigmaeps2"))

names(trueval)=colnames(rr %>% select(-time))

rr=rr %>% reshape2::melt(.,id.var="time",
                         value.name="value",variable.name="parameter")


# Create a separate dataset for vlines with true values of parameters
vline_data <- data.frame(
  parameter = names(trueval), # Ensure factors match facet levels
  xintercept = trueval
)

# Time series with true values
p.tracechainmulti=ggplot(data=rr,aes(x=time,y=value))+
  geom_line()+
      # Add vlines using separate dataset
  geom_hline(data = vline_data, aes(yintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +

  facet_wrap(~parameter,scales="free_y")+
  theme_minimal()+
  ggtitle("Gibbs sampler")

p.postmulti=ggplot(data=rr,aes(x=value))+
  geom_histogram()+
      # Add vlines using separate dataset
  geom_vline(data = vline_data, aes(xintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +
  facet_wrap(~parameter,scales="free")+
  theme_minimal()+
  ggtitle("Posterior approximation")



## ----fig.cap="Approximation of posterior for $\\theta_1$, $\\theta_2$ and $\\theta_3$ in the toy example $J_N=1,2,3,4$.  \\label{fig-toyexpthetavarjn}"----

rr1=res.gibbsapproxuni %>% select(starts_with("theta")) %>% 
  mutate(J=1,time=1:nrow(res.gibbsapproxuni))
rr2=res.gibbsapproxmulti.ind.2 %>% data.frame() %>% 
  select(starts_with("theta")) %>% 
  mutate(J=2,time=1:nrow(res.gibbsapproxmulti.ind.2)) 
rr3=res.gibbsapproxmulti.ind.3 %>% data.frame() %>% 
  select(starts_with("theta")) %>% 
  mutate(J=3,time=1:nrow(res.gibbsapproxmulti.ind.3)) 
rr4=res.gibbsapproxmulti.ind.4 %>% data.frame() %>% 
  select(starts_with("theta")) %>%
  mutate(J=4,time=1:nrow(res.gibbsapproxmulti.ind.4)) 

rr = rbind(rr1,rr2,rr3,rr4) %>%
  reshape2::melt(.,id.var=c("time","J"),value.name="value",
                 variable.name="parameter") 


rr$parameter = factor(rr$parameter, levels = c("theta1", "theta2", "theta3"),
                      labels = c("theta[1]", 
                              "theta[2]", 
                              "theta[3]"))

vline_data$parameter = factor(vline_data$parameter, levels = c("theta1", "theta2", "theta3"),
                      labels = c("theta[1]", 
                              "theta[2]", 
                              "theta[3]"))

# Comparison Graph
ggplot(data=rr,aes(x=factor(J),y=value))+
  geom_boxplot(outliers = FALSE)+
  # Add vlines using separate dataset
  geom_hline(data = vline_data[1:nrow(data),], aes(yintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +
  facet_wrap(~parameter, labeller = label_parsed)+xlab(TeX("Number of phenotypes $J_N$"))+ylab(TeX("$\\theta$"))+
  theme_minimal()+
  ggtitle(TeX("Posterior distribution of $\\theta$ in the toy example (one to four phenotypes)"))



## ----fig.cap="Approximation of posterior SNR for each model, $\\beta^2_{1j}/\\sigma^2_{\\epsilon_{ij}}$ along with the overall $\\sum_{j=1}^{J_N} \\beta^2_{1j}/\\sigma^2_{\\epsilon_{ij}}$. \\label{fig-toyexsnrjn}"----

rrsnr=res.gibbsapproxmulti.ind.4 %>% data.frame() %>%
    select(starts_with(c("sigmaeps","beta1"))) %>%
    mutate(across(starts_with("beta1"),
                  ~ .^2 /get(str_replace(cur_column(), "beta1", "sigmaeps2")),
                  .names = "snr_{.col}")) %>% select(starts_with("snr")) %>%
  mutate(time=1:nrow(res.gibbsapproxmulti.ind.4)) %>%
  rowwise() %>%
  mutate(overallSNR = sum(c_across(starts_with("snr_")))) %>%
  ungroup() %>%
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")

  rrsnr$parameter=factor(rrsnr$parameter,levels=c("snr_beta11","snr_beta12","snr_beta13","snr_beta14","overallSNR"),
                      labels=c("$y_{\\cdot 1}$","$y_{\\cdot 2}$","$y_{\\cdot 3}$",
                               "$y_{\\cdot 4}$","Overall SNR"))

# Comparison Graph
p5=ggplot(data=rrsnr,aes(y=value,x=parameter))+
  geom_boxplot(outliers=FALSE)+
  theme_minimal()+ylab(TeX("SNR=$\\beta_{1j}^2 / \\sigma_{\\epsilon_{ij}}^2$"))+
  xlab("Phenotype")+scale_x_discrete(labels = function(x) TeX(x))+
  ggtitle("Posterior distribution of SNR in the toy example (four phenotypes)")

print(p5)


## ---------------------------------------------------------------
vars_etabiol_AFFLY67 <- data.frame(
  numero = 1+c(2, 13, 69, 70, 104, 100, 66, 4, 3,1),
  nome_completo = c(
    "WEIGHT",
    "STATURE",
    "CHEST CIRC",
    "WAIST CIRC-OMPHALION",
    "ARM CIRC-BICEPS RELAXED",
    "CALF CIRC",
    "NECK CIRC OVER LARYNX",
    "TRICEPS SKF",
    "SUBSCAPULAR SKF",
    "AGE"
  ),
  abbreviazione_R = c(
    "weight",
    "stature",
    "chest_circ",
    "waist_circ",
    "arm_bic_circ",
    "calf_circ",
    "neck_circ",
    "triceps_skinfold",
    "subscap_skinfold",
    "cage"
  )
)


datex1=USAFSurvey[,vars_etabiol_AFFLY67$numero] %>% log()
colnames(datex1)=vars_etabiol_AFFLY67$abbreviazione_R

datex1=datex1 %>% mutate(cage=exp(cage)/12)

pc=princomp(datex1 %>% select(-cage),cor=TRUE)
datex1pc=pc$scores %>% scale() %>% as.data.frame() %>%
  mutate(cage=datex1$cage)

datex1scaled=scale(datex1 %>% select(-cage)) %>% 
  as.data.frame() %>% mutate(cage=datex1$cage)


## ---------------------------------------------------------------
useddat=datex1pc
if (file.exists("res.gibbsapproxmulti.ind.RData")) {
  load("res.gibbsapproxmulti.ind.RData")
} else {
  res.gibbsapproxmulti.ind=gibbsapproxmulti.ind(data=useddat)
  save(res.gibbsapproxmulti.ind, file = "res.gibbsapproxmulti.ind.RData")
}

rr=res.gibbsapproxmulti.ind %>% as.data.frame() %>% select(starts_with("theta"))
post.bioagemean=apply(rr,2,mean)+useddat$cage
post.bioage.sd=apply(rr,2,sd)
rrs=res.gibbsapproxmulti.ind %>% data.frame() %>% 
  mutate(time=1:nrow(res.gibbsapproxmulti.ind)) %>% 
  mutate(across(starts_with("sigmaeps2"), ~ sqrt(.)))%>% 
  rename_with(~ gsub("^sigmaeps2", "sigmaeps", .), starts_with("sigmaeps2")) %>%
  select("time",starts_with("sigmaeps")) %>%
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")

rrsnr=res.gibbsapproxmulti.ind %>% data.frame() %>%
    select(starts_with(c("sigmaeps","beta1"))) %>%
    mutate(across(starts_with("beta1"),
                  ~ .^2 /get(str_replace(cur_column(), "beta1", "sigmaeps2")),
                  .names = "snr_{.col}")) %>% select(starts_with("snr")) %>%
  mutate(time=1:nrow(res.gibbsapproxmulti.ind)) %>%
  rowwise() %>%
  mutate(overallSNR = sum(c_across(starts_with("snr_")), na.rm = TRUE)) %>%
  ungroup() %>%
  
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")


## ----fig.cap="Approximation of posterior SNR for each model for each PC in the USAF 1967 data along with the overall SNR. \\label{fig-usaf1967exsnr}"----

  rrsnr$parameter=factor(rrsnr$parameter,labels=c(paste("PC",1:9,sep=""),"Overall SNR"))

# Comparison Graph
p.snrex1=ggplot(data=rrsnr,aes(y=value,x=parameter))+
  geom_boxplot(outliers = FALSE)+
  theme_minimal()+ylab(TeX("SNR=$\\beta_{1j}^2 / \\sigma_{\\epsilon_{ij}}^2$"))+
  xlab("Phenotype")+scale_x_discrete(labels = function(x) TeX(x))+
  ggtitle("Posterior distribution of SNR in the USAF 1967 data set")

print(p.snrex1)


## ----fig.cap="Approximation of posterior for the USAF 1967 data set. The posterior mean and the 95\\% credible interval for the BA are reported in comparison with the chronological age. \\label{fig-usaf1967postba}"----

rrs$parameter=factor(rrs$parameter,
                      labels=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9"))
p.postsigmaex1=ggplot(data=rrs,aes(y=value,x=parameter))+
  geom_boxplot(outliers = FALSE)+
  theme_minimal()+ylab(TeX("$\\sigma_{\\epsilon_{ij}}^2$"))+xlab("Principal Component")+
  ggtitle(TeX("Posterior approximation of $\\sigma_{\\epsilon_{ij}}^2$ in the USAF 1967 data set"))



p.biovschronex1=ggplot(data=data.frame(cage=useddat$cage,bioage=post.bioagemean),
          aes(x=cage,y=bioage))+
  geom_errorbar(aes(ymin=post.bioagemean-2*post.bioage.sd,
                    ymax=post.bioagemean+2*post.bioage.sd),color="grey")+
  geom_point()+
  geom_abline(intercept = 0,slope=1,color="red")+
  xlab("Chronological age")+
  ylab("Biological age")+
  ggtitle("Biological age estimation")

print(p.biovschronex1)



## ---------------------------------------------------------------
rmpostfullcondtheta <- function(J, n, y, x, beta, Sigma, sigma.theta = 7.5) {
  # Compute the inverse of Sigma once (Sigma is JxJ)
  Sigma_inv <- solve(Sigma)
  # Extract beta0 and beta1 from the beta matrix
  beta0 <- beta[, 1]
  beta1 <- beta[, 2]
  theta <- numeric(n)
  
  for (i in 1:n) {
    # For the i-th observation:
    # y[i,] is a vector of length J and x[i] is c_i (a scalar)
    # Compute the residual: y_i - beta0 - beta1 * c_i
    resid_i <- y[i, ] - beta0 - beta1 * x[i]
    # Numerator: beta1' * Sigma^{-1} * (y_i - beta0 - beta1*c_i)
    num <- t(beta1) %*% Sigma_inv %*% resid_i
    # Denom: beta1' * Sigma^{-1} * beta1 + 1/sigma.theta^2
    denom <- t(beta1) %*% Sigma_inv %*% beta1 + 1/(sigma.theta^2)
    mean_theta <- num / denom
    var_theta  <- 1 / denom
    theta[i] <- rnorm(1, mean = mean_theta, sd = sqrt(var_theta))
  }
  return(theta)
}

rmpostfullcondbeta <- function(J, n, y, x, theta, Sigma) {
  # Construct the design matrix X (n x 2)
  X <- cbind(rep(1, n), x + theta)
  XtX <- t(X) %*% X       # 2 x 2 matrix
  invXtX <- solve(XtX)    # (X^T X)^{-1}
  
  # Compute the posterior mean for beta (J x 2)
  beta_hat <- t(y) %*% X %*% invXtX  # dimensions: J x 2
  
  # To sample from MN_{J x 2}(beta_hat, Sigma, invXtX), we can do:
  #   beta = beta_hat + L_Sigma %*% Z %*% t(L_V)
  # where L_Sigma is the Cholesky factor of Sigma (upper triangular)
  # and L_V is the Cholesky factor of invXtX.
  L_Sigma <- chol(Sigma)       # J x J
  L_V <- chol(invXtX)          # 2 x 2
  Z <- matrix(rnorm(J * 2), nrow = J, ncol = 2)
  beta_sample <- beta_hat + L_Sigma %*% Z %*% t(L_V)
  return(beta_sample)
}

rmpostfullcondSigma <- function(J, n, y, x, beta, theta, eta0, S0) {
  # Construct design matrix X (n x 2)
  X <- cbind(rep(1, n), x + theta)
  # Initialize the sum of squares matrix
  S <- matrix(0, nrow = J, ncol = J)
  for (i in 1:n) {
    # For each i, compute the mean vector: mu_i = beta %*% X[i, ]
    mu_i <- beta %*% X[i, ]
    resid <- y[i, ] - as.vector(mu_i)
    S <- S + resid %*% t(resid)
  }
  S_post <- S0 + S
  df_post <- eta0 + n
  # Sample Sigma from the Inverse-Wishart distribution
  Sigma_sample <- rinvwishart(df_post, S_post)
  return(Sigma_sample)
}


## ---------------------------------------------------------------
gibbsapproxmulti=function(R=10000,burnin=1000,data,eta0=1,priorS0=1,thinning=10){
  J=ncol(data)-1
  S0=diag(rep(priorS0,J))
  y=as.matrix(data[,1:J])
  x=data$cage
  n=nrow(y)
  # Initialize the chain
  beta.init=t(array(c(0,0),dim=c(2,J)))
  colnames(beta.init)=c("beta0","beta1")
  theta.init=rep(0,n)
  names(theta.init)=paste0("theta",1:n)
  Sigma.init=diag(J)
  dimnames(Sigma.init)=list(colnames(y),colnames(y))
  chain.beta=array(NA,dim=c(J,2,R+burnin))
  dimnames(chain.beta)=list(colnames(y),c("beta0","beta1"),NULL)
  chain.theta=matrix(NA,nrow=R+burnin,ncol=n)
  colnames(chain.theta)=paste0("theta",1:n)
  chain.Sigma=array(NA,dim=c(J,J,R+burnin))
  dimnames(chain.Sigma)=list(colnames(y),colnames(y),NULL)
  chain.beta[,,1]=beta.init
  chain.theta[1,]=theta.init
  chain.Sigma[,,1]=Sigma.init
  for(r in 2:(R+burnin)){
  # Sample theta given current beta and Sigma:
  chain.theta[r,] <-
    rmpostfullcondtheta(J, n, y, x, beta=chain.beta[,,r-1],
                        Sigma=chain.Sigma[,,r-1])
  # Sample beta given current theta and Sigma:
  chain.beta[,,r] <- 
    rmpostfullcondbeta(J, n, y, x, theta = chain.theta[r,],
                       Sigma=chain.Sigma[,,r-1])
  # Sample Sigma given current beta and theta:
  chain.Sigma[,,r] <- 
    rmpostfullcondSigma(J, n, y, x, beta=chain.beta[,,r],
                        theta = chain.theta[r,], eta0, S0)
  }
  idx_keep <- seq(burnin + 1, R + burnin, by = thinning)
  chain.beta=chain.beta[,,idx_keep]
  chain.theta=chain.theta[idx_keep,]
  chain.Sigma=chain.Sigma[,,idx_keep]
  
  return(list(beta=chain.beta,
              theta=chain.theta,
              Sigma=chain.Sigma))
}


## ---------------------------------------------------------------
if (file.exists("res_gibbsapproxmulti.RData")) {
  load("res_gibbsapproxmulti.RData")
} else {
  res.gibbsapproxmulti <- gibbsapproxmulti(data = datex1scaled)
  save(res.gibbsapproxmulti, file = "res_gibbsapproxmulti.RData")
}

rr=res.gibbsapproxmulti$beta
rr=rr %>% matrix(aperm(., c(3, 1, 2)), nrow = dim(rr)[3], 
                 ncol = prod(dim(rr)[1:2])) %>% data.frame()

# Extract row names (variable names)
var_names <- rownames(res.gibbsapproxmulti$beta[,,1])

# Define column suffixes corresponding to beta0 and beta1
suffixes <- colnames(res.gibbsapproxmulti$beta[,,1])

# Generate new column names: "variable_beta0" and "variable_beta1"
col_names <- as.vector(outer(var_names, suffixes, paste, sep = "_"))

# Assign column names to rr
colnames(rr) <- col_names

rr=rr %>% mutate(time=1:nrow(rr)) %>% 
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")



## ----fig.cap="Approximation of posterior distribution of the correlation among the phenotypes in the linear model for the USAF 1967 data set. \\label{fig-usaf1967postSigma}"----


rr=res.gibbsapproxmulti$Sigma
rr=apply(rr,3,cov2cor)
rr=array(rr, dim = dim(res.gibbsapproxmulti$Sigma),
         dimnames = dimnames(res.gibbsapproxmulti$Sigma))

rr=rr %>% matrix(aperm(., c(3, 1, 2)), nrow = dim(rr)[3], 
                 ncol = prod(dim(rr)[1:2])) %>% data.frame()

# Extract row names (variable names)
var_names <- rownames(res.gibbsapproxmulti$Sigma[,,1])

# Define column suffixes corresponding to beta0 and beta1
suffixes <- colnames(res.gibbsapproxmulti$Sigma[,,1])

# Generate new column names: "variable_beta0" and "variable_beta1"
col_names <- as.vector(outer(var_names, suffixes, paste, sep = "."))

# Assign column names to rr
colnames(rr) <- col_names

rr=rr %>% mutate(time=1:nrow(rr)) %>% 
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")

rr$par1=strsplit(as.character(rr$parameter),".",fixed=TRUE) %>% sapply(.,function(x) x[1])
rr$par2=strsplit(as.character(rr$parameter),".",fixed=TRUE) %>% sapply(.,function(x) x[2])

rr = rr %>% filter(par1 != par2)

tt=table(rr$par1,rr$par2)
dd=data.frame(par1=colnames(tt)) %>% mutate(par2=par1)
nn=which(lower.tri(tt),arr.ind = TRUE)
exclude=cbind(dd[nn[,1],1],dd[nn[,2],2])

exclude=c(apply(exclude,1,function(x) which((rr$par1==x[1]) & (rr$par2==x[2]))))

rr=rr[-exclude,]

# Posterior distributions
p.postcorrusaf1967 <- ggplot(data = rr, aes(x = value)) +
#  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  geom_density() +  # Remove outliers
#  geom_histogram() +  # Remove outliers
#  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  xlab("Correlation") +  # Remove x-axis label
  facet_grid(par1 ~ par2) + xlim(c(0,0.9))+
  theme_minimal() +
  ggtitle("Posterior distributions of correlations among phenotypes in the USAF 1967 data set.") +
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1),  # Rotate x-axis labels 90 degrees
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.y = element_blank(),   # Remove y-axis text (numbers)
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    strip.text = element_text(size = 6),  # Make facet titles smaller
    plot.title = element_text(size = 12),  # Adjust main title size
    strip.text.y = element_text(angle = 0)  # Set y-axis facet labels horizontal
 
  )

print(p.postcorrusaf1967)


# # Time series
# p.postcorrusaf1967trace=ggplot(data=rr,aes(x=time,y=value))+
#   geom_line()+
#   facet_wrap(~parameter)+ylab(c(-1,1))+
#   theme_minimal()+
#   ggtitle("Gibbs sampler for Correlations")



## ----fig.cap="Approximation of the posterior of BA for the USAF 1967 data set when accounting for dependence among phenotypes and a linear model. The posterior mean and the 95\\% credible interval for the BA are reported in comparison with the chronological age. \\label{fig-usaf1967postbadep}"----


rr=res.gibbsapproxmulti$theta %>% as.data.frame()
post.bioagemean=apply(rr,2,mean)+datex1$cage
post.bioage.sd=apply(rr,2,sd)

p.usaf1967postba.multi=ggplot(data=data.frame(cage=datex1$cage,bioage=post.bioagemean),
          aes(x=cage,y=bioage))+
  geom_errorbar(aes(ymin=post.bioagemean-2*post.bioage.sd,
                    ymax=post.bioagemean+2*post.bioage.sd),color="grey")+
  geom_point()+
  geom_abline(intercept = 0,slope=1,color="red")+
  xlab("Chronological age")+
  ylab("Biological age")+
  ggtitle("Biological age estimation")

print(p.usaf1967postba.multi)



## ---------------------------------------------------------------
# This calculates the observational log-likelihood 
# with mixed phenotypes and covariates X
calcllk=function(dlmodel,dlmodinput,obspheno,invSigma){
  preds=predict(dlmodel,dlmodinput,verbose = 0)
  mus=yN=NULL
  contrib.multinom=NULL
  for(i in 1:length(preds)){
    rr=preds[[i]]
    if(ncol(rr)==1){
      mus=cbind(mus,rr[,1])
      yN=cbind(yN,obspheno[[i]])}else{
        ccmn=rowSums(log(rr)*obspheno[[i]])
        contrib.multinom=cbind(contrib.multinom,ccmn)
      }
  }
  multinompart=ifelse(is.null(contrib.multinom),0,rowSums(contrib.multinom))    
  normalpart=-0.5*rowSums(((yN - mus) %*% invSigma) * (yN - mus))
  return(list(ll=multinompart+normalpart,mus=mus))
}


mhgibbsdl <- function(R = 1000, burnin = 100, data,obspheno,dlmodinput, dlmodel,
                      indJN,indJC, eta0 = 1, priorS0 = 1, 
                      thinning = 1) {

  JN <- length(indJN)      # number of numerical phenotype columns
  yN=as.matrix(data[,indJN])
  S0 <- diag(rep(priorS0, JN))
  x <- data$cage
  n <- nrow(data)
  # Fixed prior variance for theta: sigma_theta = 7.5
  sigma_theta <- 7.5
  
  # Initialize the chain
  theta <- rep(0, n)         # latent offsets
  Sigma <- diag(JN)           # covariance matrix
  chain.theta <- matrix(NA, nrow = R + burnin, ncol = n)
  chain.Sigma <- array(NA, dim = c(JN, JN, R + burnin),
                       dimnames = list(colnames(yN),colnames(yN),NULL))
  
  chain.theta[1, ] <- theta
  chain.Sigma[,  1] <- Sigma
  
  # Tuning parameter for MH proposal for theta
  tau <- 0.1
  idcage=ncol(dlmodinput)
  dlmodinput[,idcage]=x+theta
  llcurrent=calcllk(dlmodel,dlmodinput,obspheno,invSigma = solve(Sigma))
  f_current=llcurrent$mus
  input.prop=dlmodinput
  
  for (r in 2:(R + burnin)) {
    invSigma=solve(Sigma)
    # Update theta for each subject using MH
    theta_prop <- theta + rnorm(n, mean = 0, sd = tau)
    # Predict the phenotype means using dlmodel for proposed theta
    input.prop[,idcage]=x+theta_prop
    llprop=calcllk(dlmodel,input.prop,obspheno,invSigma = invSigma)
    prior_current <- (- (theta^2) / (2 * sigma_theta^2))
    prior_prop    <- (- (theta_prop^2) / (2 * sigma_theta^2))
    alpha <- (llprop$ll + prior_prop) - (llcurrent$ll + prior_current)
    ii=which(log(runif(n))<alpha)
    theta[ii]=theta_prop[ii]
    chain.theta[r, ] <- theta
    llcurrent$ll[ii]=llprop$ll[ii]
    llcurrent$mus[ii,]=llprop$mus[ii,]

    # Update Sigma using its full conditional (Inverse-Wishart)
    # Compute the sum of squared deviations: 
    # S_update = S0 + sum_i (y_i - f(c_i+theta_i))(y_i - f(c_i+theta_i))^T
    S_update <- S0
    diff <- yN - llcurrent$mus
    for (i in 1:n) S_update <- S_update + diff[i,] %*% t(diff[i,])
    df_post <- eta0 + n
    Sigma <- rinvwishart(df_post, S_update)
    chain.Sigma[,  r] <- Sigma
  }
  
  # Apply burn in and thinning
  idx_keep <- seq(burnin + 1, R + burnin, by = thinning)
  chain.theta <- chain.theta[idx_keep, ]
  chain.Sigma <- chain.Sigma[,  idx_keep]
  
  return(list(theta = chain.theta, Sigma = chain.Sigma))
}



## ----warning=FALSE----------------------------------------------
create.dlmodel=function(data,indJN,indJC,indX=NULL){
  # Assume 'data' is a data.frame that includes:
  # - The indJ columns: the phenotypes (outputs)
  # - the indX columns: covariates of phenotypes
  # - A column 'cage' giving the chronological age (input)
  
  # Determine the number of phenotype outputs (J)
  J <- length(indJN)+length(indJC)
  # Identify the output variable names (first J columns)
  output_vars <- names(data)[c(indJN,indJC)]
  
  if(!is.null(indX)){
  X = data[indX] 
  X = cbind(X %>% select(where(is.numeric)),X %>% 
              select(-where(is.numeric)) %>%
              dummy_cols(remove_first_dummy = FALSE, 
                         remove_selected_columns = TRUE)) %>% as.matrix()
  }
  # Prepare lists to store output layers, loss functions, and metrics per output
  output_layers <- list()
  loss_list <- list()
  metrics_list <- list()
  
  # Use index to ensure unique names even if output_vars have duplicates
  for (i in seq_along(output_vars)) {
    var <- output_vars[i]
    if (is.factor(data[[var]])) {
      # For categorical outputs: determine number of classes 
      # and use softmax activation
      n_classes <- length(levels(data[[var]]))
      # Unique name for the output layer and metric
      output_name <- paste0(var, "_output_", i)
      metric_name <- paste0(var, "_accuracy_", i)
      
      output_layers[[var]] <- layer_dense(units = n_classes,
                                          activation = "softmax", 
                                          name = output_name)
      loss_list[[var]] <- "categorical_crossentropy"
      metrics_list[[var]] <- metric_accuracy(name = metric_name)
    } else {
      # For continuous outputs: use linear activation
      output_name <- paste0(var, "_output_", i)
      metric_name <- paste0(var, "_mae_", i)
      
      output_layers[[var]] <- layer_dense(units = 1, activation = "linear", 
                                          name = output_name)
      loss_list[[var]] <- "mse"
      metrics_list[[var]] <- metric_mean_squared_error(name = metric_name)
    }
  }
  
  # Build the base network using the functional API:
  # Input layer: scalar input (chronological age)
#  input_layer <- layer_input(shape = 2, name = "cage_input")
  # Input: covariates and chronological age as a matrix.
  if(!is.null(indX)){
  x_train <- cbind(X,cage=as.matrix(data$cage)) %>% as.matrix()}else{
    x_train = as.matrix(data$cage)
  }

  input_layer <- layer_input(shape = ncol(x_train), name = "input_layer")
  # Three hidden dense layers with ReLU activation
  # Three hidden dense layers with ReLU activation
  dense1 <- input_layer %>% layer_dense(units = 3*J, activation = "relu")
  dense2 <- dense1 %>% layer_dense(units = 3*J, activation = "relu")
  dense2d <- dense2 %>% layer_dropout(rate = 0.3)
  densefinal <- dense2d %>% layer_dense(units = 3*J, activation = "relu")
  
  # Branch out to the output layers for each phenotype variable.
  # Each output branch takes the common hidden representation.
  output_list <- lapply(names(output_layers), function(var) {
    output_layers[[var]](densefinal)
  })
  names(output_list) <- names(output_layers)
  
  # Create the Keras model with the specified input and outputs.
  dlmodel <- keras_model(inputs = input_layer, outputs = output_list)
  
  # Compile the model with the appropriate loss functions and optimizer.
  dlmodel %>% compile(
    optimizer = "adam",
    loss = loss_list,
    metrics = metrics_list
  )
  
  # Prepare the training data.
  
  # Outputs: a named list, one element per phenotype.
  y_train <- list()
  for (var in output_vars) {
    if (is.factor(data[[var]])) {
      # Convert factor to one-hot encoded matrix
      n_classes <- length(levels(data[[var]]))
      y_train[[var]] <- to_categorical(as.integer(data[[var]]) - 1, 
                                       num_classes = n_classes)
    } else {
      y_train[[var]] <- as.matrix(data[[var]])
    }
  }
  return(list(dlmodel=dlmodel,data=data,
              x_train=x_train,
              y_train=y_train,
              indJN=indJN,
              indJC=indJC,indX=indX))
}


## ----warning=FALSE----------------------------------------------
tt=create.dlmodel(datex1,indJN = 1:8,indJC=NULL)

if (file.exists("dlmodel.tf")) {
  dlmodel=load_model_tf("dlmodel.tf")
} else {
  dlmodel=tt$dlmodel
  # Fit the model.
  history <- dlmodel %>% fit(
    x = tt$x_train,
    y = tt$y_train,
    epochs = 100,
    batch_size = 16,
    validation_split = 0.2,verbose = 0)
  save_model_tf(dlmodel,file="dlmodel.tf")
}





## ---------------------------------------------------------------
preds=matrix(unlist(predict(dlmodel,tt$x_train,verbose = 0)),
             ncol=ncol(tt$data)-1,nrow=nrow(tt$data))
colnames(preds)=colnames(tt$data[,-ncol(tt$data)])
obs=tt$data[,-ncol(tt$data)]
R2est=round(diag(cor(preds,obs)^2),3)


## ----results='asis'---------------------------------------------
row_names <- c("weight", "stature", "chest circumference", "waist circumference (omphalion)", "relaxed biceps arm circumference", "calf circumference", "neck circumference (over larynx)", "triceps skinfold", "subscapular skinfold")

# Creo un data frame e assegno i nomi delle righe
table_data <- data.frame(R2est)
rownames(table_data) <- row_names

# Genero la tabella in LaTeX tramite kable e kableExtra
kable(table_data, format = "latex", escape = FALSE,
      caption = "\\label{tab:usaf1967dlR2} The $R^2$ values for the fitted FFNN on the USAF 1967 data set.",
      align = "c", booktabs = TRUE, col.names = "$R^2$") %>%
  kable_styling(latex_options = c("HOLD_position"))



## ---------------------------------------------------------------
if (file.exists("res_gibbsapproxmulti_dl.RData")) {
  load("res_gibbsapproxmulti_dl.RData")
} else {
  res.gibbsapproxmulti.dl <- mhgibbsdl(
    R = 1000,
    burnin = 100, 
    data = tt$data,
    obspheno = tt$y_train,
    dlmodinput = tt$x_train,
    dlmodel = dlmodel,
    indJN = tt$indJN,
    indJC = tt$indJC,
    eta0 = 1,
    priorS0 = 1,
    thinning = 10
  )
  save(res.gibbsapproxmulti.dl, file = "res_gibbsapproxmulti_dl.RData")
}


## ---------------------------------------------------------------
rr=res.gibbsapproxmulti.dl$Sigma
rr=apply(rr,3,cov2cor)
rr=array(rr, dim = dim(res.gibbsapproxmulti$Sigma),
         dimnames = dimnames(res.gibbsapproxmulti$Sigma))

rr=rr %>% matrix(aperm(., c(3, 1, 2)), nrow = dim(rr)[3], 
                 ncol = prod(dim(rr)[1:2])) %>% data.frame()

# Extract row names (variable names)
var_names <- rownames(res.gibbsapproxmulti$Sigma[,,1])

# Define column suffixes corresponding to beta0 and beta1
suffixes <- colnames(res.gibbsapproxmulti$Sigma[,,1])

# Generate new column names: "variable_beta0" and "variable_beta1"
col_names <- as.vector(outer(var_names, suffixes, paste, sep = "."))

# Assign column names to rr
colnames(rr) <- col_names

rr=rr %>% mutate(time=1:nrow(rr)) %>% 
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")

rr$par1=strsplit(as.character(rr$parameter),".",fixed=TRUE) %>% sapply(.,function(x) x[1])
rr$par2=strsplit(as.character(rr$parameter),".",fixed=TRUE) %>% sapply(.,function(x) x[2])

rr = rr %>% filter(par1 != par2)

tt=table(rr$par1,rr$par2)
dd=data.frame(par1=colnames(tt)) %>% mutate(par2=par1)
nn=which(lower.tri(tt),arr.ind = TRUE)
exclude=cbind(dd[nn[,1],1],dd[nn[,2],2])

exclude=c(apply(exclude,1,function(x) which((rr$par1==x[1]) & (rr$par2==x[2]))))

rr=rr[-exclude,]

# Posterior distributions
p.postcorrusaf1967dl <- ggplot(data = rr, aes(x = value)) +
#  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  geom_density() +  # Remove outliers
#  geom_histogram() +  # Remove outliers
#  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  xlab("Correlation") +  # Remove x-axis label
  facet_grid(par1 ~ par2) + xlim(c(-0.5,0.9))+
  theme_minimal() +
  ggtitle("Posterior distributions of correlations among phenotypes in the USAF 1967 data set.") +
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1),  # Rotate x-axis labels 90 degrees
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.y = element_blank(),   # Remove y-axis text (numbers)
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    strip.text = element_text(size = 6),  # Make facet titles smaller
    plot.title = element_text(size = 12),  # Adjust main title size
    strip.text.y = element_text(angle = 0)  # Set y-axis facet labels horizontal
 
  )



## ----fig.cap="Approximation of the posterior of BA for the USAF 1967 data set when accounting for dependence among phenotypes and a the DL model. The posterior mean and the 95\\% credible interval for the BA are reported in comparison with the chronological age. \\label{fig-usaf1967postbadl}"----


rr=res.gibbsapproxmulti.dl$theta %>% as.data.frame()
post.bioagemean=apply(rr,2,mean)+datex1$cage
post.bioage.sd=apply(rr,2,sd)

p.usaf1967postba.dl=ggplot(data=data.frame(cage=datex1$cage,bioage=post.bioagemean),
          aes(x=cage,y=bioage))+
  geom_errorbar(aes(ymin=post.bioagemean-2*post.bioage.sd,
                    ymax=post.bioagemean+2*post.bioage.sd),color="grey")+
  geom_point()+
  geom_abline(intercept = 0,slope=1,color="red")+
  xlab("Chronological age")+
  ylab("Biological age")+
  ggtitle("Biological age estimation")

print(p.usaf1967postba.dl)



## ---------------------------------------------------------------
load("Datos/dat_con_Bioage.RData")
rownames(dat.used)=dat.used$n.inclusion
datamutua=dat.used %>% 
  select(-n.inclusion,-edad.biologica,-edad.diff,-edad.rango) %>%
  rename(cage=edad) %>%
  select(where(is.numeric), everything()) %>%
  select(-cage, cage)
data=datamutua
indJN=which(colnames(datamutua)%in%c("actividad.fisica.intensa.met",
                                     "ejercicio.fisico.moderado..dias.",
                                     "ejercicio.fisico.moderado..minutos.dia.",
                                     "cuanto.tiempo.real.dedico.a.caminar.en.los.ultimos.7.dias..minutos.",
                                     "sumatoria"))
indJC=which(colnames(datamutua)%in%c("IMC.fuera.normal","comorbilidades",
                                     "tiene.stress",
                                     "cuantas.horas.duerme.al.dia"))

indX=which(!(colnames(datamutua)%in%c("cage",
                                      colnames(datamutua)[c(indJN,indJC)])))

colnames(datamutua)=c(
  "Animal fat consumption",
  "Fish consumption",
  "Fish consumption score (9-point scale)",
  "Fruit consumption",
  "Intense physical activity (METs)",
  "Moderate physical activity (days/week)",
  "Moderate physical activity (minutes/day)",
  "Time spent walking in the last 7 days (minutes)",
  "Sum of lifestyle scores",
  "Sex",
  "BMI outside normal range",
  "Presence of comorbidities",
  "Presence of stress",
  "Sleep duration (hours/day)",
  "Educational level",
  "Consumption of distilled alcohol",
  "Consumption of fermented alcohol",
  "Gene rs8192678: Count of A",
  "Gene rs1726866: Count of a",
  "Gene rs660839: Count of A",
  "Gene rs480902: Count of A",
  "Gene rs5030980: Count of A",
  "Gene rs4994: Count of a",
  "Gene rs738409: Count of A",
  "Gene rs12033832: Count of a",
  "Gene rs11549467: Count of a",
  "Gene rs1137100: Count of a",
  "Gene rs571312: Count of a",
  "Gene rs11549465: Count of A",
  "Gene rs762551: Count of a",
  "Gene rs429358: Count of A",
  "Gene rs35874116: Count of A",
  "Gene rs713598: Count of A",
  "Gene rs6821591: Count of A",
  "Gene rs1042713: Count of a",
  "Gene rs17782313: Count of A",
  "Gene rs10246939: Count of A",
  "Gene rs7412: Count of A",
  "Gene rs4646116: Count of A",
  "Chronological age"
)

datamutuaused=datamutua[,c(indJN,indJC,indX)]



## ----warning=FALSE----------------------------------------------
tt=create.dlmodel(data,indJN=indJN,indJC=indJC,indX=indX)

if (file.exists("dlmodelmutua.tf")) {
  dlmodel=load_model_tf("dlmodelmutua.tf")
} else {
  dlmodel=tt$dlmodel
  # Fit the model.
  history <- dlmodel %>% fit(
    x = tt$x_train,
    y = tt$y_train,
    epochs = 1000,
    batch_size = 8,
    validation_split = 0.2,verbose = 0)
  save_model_tf(dlmodel,file="dlmodelmutua.tf")
}


## ---------------------------------------------------------------
preds=predict(dlmodel,tt$x_train,verbose = 0)
fittstat=NULL
for(i in 1:length(preds)){
  rr=preds[[i]]
  if(ncol(rr)==1){
    fittstat=c(fittstat,R2=cor(tt$y_train[[i]],rr[,1])^2)}else{
      fittstat=c(fittstat,probs.corr=mean(rr*tt$y_train[[i]]))}
}
fittstat=round(fittstat,2)
names(fittstat)=colnames(datamutua)[c(indJN,indJC)]


## ---------------------------------------------------------------
if (file.exists("res_gibbsmutua.RData")) {
  load("res_gibbsmutua.RData")
} else {
  res.gibbsmutua <- mhgibbsdl(
    R = 1000, 
    burnin = 100, 
    data = datamutua,
    obspheno = tt$y_train,
    dlmodinput = tt$x_train, 
    dlmodel = dlmodel,
    indJN = indJN,
    indJC = indJC, 
    eta0 = 1, 
    priorS0 = 1, 
    thinning = 10
  )
  save(res.gibbsmutua, file = "res_gibbsmutua.RData")
}


## ---------------------------------------------------------------
rr=res.gibbsmutua$Sigma
rr=apply(rr,3,cov2cor)
rr=array(rr, dim = dim(res.gibbsapproxmulti$Sigma),
         dimnames = dimnames(res.gibbsapproxmulti$Sigma))

rr=rr %>% matrix(aperm(., c(3, 1, 2)), nrow = dim(rr)[3], 
                 ncol = prod(dim(rr)[1:2])) %>% data.frame()

# Extract row names (variable names)
var_names <- rownames(res.gibbsapproxmulti$Sigma[,,1])

# Define column suffixes corresponding to beta0 and beta1
suffixes <- colnames(res.gibbsapproxmulti$Sigma[,,1])

# Generate new column names: "variable_beta0" and "variable_beta1"
col_names <- as.vector(outer(var_names, suffixes, paste, sep = "."))

# Assign column names to rr
colnames(rr) <- col_names

rr=rr %>% mutate(time=1:nrow(rr)) %>% 
  reshape2::melt(.,id.var="time",value.name="value",variable.name="parameter")

rr$par1=strsplit(as.character(rr$parameter),".",fixed=TRUE) %>% sapply(.,function(x) x[1])
rr$par2=strsplit(as.character(rr$parameter),".",fixed=TRUE) %>% sapply(.,function(x) x[2])

rr = rr %>% filter(par1 != par2)

tt=table(rr$par1,rr$par2)
dd=data.frame(par1=colnames(tt)) %>% mutate(par2=par1)
nn=which(lower.tri(tt),arr.ind = TRUE)
exclude=cbind(dd[nn[,1],1],dd[nn[,2],2])

exclude=c(apply(exclude,1,function(x) which((rr$par1==x[1]) & (rr$par2==x[2]))))

rr=rr[-exclude,]

# Posterior distributions
p.postcorrmutua <- ggplot(data = rr, aes(x = value)) +
#  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  geom_density() +  # Remove outliers
#  geom_histogram() +  # Remove outliers
#  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  xlab("Correlation") +  # Remove x-axis label
  facet_grid(par1 ~ par2) + xlim(c(-0.5,0.9))+
  theme_minimal() +
  ggtitle("Posterior distributions of correlations among phenotypes in the USAF 1967 data set.") +
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1),  # Rotate x-axis labels 90 degrees
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.y = element_blank(),   # Remove y-axis text (numbers)
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    strip.text = element_text(size = 6),  # Make facet titles smaller
    plot.title = element_text(size = 12),  # Adjust main title size
    strip.text.y = element_text(angle = 0)  # Set y-axis facet labels horizontal
 
  )



## ----fig.cap="Approximation of the posterior of BA for the SDA data using the general DL model. The posterior mean and the 95\\% credible interval for the BA are reported in comparison with the chronological age. \\label{fig:mutuapostba}"----


rr=res.gibbsmutua$theta %>% as.data.frame()
post.bioagemean=apply(rr,2,mean)+datex1$cage
post.bioage.sd=apply(rr,2,sd)

p.sdapostba=ggplot(data=data.frame(cage=datex1$cage,bioage=post.bioagemean),
          aes(x=cage,y=bioage))+
  geom_errorbar(aes(ymin=post.bioagemean-2*post.bioage.sd,
                    ymax=post.bioagemean+2*post.bioage.sd),color="grey")+
  geom_point()+
  geom_abline(intercept = 0,slope=1,color="red")+
  xlab("Chronological age")+
  ylab("Biological age")+
  ggtitle("Biological age estimation for SDA data set")

print(p.sdapostba)



## ----fig.cap="The posterior mean and the 95\\% credible interval for the BA are reported in comparison with the previous estimation of the BA for the individuals in the SDA data set. \\label{fig:mutuaprevbioage}"----


rr=res.gibbsmutua$theta %>% as.data.frame()
post.bioagemean=apply(rr,2,mean)+datamutua$`Chronological age`
post.bioage.sd=apply(rr,2,sd)

p.mutuacompbaest=ggplot(data=data.frame(prevbioage=dat.used$edad.biologica,bioage=post.bioagemean),aes(x=prevbioage,y=bioage))+
  geom_errorbar(aes(ymin=post.bioagemean-2*post.bioage.sd,ymax=post.bioagemean+2*post.bioage.sd),color="grey")+   geom_point()+
  geom_abline(intercept = 0,slope=1,color="red")+
  xlab("Previous estimation of Biological age")+
  ylab("Biological age")+
  ggtitle("Comparison of different estimations of the Biological age")

print(p.mutuacompbaest)



## ----fig.cap="Point estimation of BA in the inverted approach with simple linear regression.\\label{fig:invertedapplm}"----

olsmodinv=lm(cage~phenotype,data=datainv)
theta=datainv$theta
datainv <- data.frame(
  phenotype = c(datainv$phenotype),
  cage = c(datainv$cage),
  true.bioage = cage + theta,
  beta1m0=1,
  beta1m1=1/2,
  beta1m2=2,
  beta1mols=coef(olsmodinv)[2],
  beta0mols=coef(olsmodinv)[1],
  biom1 = datainv$phenotype * 1/2, 
  biom2 = datainv$phenotype * 2) %>%
  mutate(biomols= phenotype*beta1mols+beta0mols,
         sigmam0=sqrt(mean(phenotype^2)+mean((cage-beta1m0*phenotype)^2)),
         sigmam1=sqrt(mean(phenotype^2)+mean((cage-beta1m1*phenotype)^2)),
         sigmam2=sqrt(mean(phenotype^2)+mean((cage-beta1m2*phenotype)^2)),
         sigmay=sqrt(mean(phenotype^2)),
         sigmamols=sqrt(mean(phenotype^2)+mean((cage-beta0mols-beta1mols*phenotype)^2)))



modnames=c("Correct (True)", "Underestimation", "Overestimation","OLS")
# Create a long format data frame for correct mapping
plot_data <- data.frame(
  phenotype = rep(datainv$phenotype, 4),
  age = c(datainv$true.bioage, datainv$biom1, datainv$biom2,datainv$biomols),
  category = rep(modnames, each = nrow(datainv))
)

# Define colors for legend
colors <- c("Correct (True)" = "grey", "Underestimation" = "red", 
            "Overestimation" = "orange","OLS"="black")

# Plot
p1 <- ggplot() +
  # Add main scatter points
  geom_point(data = datainv, aes(y = cage, x = phenotype), size = 2) + ylim(c(0,80))+
  
  # Add legend-associated points
  geom_point(data = plot_data, aes(y = age, x = phenotype, 
                                   color = category, shape = category), 
             size = 4) +
  
  # Add identity reference lines with legend colors
  geom_abline(aes(intercept = 0, slope = 1, color = "Correct (True)"), size = 1) +
  geom_abline(aes(intercept = 0, slope = 1/2, color = "Underestimation"), size = 1) +
  geom_abline(aes(intercept = 0, slope = 1/0.5, color = "Overestimation"), size = 1) +
  geom_abline(aes(intercept = datainv$beta0mols, slope = datainv$beta1mols, color = "OLS"), size = 2) +
  
  # Set colors and shapes for legend
  scale_color_manual(values = colors, name = "Conditional Estimations") +
  scale_shape_manual(values = c(2, 3, 3,3), name = "Conditional Estimations") +

  # Labels and theme
  ylab("Chronological and Biological Age") +
  xlab("Phenotype") +
  ggtitle("Observed Chronological Age vs Phenotype") +
  theme_minimal()

# Print the plot
print(p1)


## ----fig.cap="Conditional likelihood of $\\boldsymbol\\theta$ with measurament errors on the phenotypes in the inverted approach with simple linear model.\\label{fig:invertedappclme}"----
# Define the conditional likelihood function with measurament errors
llprofileinv <- function(theta, y,x, beta0=0, beta1,sigma=1) 
  dnorm(theta, mean = beta0 + beta1 * y-x, sd = sigma, log = TRUE)

# Lets draw for each one of the three subject the corresponding 
# log likelihood profile for each of the three models
ss=seq(-50,50,length=10000)
llres=data.frame(NULL)
for(i in 1:3) for(j in 1:4){
  nll=llprofileinv(theta=ss,y = datainv$phenotype[i],x = datainv$cage[i],
                beta1 = datainv[i,grep("beta1",colnames(datainv))[j]],
                beta0 = datainv$beta0mols[1]*(i==4),
                sigma = datainv[i,grep("sigmam",colnames(datainv))[j]])
  nll=nll-max(nll)
  llres=rbind(llres,data.frame(theta=ss,loglik=nll,model=paste0("Model",j-1),
                               subject=paste0("Subject",i)))
}  
                                        
llres %>% 
  mutate(model=factor(model,labels=modnames),subject=factor(subject)) -> llres

# Create a separate dataset for vlines with unique bioage per subject
vline_data <- data.frame(
  subject = factor(paste0("Subject", 1:3)), # Ensure factors match facet levels
  xintercept = theta
)

# Plot with vertical lines appearing only once per subject
ggplot(data = llres, aes(x = theta, y = exp(loglik), color = factor(model))) +
  geom_line() +
  facet_grid(subject ~ ., scales = "free_y") +  # Allow different y-scales per subject
  theme_minimal() +

  # Add vlines using separate dataset
  geom_vline(data = vline_data, aes(xintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +

  # Define color legend
  scale_color_manual(values = colors, name = "Conditional Estimation") +

  # Labels and aesthetics
  xlab(TeX("$\\theta$")) + 
  ylab("Normalized conditional likelihood") + 
  ggtitle("Conditional likelihoods for the three Subjects") +
  
  # Remove duplicate y-axis labels and improve facet formatting
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1), # Align facet labels
    axis.title.y = element_blank(),                   # Remove global y-axis title
    strip.placement = "outside",                      # Move facet labels outside
    panel.spacing = unit(1, "lines"),                 # Space panels properly
    axis.text.y = element_blank(),                    # Remove individual y-axis labels
    axis.ticks.y = element_blank()                    # Remove y-axis ticks
  )


## ----fig.cap="Conditional likelihood of $\\boldsymbol\\theta$ with no measurament errors on the phenotypes in the inverted approach with simple linear model.\\label{fig:invertedappcl}"----
# Define the conditional likelihood function with no measurament errors
llprofileinv <- function(theta, y,x, beta0=0, beta1,sigma=1) 
  dnorm(theta, mean = beta0 + beta1 * y-x, sd = sqrt(sigma^2-mean(datainv$phenotype^2)), log = TRUE)

# Lets draw for each one of the three subject the corresponding 
# log likelihood profile for each of the three models
ss=seq(-50,50,length=10000)
llres=data.frame(NULL)
for(i in 1:3) for(j in 1:4){
  nll=llprofileinv(theta=ss,y = datainv$phenotype[i],x = datainv$cage[i],
                beta1 = datainv[i,grep("beta1",colnames(datainv))[j]],
                beta0 = datainv$beta0mols[1]*(i==4),
                sigma = datainv[i,grep("sigmam",colnames(datainv))[j]])
  nll=nll-max(nll)
  llres=rbind(llres,data.frame(theta=ss,loglik=nll,model=paste0("Model",j-1),
                               subject=paste0("Subject",i)))
}  
                                        
llres %>% 
  mutate(model=factor(model,labels=modnames),subject=factor(subject)) -> llres

# Create a separate dataset for vlines with unique bioage per subject
vline_data <- data.frame(
  subject = factor(paste0("Subject", 1:3)), # Ensure factors match facet levels
  xintercept = theta
)

# Plot with vertical lines appearing only once per subject
ggplot(data = llres, aes(x = theta, y = exp(loglik), color = factor(model))) +
  geom_line() +
  facet_grid(subject ~ ., scales = "free_y") +  # Allow different y-scales per subject
  theme_minimal() +

  # Add vlines using separate dataset
  geom_vline(data = vline_data, aes(xintercept = xintercept), 
             linetype = 2, linewidth = 1, color = "grey") +

  # Define color legend
  scale_color_manual(values = colors, name = "Conditional Estimation") +

  # Labels and aesthetics
  xlab(TeX("$\\theta$")) + 
  ylab("Normalized conditional likelihood") + 
  ggtitle("Conditional likelihoods for the three Subjects (no measurament errors)") +
  
  # Remove duplicate y-axis labels and improve facet formatting
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1), # Align facet labels
    axis.title.y = element_blank(),                   # Remove global y-axis title
    strip.placement = "outside",                      # Move facet labels outside
    panel.spacing = unit(1, "lines"),                 # Space panels properly
    axis.text.y = element_blank(),                    # Remove individual y-axis labels
    axis.ticks.y = element_blank()                    # Remove y-axis ticks
  )


## ----fig.cap="Root Mean Squared error (vertical axis) in point estimating the BA under the approach here proposed and the inverted approach with different sample sizes (horizonthal axis) and a different informative phenotypes (values of $\\beta_{1j}$). The prior standard deviation of the BA is the dashed horizontal segment, $\\sigma_\\theta=7.5$.\\label{fig:invertedappsimulation}"----
if (file.exists("simstudycomp.RData")) {
  load("simstudycomp.RData")
} else {
  set.seed(17)
  nrep=1000
  simres=expand.grid(beta1=c(0.1,1,2),n=c(20,50,100,150,200),sigmaY=c(1/4,1,4)) %>% 
    mutate(Proposed=NA,Inverted=NA)
  
  for(i in 1:nrow(simres)){ 
    mse.hatBA=mse.hatBAprime=0
    for(j in 1:nrep){
      cs=runif(simres$n[i],20,50)
      BA=cs+rnorm(simres$n[i],0,sd=7.5)
      ys=BA*simres$beta1[i]+rnorm(simres$n[i],0,sd=simres$sigmaY[i])
      olsmy=lm(ys~cs)
      olsinv=lm(cs~ys)
      hatBA=(ys-olsmy$coefficients[1])/olsmy$coefficients[2]
      hatBAprime=olsinv$fitted.values
      mse.hatBA=mse.hatBA+mean((BA-hatBA)^2)
      mse.hatBAprime=mse.hatBAprime+mean((BA-hatBAprime)^2)
    }
    simres$Proposed[i]=round(sqrt(mse.hatBA/nrep),0)
    simres$Inverted[i]=round(sqrt(mse.hatBAprime/nrep),0)
  }
  save(simres,file="simstudycomp.RData")
}

rr=simres %>% reshape2::melt(.,id.var=c("beta1", "n", "sigmaY"),value.name="rMSE",variable.name="estimate")

# Convert sigmaY and beta1 to factors with expression labels
rr$sigmaY <- factor(rr$sigmaY, 
  levels = c("0.25", "1", "4"),
  labels = c("sigma[Y] == 1/4", "sigma[Y] == 1", "sigma[Y] == 4")
)
rr$beta1 <- factor(rr$beta1, 
  levels = c("0.1", "1", "2"),
  labels = c("beta[1] == 0.1", "beta[1] == 1", "beta[1] == 2")
)

# Define colors for legend
colors <- c("Proposed" = "red", "Inverted" = "blue")

ggplot(data=rr,aes(x=n,y=rMSE,color=estimate)) + geom_point()+geom_line()+
  facet_grid(sigmaY ~ beta1, labeller = label_parsed) +
  ylim(c(0,30))+
  geom_hline(yintercept=7.5,linetype=2)+ylab("Average rMSE (in years)")+xlab("Sample size (n)")+theme_minimal()+ggtitle("Root Mean Squared Error for both approaches")+scale_color_manual(values = colors, name = "Method")


## ----fig.cap="Trace of all parameters in the toy example (one phenotype) \\label{fig-toyextracechainuni}"----

print(p.tracechainuni)


## ----fig.cap="Posterior approximation of all parameters in the toy example (two phenotypes) \\label{fig-toyexpostmulti}"----
print(p.postmulti)


## ----fig.cap="Trace of all parameters in the toy example (two phenotypes) \\label{fig-toyextracechainmulti}"----
print(p.tracechainmulti)


## ----fig.cap="Posterior distribution of $\\sigma_{\\epsilon{ij}}$ for the PC analysis on the USAF 1967 data \\label{fig-usaf1967postsigma}"----
print(p.postsigmaex1)


## ----fig.cap="Approximation of posterior distribution of the correlation among the phenotypes in the DL model for the USAF 1967 data set. \\label{fig-usaf1967postSigmadl}"----

print(p.postcorrusaf1967dl)



## ---------------------------------------------------------------

# Define the FFNN structure graph using grViz with horizontal layout
graph <- grViz("
digraph ffnn {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle, style = filled, fillcolor = lightblue, fontname = Helvetica, fontsize = 10]
  
  Input [label = 'Input Layer\\n(n = # of covariates)']
  Dense1 [label = 'Dense Layer 1\\n3 * JN units\\n(ReLU)']
  Dense2 [label = 'Dense Layer 2\\n3 * JN units\\n(ReLU)']
  Dropout [label = 'Dropout Layer\\n(rate = 0.3)']
  Dense3 [label = 'Dense Layer 3\\n3 * JN units\\n(ReLU)']
  
  subgraph cluster_outputs {
    label = 'Output Branches';
    style = dashed;
    node [shape = rectangle, style = filled, fillcolor = lightyellow]
    Output1 [label = 'Output 1\\n(Quantitative)']
    Output2 [label = 'Output 2\\n(Categorical)']
  }
  
  Input -> Dense1;
  Dense1 -> Dense2;
  Dense2 -> Dropout;
  Dropout -> Dense3;
  Dense3 -> Output1;
  Dense3 -> Output2;
}
")

# Convert the DiagrammeR graph to SVG
svg <- DiagrammeRsvg::export_svg(graph)

# Write the SVG output to a file
svg_file <- "ffnn_diagram_horizontal.svg"
writeLines(svg, svg_file)

# Convert the SVG file to PNG format for inclusion in LaTeX
png_file <- "ffnn_diagram.png"
rsvg::rsvg_png(svg_file, png_file)



## ----fig.cap="Structure of the FFNN model for the multivariate phenotype \\label{fig-ffnn}",out.width='100%'----
knitr::include_graphics("ffnn_diagram.png")


## ----results='asis'---------------------------------------------
# Assumendo che 'datamutuaused' sia il tuo dataset originale
sumdata_quant <- datamutuaused %>% 
  select(where(is.numeric)) %>% 
  summarise_all(list(mean = ~mean(., na.rm = TRUE), 
                     sd = ~sd(., na.rm = TRUE),
                     min = ~min(., na.rm = TRUE),
                     max = ~max(., na.rm = TRUE))) %>%
  pivot_longer(everything(), 
               names_to = c("Variable", ".value"), 
               names_pattern = "(.*)_(.*)")

sumdata_cat <- datamutuaused %>% 
  select(where(~!is.numeric(.))) %>% 
  summarise_all(~paste(names(sort(table(.), decreasing = TRUE))[1], 
                       "(", round(max(prop.table(table(.))) * 100, 1), "%)", sep="")) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Mode (Freq%)")

# Unione delle sintesi numeriche e categoriche
sumdataset <- bind_rows(
  sumdata_quant %>% 
    mutate(`Summary` = paste0("Mean=", round(mean,2), ", SD=", round(sd,2),
                              ", Min=", round(min,2), ", Max=", round(max,2))) %>% 
    select(Variable, Summary),
  sumdata_cat %>%
    rename(Summary = `Mode (Freq%)`)
)

kable(sumdataset, format = "latex", escape = TRUE,
      caption = "\\label{tab:smadatasummary} Summary of SMA data set used in the analysis ($n=911$ individuals).",
      row.names = FALSE, align = "c", booktabs = TRUE) %>%
  kable_styling(latex_options = c("HOLD_position"))


## ----results='asis'---------------------------------------------
# Creo un data frame e assegno i nomi delle righe
table_data <- data.frame(fittstat)

# Genero la tabella in LaTeX tramite kable e kableExtra
kable(table_data, format = "latex", escape = FALSE,
      caption = "\\label{tab:mutuar2} The $R^2$ values for the fitted FFNN on the USAF 1967 data set.",
      align = "c", booktabs = TRUE, col.names = "$R^2$") %>%
  kable_styling(latex_options = c("HOLD_position"))



## ----fig.cap="Approximation of posterior distribution of the correlation among the phenotypes in the DL model for SDA data set. \\label{fig-mutuapostSigmadl}"----

print(p.postcorrmutua)



## ----echo=FALSE, include=FALSE,eval=FALSE-----------------------
# knitr::purl(input = "bayes-enfoque-ejemplo-v6.Rmd", output = "bayes-bioage.R", documentation = 1)

