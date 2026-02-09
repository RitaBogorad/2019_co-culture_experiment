# packages
library(ggplot2)
library(nls2)

#script for analyses of the fatty acid 13C data for three selected fatty acids. 
#the dataset is loaded in its entirety but each fatty acid is analyses separately
#based on the model comparisons: model 2 (fit2) is selected as best model and used to compare max uptake for each FA

#output of fit 2 at the end of the script - order of levels for a intercept: AX DPS DVM

# References
#Webb WL, Newton M, Starr D (1974) Carbon dioxyde exchange of Alnus rubra: a mathematical model. Oecologia 17:281-291
#->saturating exponential = specific case of logistic regression
# Exponential Decay (increasing form) - https://people.richland.edu/james/lecture/m116/logs/models.html
# non linear regression
setwd("~/phd/metatransccr exp")
data<-read.table("all_known_FAs_13C_and_rel_abun.txt",h=TRUE)
data$time.code<-c("20.01.2022 12:50","20.01.2022 16:30","20.01.2022 19:30","21.01.2022 07:00")[factor(data$Time)]
data$time.code<-as.POSIXct(strptime(data$time.code,format='%d.%m.%Y %H:%M'))
data$time.h<-c(0,3.67,6.67,18.17)[factor(data$Time)]


# 14:0####
dac14<-data[data$Treat=="DPS"&data$Peak=="C14:0"|
              data$Treat=="DVM"&data$Peak=="C14:0"|
              data$Treat=="AX"&data$Peak=="C14:0",]

my.dac14<-data.frame(X=dac14$time.h,Y=dac14$AT.13C.12C,Treat=dac14$Treat)
my.dac14$Treat=as.factor(my.dac14$Treat)
my.data=my.dac14

##common analysis block


# Model 1: Y ~ a * (1 - exp(-X/b))
# Model 2: Y ~ a[Treat] * (1 - exp(-X/b))
# Model 3: Y ~ a * (1 - exp(-X/b[Treat]))
# Model 4: Y ~ a[Treat] * (1 - exp(-X/b[Treat]))
# Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
# 1     27     85.466                                
# 2     25     34.657  2 50.809 18.3253 1.259e-05 ***
#   3     25     61.390  0  0.000                      
# 4     23     33.388  2 28.002  9.6448 0.0009082 ***

#anova(fit2,fit4) #not sign different


#fit3
#Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a1   7.6747     0.4901  15.658 1.98e-14 ***
#   a2   7.7734     0.4894  15.885 1.43e-14 ***
#   a3  11.2987     0.5846  19.327  < 2e-16 ***

# 16:1w7####
dac16<-data[data$Treat=="DPS"&data$Peak=="C16:1w7"|
              data$Treat=="DVM"&data$Peak=="C16:1w7"|
              data$Treat=="AX"&data$Peak=="C16:1w7",]

dac16.2<-data.frame(X=dac16$time.h,Y=dac16$AT.13C.12C,Treat=dac16$Treat)
dac16.2$Treat=as.factor(dac16.2$Treat)
my.data=dac16.2

# Model 1: Y ~ a * (1 - exp(-X/b))
# Model 2: Y ~ a[Treat] * (1 - exp(-X/b))
# Model 3: Y ~ a * (1 - exp(-X/b[Treat]))
# Model 4: Y ~ a[Treat] * (1 - exp(-X/b[Treat]))
# Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
# 1     36     283.14                                
# 2     34     183.17  2 99.975  9.2789 0.0006086 ***
#   3     34     218.25  0  0.000                      
# 4     32     175.51  2 42.733  3.8956 0.0306053 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > anova(fit2,fit4)
# Analysis of Variance Table
# 
# Model 1: Y ~ a[Treat] * (1 - exp(-X/b))
# Model 2: Y ~ a[Treat] * (1 - exp(-X/b[Treat]))
# Res.Df Res.Sum Sq Df Sum Sq F value Pr(>F)
# 1     34     183.17                         
# 2     32     175.51  2  7.653  0.6977 0.5052
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a1  13.0936     0.8813  14.858  < 2e-16 ***
#   a2  10.1614     0.8048  12.626 2.18e-14 ***
#   a3  14.5129     0.8780  16.529  < 2e-16 ***
#   b    2.2533     0.4797   4.697 4.23e-05 ***

# 17:0 ####   
#NOTE: no T0 values available for this FA???

dac17<-data[data$Treat=="DPS"&data$Peak=="C17:0"|
              data$Treat=="DVM"&data$Peak=="C17:0",]
dac17.2<-data.frame(X=dac17$time.h,Y=dac17$AT.13C.12C,Treat=dac17$Treat)
dac17.3<-data.frame(X=c(0,0),Y=c(0,0),Treat=c("DPS","DVM"))
dac17.2=rbind(dac17.2,dac17.3)
dac17.2$Treat=as.factor(dac17.2$Treat)
my.data=dac17.2

# Analysis of Variance Table
# 
# Model 1: Y ~ a * (1 - exp(-X/b))
# Model 2: Y ~ a[Treat] * (1 - exp(-X/b))
# Model 3: Y ~ a * (1 - exp(-X/b[Treat]))
# Model 4: Y ~ a[Treat] * (1 - exp(-X/b[Treat]))
# Res.Df Res.Sum Sq Df  Sum Sq F value    Pr(>F)    
# 1     18    182.079                                 
# 2     17     68.582  1 113.497  28.133 5.829e-05 ***
#   3     17    121.953  0   0.000                      
# 4     16     68.130  1  53.823  12.640  0.002637 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > anova(fit2,fit4) 
# Analysis of Variance Table
# 
# Model 1: Y ~ a[Treat] * (1 - exp(-X/b))
# Model 2: Y ~ a[Treat] * (1 - exp(-X/b[Treat]))
# Res.Df Res.Sum Sq Df  Sum Sq F value Pr(>F)
# 1     17     68.582                          
# 2     16     68.130  1 0.45234  0.1062 0.7487

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a1  12.2171     0.8036  15.202 2.50e-11 ***
#   a2  17.7476     0.9723  18.253 1.33e-12 ***
#   b    2.2969     0.4681   4.906 0.000133 ***

#common analysis block####

#substrat t0 measurements from all measurements (functions fits through the origin so we have substract our t0)
my.data_zerocor=data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, colnames(my.data))))
for (treat in levels(my.data$Treat)) {
  my.data_subset=my.data[my.data$Treat==treat,]
  my.data_subset$Y=my.data_subset$Y-as.numeric(my.data_subset[my.data_subset$X==0,]['Y'])
  my.data_zerocor=rbind(my.data_zerocor,my.data_subset)
}

#exploratorive plot
ggplot(data=my.data_zerocor, aes(x=X,y=Y,color=Treat))+geom_point()

#no distinction between treatments
fit1<-nls(Y~a*(1-exp(-X/b)),
          data=my.data_zerocor,
          start=list(a=15,b=2),
          nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1e-05,
                      printEval=F, warnOnly=F))

#distinction between treatments only for a 
fit2<-nls(Y~a[Treat]*(1-exp(-X/b)),
          data=my.data_zerocor,
          start=list(a=rep(10,nlevels(my.data_zerocor$Treat)),b=2),
          nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1e-05,
                      printEval=F, warnOnly=F))

#distinction between treatments only for b 
fit3<-nls(Y~a*(1-exp(-X/b[Treat])),
          data=my.data_zerocor,
          start=list(a=10,b=rep(2,nlevels(my.data_zerocor$Treat))),
          nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1e-05,
                      printEval=F, warnOnly=F))


#distinction between treatments for both parameters
fit4<-nls(Y~a[Treat]*(1-exp(-X/b[Treat])),
          data=my.data_zerocor,
          start=list(a=rep(10,nlevels(my.data_zerocor$Treat)),b=rep(2,nlevels(my.data_zerocor$Treat))),
          nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1e-05,
                      printEval=F, warnOnly=F))

levels(my.data_zerocor$Treat)
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)

anova(fit1,fit2,fit3,fit4) #fit 2 and fit4 are better 
anova(fit2,fit4) #fit2 and fit4 are not significantly different, therefore most simple model is preferred (fit2)
#retain the fit 2 and look at the summary p values


ggplot(data=my.data_zerocor,mapping=aes(x=X,y=Y,color=Treat))+geom_point()+
  geom_smooth(method = "nls",formula = y ~ a*(1-exp(-x/b)),se = FALSE,method.args = list(start = c(a = 10, b = 2)))+
  theme_minimal()

#making a ggplot for all data with the webb function fitted to the data
my.data=dac14
dac16.2
dac17.2

my.data_zerocor=data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, colnames(my.data))))
for (treat in levels(my.data$Treat)) {
  my.data_subset=my.data[my.data$Treat==treat,]
  my.data_subset$Y=my.data_subset$Y-as.numeric(my.data_subset[my.data_subset$X==0,]['Y'])
  my.data_zerocor=rbind(my.data_zerocor,my.data_subset)
}

#adjust code to fit each FA dataset indiv, otherwise duplicated data!
graph_data_14=my.data_zerocor;graph_data_14$FA='14'
graph_data_17=my.data_zerocor;graph_data_17$FA='17'
graph_data_16=my.data_zerocor;graph_data_16$FA='16'
graph_data = rbind(graph_data_14,graph_data_17,graph_data_16)

ggplot(data=graph_data,mapping=aes(x=X,y=Y,color=Treat))+geom_point()+
  geom_smooth(method = "nls",formula = y ~ a*(1-exp(-x/b)),se = FALSE,method.args = list(start = c(a = 10, b = 2)))+
  theme_light()+facet_grid(FA~.)

#order of levels for a intercept: AX DPS DVM

#fit 2 output for 14:0####
# Formula: Y ~ a[Treat] * (1 - exp(-X/b))
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a1   7.6747     0.4901  15.658 1.98e-14 ***
#   a2   7.7734     0.4894  15.885 1.43e-14 ***
#   a3  11.2987     0.5846  19.327  < 2e-16 ***
#   b    3.1129     0.4551   6.839 3.61e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.177 on 25 degrees of freedom
# 
# Number of iterations to convergence: 6 
# Achieved convergence tolerance: 2.246e-06


#fit 2 output for 16:1w7####
# Formula: Y ~ a[Treat] * (1 - exp(-X/b))
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a1   7.6747     0.4901  15.658 1.98e-14 ***
#   a2   7.7734     0.4894  15.885 1.43e-14 ***
#   a3  11.2987     0.5846  19.327  < 2e-16 ***
#   b    3.1129     0.4551   6.839 3.61e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.177 on 25 degrees of freedom
# 
# Number of iterations to convergence: 6 
# Achieved convergence tolerance: 2.246e-06

#fit 2 output for 17:0####
# Formula: Y ~ a[Treat] * (1 - exp(-X/b))
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a1   7.6747     0.4901  15.658 1.98e-14 ***
#   a2   7.7734     0.4894  15.885 1.43e-14 ***
#   a3  11.2987     0.5846  19.327  < 2e-16 ***
#   b    3.1129     0.4551   6.839 3.61e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.177 on 25 degrees of freedom
# 
# Number of iterations to convergence: 6 
# Achieved convergence tolerance: 2.246e-06