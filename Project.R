rm(list=ls()) 


#REmember to import package, use library(deSolve)

# load required packages 

install.packages("deSolve")
install.packages("scatterplot3d")
library("deSolve")
require(deSolve)
library("scatterplot3d")
require(scatterplot3d)



#Model Equations: 

# also can hard code this into Model_LVE
#f <- function(params, x){
#  
#}

#g <- function(params, y){
#  
#}



RunModel<-function(params){
  Model_LVE <- function(t, variables, parameters) 
  {with( as.list(c(variables,parameters)),  #lets us access variables and parameters by name 
         { 
           #System of ODES:
           
           dx = b*x* (( 1 - (x)/(min(k, (p-theta_y*y-theta_z*z-theta_w*w)/q)))) - ((m_y*x)/(x+a_y))*y
           
           dy = e_y * min(1,((p-theta_y*y-theta_z*z-theta_w*w)/x)/theta_y) * 
                         ((m_y*x)/(x+a_y))*y - ((m_z*y)/(y+a_z))*z - ((m_w1*y)/(y+a_w1))*w-delta_y*y
           
           dz = e_z * min(1,theta_y/theta_z)*((m_z*y)/(y+a_z))*z - ((m_w2*y)/(z+a_w2))*w - delta_z*z
           
           dw = e_w1 * min(1,theta_y/theta_w) * ((m_w1*y)/(y+a_w1))*w + e_w2* min(1,theta_z/theta_w) * ((m_w2*z)/(z+a_w2))*w-delta_w*w
           
           #End of ODEs
           return(list(c(dx,dy,dz,dw)))
         }) 
  } #End Model Equations


#Initialization: load parameter values:

b = params[1]
k = params[2]
p = params[3]
q = params[4]

delta_y = params[5]
delta_z = params[6]
delta_w = params[7]

theta_y = params[8]
theta_z = params[9]
theta_w = params[10]

e_y = params[11]
e_z = params[12]
e_w1 = params[13]
e_w2 = params[14]

m_y = params[15]
m_z = params[16]
m_w1 = params[17]
m_w2 = params[18]

a_y = params[19]
a_z = params[20]
a_w1 = params[21]
a_w2 = params[22]



#set parameter list:
parameter_values<-c(b=b, k=k, p=p, e_y=e_y, e_z=e_z, delta_y = delta_y, delta_z=delta_z, q=q, theta_y=theta_y, 
                    theta_z=theta_z, e_w1=e_w1, e_w2=e_w2, m_y=m_y, m_z=m_z, m_w1=m_w1, m_w2=m_w2, a_y=a_y,
                    a_z=a_z, a_w1=a_w1, a_w2=a_w2, delta_w=delta_w, theta_w=theta_w)
  
  
tend=200 #Days for simulation
time_values<-seq(0,tend)

#Initialization: Set initial conditions:  
initial_values <- c(
  x=.5,
  y=.25,
  z=.25, 
  w=.25
)

#Solve Model 
Model_solutions <- ode(
  y = initial_values,
  times = time_values,
  func = Model_LVE,
  parms = parameter_values 
)

Model_solutions<-as.data.frame(Model_solutions) 
print(Model_solutions)





plot(time_values,Model_solutions$x,type="l",col=2,lwd=2,ylim=c(0,max(.5,k)),
     main=paste("Plot for K=",k),ylab="Population Density",xlab="Time", lty=3)
lines(time_values,Model_solutions$y,type="l",col=4,lwd=2, lty=2)
lines(time_values,Model_solutions$z,type="l",col=1,lwd=2)
lines(time_values,Model_solutions$w,type="l",col=3,lwd=2)
legend(x = "topright", legend=c("x", "y", "z", "w"), 
       fill = c("red","blue", "black", "green"))

if (require(scatterplot3d))
  scatterplot3d(Model_solutions[,-1], type = "l")
print(Model_solutions[,2:4])

}


# order: b,k,p,q
# delta_y, delta_z, delta_w
# theta_y, theta_z, theta_w, 
# e_y, e_z, e_w1, e_w2, 
# m_y, m_z, m_w1, m_w2
# a_y, a_z, a_w1, a_w2 

# .7-1 for e_w
# m_w1 0.8, m_w2 0.01
# a_w1   0.5, a_w2  0.5, 
# theta_w between 0 and less than other theta (0.006)
# theta_w 0.004
# delta_w needs to be small 

P <- c(1.2, 1, 0.12, 0.0038,          # b,k,p,q
       0.25, 0.003, 0.001,            # delta_y, delta_z, delta_w
       0.03, 0.013, 0.004,            # theta_y, theta_z, theta_w
       0.8, 0.75, 0.7, 0.7,           # e_y, e_z, e_w1, e_w2
       0.81, 0.03, 0.8, 0.01,         # m_y, m_z, m_w1, m_w2
       0.25, 0.75, 0.5, 0.5)          # a_y, a_z, a_w1, a_w2


ks <- c(0.1, 0.18, 0.4, 0.6, 1, 4, 8, 10)
      

densities <- c()

for (k in ks){
  P <- c(1.2, k, 0.12, 0.0038,          # b,k,p,q
       0.25, 0.003, 0.001,            # delta_y, delta_z, delta_w
       0.03, 0.013, 0.004,            # theta_y, theta_z, theta_w
       0.8, 0.75, 0.5, 0.01,           # e_y, e_z, e_w1, e_w2
       0.81, 0.03, 0.0001, 0.01,         # m_y, m_z, m_w1, m_w2
       0.25, 0.75, 0.5, 0.5)          # a_y, a_z, a_w1, a_w2
  RunModel(P)
  #model_solutions <- RunModel(P)
  #cur_dens <- model_solutions[,-1]
  #densities <- c(densities, cur_dens)
}

  
