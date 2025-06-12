rm(list=ls()) 


#REmember to import package, use library(deSolve)

# load required packages 

install.packages("deSolve")
install.packages("scatterplot3d")
library("deSolve")
require(deSolve)
library("scatterplot3d")
require(scatterplot3d)



RunModel<-function(params){
  Model_LVE <- function(t, variables, parameters) 
  {with( as.list(c(variables,parameters)),  #lets us access variables and parameters by name 
         { 
           #System of ODES:
           
           dx = b*x* (( 1 - (x)/(min(k, (p-theta_y*y-theta_z*z)/q)))) - ((m_y*x)/(x+a_y))*y
           
           dy = e_y * min(1,((p-theta_y*y-theta_z*z)/x)/theta_y) * 
             ((m_y*x)/(x+a_y))*y - ((m_z*y)/(y+a_z))*z -delta_y*y
           
           dz = e_z * min(1,theta_y/theta_z) * ((m_z*y)/(y+a_z))*z - delta_z*z
           
           #End of ODEs
           return(list(c(dx,dy,dz)))
         }) 
  } #End Model Equations
  
  
  #Initialization: load parameter values:
  
  b = params[1]
  k = params[2]
  p = params[3]
  q = params[4]
  
  delta_y = params[5]
  delta_z = params[6]
  
  theta_y = params[7]
  theta_z = params[8]
  
  e_y = params[9]
  e_z = params[10]
  
  m_y = params[11]
  m_z = params[12]
  
  a_y = params[13]
  a_z = params[14]
  
  
  
  #set parameter list:
  parameter_values<-c(b=b, k=k, p=p, e_y=e_y, e_z=e_z, delta_y = delta_y, delta_z=delta_z, 
                      q=q, theta_y=theta_y, theta_z=theta_z, m_y=m_y, m_z=m_z, a_y=a_y, a_z=a_z)
  
  
  tend=2000 #Days for simulation
  time_values<-seq(0,tend)
  
  #Initialization: Set initial conditions:  
  initial_values <- c(
    x=.5,
    y=.25,
    z=.25
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
  par(mar=c(5,6,4,1)+.1)
  plot(time_values,Model_solutions$x,type="l",col=2,lwd=2,ylim=c(0,max(.5,k+10)),
       main=paste("Plot for K=",k),ylab="Densities (mgC/L)",xlab="Days", lty = 3)
  lines(time_values,Model_solutions$y,type="l",col=4,lwd=2,lty=2 )
  lines(time_values,Model_solutions$z,type="l",col=1,lwd=2)
  
  legend(x = "topright", legend=c("x", "y", "z"), 
         fill = c("red","blue", "black"))
  
  #plot(Model_solutions)
  if (require(scatterplot3d))
    scatterplot3d(Model_solutions[,-1], type = "l") #, col = 4)
    print(Model_solutions)
  #plane(Model_solutions)
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

ks <- c(0.1, 0.18, 0.4, 0.6, 1, 4, 8, 10)


P <- c(1.2, 1, 0.12, 0.0038,          # b,k,p,q
       0.25, 0.003,            # delta_y, delta_z
       0.03, 0.013,            # theta_y, theta_z
       0.8, 0.75,           # e_y, e_z
       0.81, 0.03,          # m_y, m_z, 
       0.25, 0.75)          # a_y, a_z, a_w1, a_w2


for (k in ks){
  P <- c(1.2, k, 0.12, 0.0038,          # b,k,p,q
         0.25, 0.003,            # delta_y, delta_z
         0.03, 0.013,            # theta_y, theta_z
         0.8, 0.75,           # e_y, e_z
         0.81, 0.03,          # m_y, m_z, 
         0.25, 0.75) 
  RunModel(P)
}



#RunModel(P)
