#We define the Polya Tree method for continuous data called PT

#Inputs:
#x: Observations of the random variable X
#y: Observations of the random variable Y
#level_max: The maximum number of levels to analyze.
#At level k, the a priori branching probability is Dirichlet(ck^(power_m), ck^(power_m), ck^(power_m), ck^(power_m)) 
#printInfo: To print the evolution of the log-Bayes factor at each level

#Output:
#The logged Bayes Factor for evidence of H0 vs H1 (independence vs dependence). 

PT = function (x, y, level_max = 52, c=1, power_m=2, printInfo=0)
{
  #We check that the number of observations of X and Y are the same
  if(NROW(x)!=NROW(y))
  {
    return("x and y need to be the same length")
  }
  
  #We check that there is at least one level.
  if(level_max<1)
  {
    return("At least 1 level is needed")
  }
  
  #If a quadrant has one or zero observations, then there is no evidence favouring dependence or independence. 
  if(NROW(x)<=1)
  {
    return(0)
  }
  
  #We initialize an empty string for each datapoint. We will later fill the string with the subquadrants that each datapoint is in.
  position_point = rep("", NROW(x))
  
  #We order our datapoints
  order_x = order(x)
  x = x[order_x]
  y = y[order_x]
  
  order_y = order(y)
  y = y[order_y]
  x = x[order_y]
  
  res_x = x
  res_y = y
  
  logb_numerator = 0
  logb_denominator = 0
  

  Indicator2 = rep(FALSE, NROW(x)) #This tells us if the quadrant the datapoint belongs to has just 1 point
  
  for(i in 1:level_max)
  {
    unique_prev_quadrants = unique(position_point[!Indicator2])
    
    if(length(unique_prev_quadrants)==0)
    {
      return(logb)
    }
    
    #We see if the data goes left or right in the i'th level
    position_x = ifelse((2*res_x)>=1, 1, 0)
    res_x = ifelse((2*res_x)>=2, 1, (2*res_x) %% 1) #The first condition is necessary for the case that you have an observation exactly equal to 1.
    
    #We see if the data goes up or down in the i'th level
    position_y = ifelse((2*res_y)>=1, 1, 0)
    res_y = ifelse((2*res_y)>=2, 1, (2*res_y) %% 1) #The first condition is necessary for the case that you have an observation exactly equal to 1.   
    
    alpha = c*i^power_m
    
    
    #For each quadrant we see how many points there are in the 0, 1, 2, 3 partition
    for(j in 1: length(unique_prev_quadrants))
    {
      Indicator = position_point %in% unique_prev_quadrants[j]
      
      subquadrant = (position_x + 2*position_y)[Indicator]
      
      n0 = sum(subquadrant==0)
      n1 = sum(subquadrant==1)
      n2 = sum(subquadrant==2)
      n3 = sum(subquadrant==3)
      
      if((n0+n1+n2+n3)>1)
      {
        logb_numerator = logb_numerator + lbeta(n0+n2+2*alpha, n1+n3+2*alpha)+lbeta(n0+n1+2*alpha, n2+n3+2*alpha) - lbeta(2*alpha, 2*alpha)-lbeta(2*alpha, 2*alpha)
        logb_denominator = logb_denominator + lgamma(n0+alpha)+lgamma(n1+alpha)+lgamma(n2+ alpha)+ lgamma(n3+alpha)+lgamma(4*alpha)-lgamma(n0+n1+n2+n3+4*alpha) - 4* lgamma(alpha)
        
      }
      
      else
      {
        Indicator2 = (Indicator | Indicator2) #Important note: We use | instead of || to do elementwise comparison
      }
      
      logb = logb_numerator - logb_denominator
      
    }
    
    
    position_point = paste0(position_point, (position_x + 2*position_y)) 
    
    
    if(printInfo>0)
    {
      print(c(i,logb, logb_numerator, logb_denominator))  #You can print this to see how many levels the tree has, and the evidence for H0 at each level...
      
    }
    
  }
  
  print(paste0("Maximum level reached. ", length(unique(position_point[Indicator2])), " points uniquely assigned"))
  
  return(logb)
}

###############################################################################################
###############################################################################################
#We define the Polya Tree method for data with zeros

EnhancedPT = function (x, y, level_max = 52, c=1, power_m=2, printInfo=0)
{
  if(NROW(x)!=NROW(y))
  {
    return("x and y need to be the same length")
  }

  if(level_max<1)
  {
    return("At least 1 level is needed")
  }

  if(NROW(x)==1)
  {
    print("Preferably at least 2 datapoints are needed")
    return(0)
  }

  Indicator_00 = (x==0)&(y==0)
  Indicator_01 = (x==0)&(y!=0)
  Indicator_10 = (x!=0)&(y==0)
  Indicator_11 = (x!=0)&(y!=0)

  n0 = sum(Indicator_00)
  n1 = sum(Indicator_01)
  n2 = sum(Indicator_10)
  n3 = sum(Indicator_11)

  if(n3==NROW(x))
  {
    logb = PT(x[Indicator_11], y[Indicator_11], level_max, c, power_m, printInfo)
  }

  else
  {
    alpha = c
    logb_numerator = lbeta(n0+n2+2*alpha, n1+n3+2*alpha)+lbeta(n0+n1+2*alpha, n2+n3+2*alpha) - lbeta(2*alpha, 2*alpha)-lbeta(2*alpha, 2*alpha)
    logb_denominator = lgamma(n0+alpha)+lgamma(n1+alpha)+lgamma(n2+ alpha)+ lgamma(n3+alpha)+lgamma(4*alpha)-lgamma(n0+n1+n2+n3+4*alpha) - 4* lgamma(alpha)
    logb = logb_numerator - logb_denominator

    logb = logb + PT(x[Indicator_11], y[Indicator_11], level_max, c, power_m, printInfo)
  }

  return(logb)
}

