# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:26:20 2018

@author: alian
"""

# Pulled from GitHub: Original file 2D-CVM-perturb-expt-1-2-2018-01-07.py 

####################################################################################################
# Alianna J. Maren
# Computing configuration variables for the Cluster Variation Method
####################################################################################################
# Import the following Python packages

import itertools
import numpy as np
import pylab
import matplotlib
from math import exp
from math import log
from matplotlib import pyplot as plt


####################################################################################################
####################################################################################################
#
# Detailed code documentation is JUST ABOVE main(), at the very end of this program. 
#
####################################################################################################
####################################################################################################
#
#
# The crucial equations are as follows (taken from AJM's 2014 paper, "The Cluster Variation Method II: 
#   2-D Grid of Zigzag Chains":
#   h = exp(beta*epsilon/4) & lambda = 0 (Beginning of Appendix B, replicating Eqn. 2-16.)
#   We can set beta = Boltzmann's constant = 1. 
#   Thus, eps1 = epsilon = 4*log(h)
#   For the equilibrium case (which is where we have an analytic solution), eps0 = 0. 
# Thus, x1 = x2 = 0; h controls the distribution among the z, w, & y values.
# At equilibrium, when eps1 = 0, z1 = z6 = z3 = z4 = 0.125; z2 = z5 = 0.125, however,
#  they are "degenerate" so their total values each = 0.25. 
#


####################################################################################################
####################################################################################################
#
# Procedure to welcome the user and identify the code
#
####################################################################################################
####################################################################################################

def welcome ():


    print()
    print(" ******************************************************************************")
    print()
    print(" Welcome to the 2-D Cluster Variation Method")
    print("   Equilibrium (analytic) calculations")    
    print("   Version 1.2, 10/29/2018, A.J. Maren")
    print(" This version computes the negative entropies for the case where ")
    print("  there is no interaction enthalpy.")
    print() 
    print(" This means that the values for the configuration variables are")
    print("   all probabilistically determined.")
    print() 
    print(" For comments, questions, or bug-fixes, contact: alianna.maren@northwestern.edu")
    print(" Alternate email address: alianna@aliannajmaren.com")
    print()
    print("   NOTE: In these calculations, x1 = A (units are at value 1),")
    print("                            and x2 = B (units are at value 0).")
    print()
    print(" ******************************************************************************")
    print()
    return()

                                                
                                                                                                
####################################################################################################
####################################################################################################
#
# Code Documentation: 
# The MAIN module comprises of calls to:
#  (1) Welcome
#  (2) Compute the values for the xArray, which will define the granularity of the work    
#  (3) Compute the negative entropies of various sorts
#  (4) Compute the activation and interactive enthalpies, and the free energy
#  (5) Print and plot variuos results.     
#
#    NOTE: The user can modify the activation and interactive enthalpy terms, eps0 and eps1, respectively,
#          via directly changing values in __main__.  
#
#
    
####################################################################################################
####################################################################################################
#
# Function to create the x-values. 
#
####################################################################################################
####################################################################################################
    
def createXValues(xArray, xTotalSteps, xStep, xIncr):
    newX = 0.0
    
    for j in range (0,xTotalSteps, xStep):
        newX = newX + xIncr
        xArray[j] = newX

    
    return (xArray)
    
####################################################################################################
####################################################################################################
#
# Function to create the x-entrlpy. 
#
####################################################################################################
####################################################################################################
    
def createNegXEntropyValues(xArray, negXEntropyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        x = xArray[j]
        negXEntropyArray[j] = x*log(x) + (1-x)*log(1-x)
    
    return (negXEntropyArray)
    
     
####################################################################################################
####################################################################################################
#
# Function to create the y-entropy. 
#
####################################################################################################
####################################################################################################
    
def createNegYEntropyValues(xArray, negYEntropyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        x = xArray[j]
        y1=x*x
        y2=x*(1.-x)
        y3=(1.-x)*(1.-x)
        negYEntropyArray[j] = y1*log(y1) + 2.*(y2)*log(y2) + y3*log(y3)
    
    return (negYEntropyArray)
    
####################################################################################################
####################################################################################################
#
# Function to create the z-entropy. 
#
####################################################################################################
####################################################################################################
    
def createNegZEntropyValues(xArray, negZEntropyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        x = xArray[j]
        q = 1.0-x
        z1=x*x*x
        z2=x*x*q
        z3=x*q*x
        z4=q*x*q
        z5=q*q*x
        z6=q*q*q
        negZEntropyArray[j] = z1*log(z1) + 2.*(z2)*log(z2) + z3*log(z3) + z4*log(z4) + 2.*z5*log(z5) + z6*log(z6)
    
    return (negZEntropyArray)
    

####################################################################################################
####################################################################################################
#
# Function to create the numerator of the total entropy. 
#
####################################################################################################
####################################################################################################
    
def createNegYWEntropyValues(xArray, negYEntropyArray, negYWEntropyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        negY = negYEntropyArray[j]
        negW= negYEntropyArray[j]        
        negYWEntropyArray[j] = 2*negY + negW
    
    return (negYWEntropyArray)     



####################################################################################################
####################################################################################################
#
# Function to create the denominator of the total entropy. 
#
####################################################################################################
####################################################################################################
    
def createNegXZEntropyValues(xArray, negXEntropyArray, negZEntropyArray, negXZEntropyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        negX = negXEntropyArray[j]
        negZ = negZEntropyArray[j]        
        negXZEntropyArray[j] = 2*negZ + negX
    
    return (negXZEntropyArray)     
    

####################################################################################################
####################################################################################################
#
# Function to create the denominator of the total entropy. 
#
####################################################################################################
####################################################################################################
    
def createNegTotEntropyValues(xArray, negYWEntropyArray, negXZEntropyArray, negTotEntropyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        negYWEntropy = negYWEntropyArray[j]
        negXZEntropy = negXZEntropyArray[j]        
        negTotEntropyArray[j] = -(negYWEntropy - negXZEntropy)
    return (negTotEntropyArray)     

 
####################################################################################################
####################################################################################################
#
# Function to create the activation enthalpy of a simplle Ising system. 
#
####################################################################################################
####################################################################################################
    
def createActivationEnthalpyValues(xArray, activEnthalpyArray, eps0, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        x = xArray[j] 
        activEnthalpyArray[j] = eps0*x
    return (activEnthalpyArray)     

 
####################################################################################################
####################################################################################################
#
# Function to create the free energy of a simplle Ising system. 
#
####################################################################################################
####################################################################################################
    
def createInteractEnthalpyValues(xArray, interactEnthalpyArray, eps1, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        x = xArray[j]
        interactEnthalpyArray[j] = -eps1*x*x 
    return (interactEnthalpyArray)     
 
    

####################################################################################################
####################################################################################################
#
# Function to create the free energy of a simplle Ising system. 
#
####################################################################################################
####################################################################################################
    
def createSimpleIsingValues(activEnthalpyArray, interactEnthalpyArray, negXEntropyArray, freeEnergyArray, xTotalSteps, xStep, xIncr):
    
    for j in range (0,xTotalSteps, xStep):
        negXEntropy = negXEntropyArray[j]
        activEnthalpy = activEnthalpyArray[j] 
        interactEnthalpy = interactEnthalpyArray[j]
        freeEnergyArray[j] = activEnthalpy + interactEnthalpy + negXEntropy
    return (freeEnergyArray)          



####################################################################################################
####################################################################################################
#
# Procedure to print the results of equilibrium thermodynamic values when the interaction enthalpy = 0. 
#
####################################################################################################
####################################################################################################


def plotAndPrintEqulibriumResults (xArray, negXEnt, negYEnt, negZEnt, 
                               negYWEnt, negXZEnt, negTotEnt, activEnthalpy, interactEnthalpy, freeEnergy, eps0, eps1, xTotalSteps, xStep):



    #  CVM 2-D Entrop terms           
    print () 
    print (' Equilibrium results, where the interaction enthalpy = 0;')
    print ()  
    print ('    x   negXEntropy negYEntropy negZEntropy negYWEntropy negXZEntropy negTotEnt' )               
    print ()         
    for k in range (0, xTotalSteps, xStep):
        print ('   %.2f' % (xArray[k]), '    %.2f' % (negXEnt[k]), '     %.2f' % (negYEnt[k]), '     %.2f' % (negZEnt[k]),  '      %.2f' % (negYWEnt[k]), '      %.2f' % (negXZEnt[k]), '      %.4f' % (negTotEnt[k]) )        
    print ()

    #  Simple Ising model free energy and other thermodynamic terms              
    print () 
    print (' The Simple Ising model;')
    print ()  
    print ('    x   Entropy  activEnthalpy  interactEnthalpy   freeEnergy' )               
    print ()         
    for k in range (0, xTotalSteps, xStep):
        print ('   %.2f' % (xArray[k]), '    %.4f' % (negXEnt[k]), '     %.4f' % (activEnthalpy[k]), '       %.4f' % (interactEnthalpy[k]),  '       %.4f' % (freeEnergy[k]) )        
    print ()


                                                                      
    pylab.figure(1)
    pylab.plot (xArray, negXEnt)          
    pylab.plot (xArray, negYEnt, 'g')
    pylab.plot (xArray, negZEnt, 'm')            
    print ()  
    print () 
    print ('  The negative xEntropyArray is in blue,' )
    print ('  The negative yEntropyArray  is in green,'  )
    print ('  The negative zEntropyArray  is in maroon,'  )  
    print ()             
    pylab.show()  


    pylab.figure(2)
    pylab.plot (xArray, negYWEnt, 'c')          
    pylab.plot (xArray, negXZEnt, 'y')
    pylab.plot (xArray, negTotEnt, 'k')        
    print ()  
    print () 
    print ('  The combined W and doubled Y negative entropies are in cyan,' )
    print ('  The combined X and doubled Z negative entropies  are in yellow, and'  )
    print ('  The negative total entropy is in black.'  ) 
    print ()             
    pylab.show()  


    pylab.figure(3)
    pylab.plot (xArray, negXEnt)         
    pylab.plot (xArray, negTotEnt, 'k')        
    print ()  
    print () 
    print ('  The negative xEntropyArray is in blue, and' )
    print ('  The negative total entropy is in black.'  )  
    print ()             
    pylab.show() 


    pylab.figure(4)
    pylab.plot (xArray, negXEnt)             
    print ()  
    print () 
    print ('  The negative xEntropyArray is in blue' ) 
    print ()             
    pylab.show()   


    pylab.figure(5)
    pylab.plot (xArray, negXEnt)          
    pylab.plot (xArray, activEnthalpy, 'm')
    pylab.plot (xArray, interactEnthalpy, 'g')    
    pylab.plot (xArray, freeEnergy, 'k')            
    print ()  
    print () 
    print ('  Epsilon0 is   %.2f' % eps0, ' and epsilon1 is  %.2f' % eps1 )    
    print ('  Note that the interaction energy is the NEGATIVE of epsilon1*x*x.  ' )
    print ('  The negative xEntropyArray is in blue,' )
    print ('  The activation enthalpy array is in maroon,'  )    
    print ('  The interaction enthalpy array is in green,'  )
    print ('  The free energy is in black.')  
    print ()             
    pylab.show()  

      
####################################################################################################
####################################################################################################


def main():

####################################################################################################
# Obtain unit array size in terms of array_length (M) and layers (N)
####################################################################################################                

    welcome()
            

####################################################################################################
# Compute configuration variables
####################################################################################################                


    xTotalSteps = 99
    xIncr = 0.01
    xStep = 1
    eps0 = 1.0
    eps1 = 1.0
    xArray              =np.zeros(xTotalSteps, dtype=np.float)
    negXEntropyArray    =np.zeros(xTotalSteps, dtype=np.float)
    negYEntropyArray    =np.zeros(xTotalSteps, dtype=np.float)
    negZEntropyArray    =np.zeros(xTotalSteps, dtype=np.float)   
    negXZEntropyArray   =np.zeros(xTotalSteps, dtype=np.float)   
    negYWEntropyArray   =np.zeros(xTotalSteps, dtype=np.float)   
    negTotEntropyArray  =np.zeros(xTotalSteps, dtype=np.float) 
    activEnthalpyArray  =np.zeros(xTotalSteps, dtype=np.float)  
    interactEnthalpyArray  =np.zeros(xTotalSteps, dtype=np.float)     
    freeEnergyArray     =np.zeros(xTotalSteps, dtype=np.float)      
    xArray              = createXValues(xArray, xTotalSteps, xStep, xIncr)
    negXEntropyArray    = createNegXEntropyValues(xArray, negXEntropyArray, xTotalSteps, xStep, xIncr)
    negYEntropyArray    = createNegYEntropyValues(xArray, negYEntropyArray, xTotalSteps, xStep, xIncr)
    negZEntropyArray    = createNegZEntropyValues(xArray, negZEntropyArray, xTotalSteps, xStep, xIncr)    
    negYWEntropyArray   = createNegYWEntropyValues(xArray, negYEntropyArray, negYWEntropyArray, xTotalSteps, xStep, xIncr)
    negXZEntropyArray   = createNegXZEntropyValues(xArray, negXEntropyArray, negZEntropyArray, negXZEntropyArray, xTotalSteps, xStep, xIncr)
    negTotEntropyArray  = createNegTotEntropyValues(xArray, negYWEntropyArray, negXZEntropyArray, negTotEntropyArray, xTotalSteps, xStep, xIncr) 
    activEnthalpyArray   = createActivationEnthalpyValues(xArray, activEnthalpyArray, eps0, xTotalSteps, xStep, xIncr)
    interactEnthalpyArray = createInteractEnthalpyValues(xArray, interactEnthalpyArray, eps1, xTotalSteps, xStep, xIncr)
    freeEnergyArray     = createSimpleIsingValues(activEnthalpyArray, interactEnthalpyArray, negXEntropyArray, freeEnergyArray, xTotalSteps, xStep, xIncr)

    plotAndPrintEqulibriumResults (xArray, negXEntropyArray, negYEntropyArray, negZEntropyArray, 
                negYWEntropyArray, negXZEntropyArray, negTotEntropyArray, activEnthalpyArray, 
                interactEnthalpyArray, freeEnergyArray, eps0, eps1, xTotalSteps, xStep)
                                                                                                                    
                                                                                                
####################################################################################################
# Conclude specification of the MAIN procedure
####################################################################################################                
    
if __name__ == "__main__": main()

####################################################################################################
# End program
####################################################################################################  