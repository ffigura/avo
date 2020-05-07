# Synthetic tests to:

## 1_avo_synthetic.ipynb

  1: Reproduce the Figures 2 and 3 from Gidlow and Smith (1987)\
  &nbsp;&nbsp;&nbsp;1.1: Implement Aki and Richard (1980) and Shuey (1985).

  2: Reproduce the Figure 3 from Chapter 4 in Chopra and Castagna (2014).\
  &nbsp;&nbsp;&nbsp;2.1: Analyse AVO with the implementations in 1.1.\
  &nbsp;&nbsp;&nbsp;2.2: Plot RC, wavelet and an offset gather (Normal incidence).\
  &nbsp;&nbsp;&nbsp;2.3: Generate and plot an angle gather from Shuey (1985) 2 terms' approximation.
    
## 2_avo_synthetic.ipynb

  1: Plot AVO for distinct classes with a model of shale/gas.\
  &nbsp;&nbsp;&nbsp;1.1: Plot reflectivity, crossplot of intercept x gradient, Normal Incidence and Angle gathers.
    
  2: Plot AVO for distinct classes with a model of shale/brine.\
  &nbsp;&nbsp;&nbsp;2.1: Plot reflectivity and crossplot of intercept x gradient.
    
 ## 3_avo_synthetic.ipynb

  1: Plot AVO for distinct classes with a model of shale/gas.\
  &nbsp;&nbsp;&nbsp;1.1: Plot reflectivity, crossplot of intercept x gradient, generate angle gathers.\
  &nbsp;&nbsp;&nbsp;1.2: Fit a second order polynomial with L1 and L2 regularization to the angle gather and compute the intercept and gradient.
    
  2: Add gaussian random noise to the angle gathers.\
  &nbsp;&nbsp;&nbsp;2.1: Perform 1.2.

 ## 4_avo_synthetic.ipynb

  1: Compute the elastic impedance (Connolly, 1999) and the normalized elastic impedance (Whitcombe, 2002).\
  &nbsp;&nbsp;&nbsp;1.1: Plot input logs from 4 AVO scenarios with shale, gas sand and brine sand.\
  &nbsp;&nbsp;&nbsp;1.2: Plot crossplots with facies defining the thrid dimension.
