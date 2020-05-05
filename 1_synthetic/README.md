# Synthetic tests to:

## 1_avo_synthetic.ipynb

  1: Reproduce the Figures 2 and 3 from Gidlow and Smith (1987)
    1.1: Implement Aki and Richard (1980) and Shuey (1985)

  2: Reproduce the Figure 3 from Chapter 4 in Chopra and Castagna (2014).
    2.1: Analyse AVO with the implementations in 1.1
    2.2: Plot RC, wavelet and an offset gather (Normal incidence)
    2.3: Generate and plot an angle gather from Shuey (1985) 2 terms' approximation
    
## 2_avo_synthetic.ipynb

  1: Plot AVO for distinct classes with a model of shale/gas.
    1.1: Plot reflectivity, crossplot of intercept x gradient, Normal Incidence and Angle gathers
  2: Plot AVO for distinct classes with a model of shale/brine.
    1.1: Plot reflectivity and crossplot of intercept x gradient
    
 ## 3_avo_synthetic.ipynb

  1: Plot AVO for distinct classes with a model of shale/gas.
    1.1: Plot reflectivity, crossplot of intercept x gradient, generate angle gathers
    1.2: Fit a second order polynomial with L1 and L2 regularization to the angle gather and compute the intercept and gradient
  2: Add gaussian random noise to the angle gathers and perform 1.2
