require(tidyverse)

#Source required files
vte_filepath <- "R"
glue::glue("{vte_filepath}/aipw.R") %>% source()
glue::glue("{vte_filepath}/tmle.R") %>% source()

### Construct some data I use the data in the example helpers file 
glue::glue("{vte_filepath}/example_helpers.R") %>% source()

set.seed(1234)
N <- 500 #size of generated data
df <- generate_data_simple(N)

#Next we will fit the working models
#We use the function `fitMods` in the example helpers file
#At the moment the best way of modifying/ using this code is to write a new fitMods function
df_fit <- fitMods(df)
print(df_fit)

#This returns a data frame that has 5 necessary columns:
#Y - the outcome of interest
#A - the treatment indicator
#pi_hat - the predicted propensity score E(A|X)
#mu1_hat - the predicted outcome regression in the treated   E(Y|A=1,X)
#mu0_hat - the predicted outcome regression in the untreated E(Y|A=0,X)

#Often it is better to use cross fitted models (see e.g. CV-TMLE notes in Levy's paper)
#This is implemented in the `crossFit' function which essentially calls fitmods multiple times
foldIDs <- getFoldIDs(N,Nfolds=5) #5 fold cross validation in this example
df_xfit <- crossFit(df,foldIDs=foldIDs)
print(df_xfit) #Note there is an additional FoldID column

#These fitted data frames can be passed used to infer ATE and VTEs
#We also infer the square-root of the VTE which is on the same scale as the ATE
#Two methods are provided:
#AIPW_VTE() uses one step bias correction estimators
#TMLE_VTE() uses one TMLE
AIPW_VTE(df_fit)
TMLE_VTE(df_fit)

AIPW_VTE(df_xfit)
TMLE_VTE(df_xfit)



