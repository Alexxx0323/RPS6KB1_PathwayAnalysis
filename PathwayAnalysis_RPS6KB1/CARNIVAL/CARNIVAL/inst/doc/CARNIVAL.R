## ----setup, echo=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE, cache=FALSE----
library(knitr)
opts_chunk$set(
  cache = FALSE,
  echo = TRUE,
  warning = FALSE,
  error = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
# ~/Applications/IBM/ILOG/CPLEX_Studio129/cplex/bin/x86-64_osx/cplex

## -----------------------------------------------------------------------------
# /opt/ibm/ILOG/CPLEX_Studio129/cplex/bin/x86-64_linux/cplex

## -----------------------------------------------------------------------------
# C:/Program Files/IBM/ILOG/CPLEX_Studio129/cplex/bin/x64_win64/cplex.exe

## -----------------------------------------------------------------------------
library(CARNIVAL)

load(file = system.file("toy_inputs_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_measurements_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_network_ex1.RData",
                        package="CARNIVAL"))

# lpSolve
result = runCARNIVAL(inputObj = toy_inputs_ex1, measObj = toy_measurements_ex1,
                     netObj = toy_network_ex1)

print(result)

## -----------------------------------------------------------------------------
library(CARNIVAL) # load CARNIVAL library

load(file = system.file("toy_measurements_ex2.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_network_ex2.RData",
                        package="CARNIVAL"))

# lpSolve
result = runCARNIVAL(measObj = toy_measurements_ex2, netObj = toy_network_ex2)

print(result)

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

