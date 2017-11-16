
###  Fengler09LV.R  ###


rm(list = ls())

source("./source/BSformulas.R")
source("./source/Fengler09/pre_smoother.R")


# load/install packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(  #  
               )

# load dataset  
load("./Rdata/data_list.RData")
source("./source/plotter.R")


Fengler09_getLV <- function (l, moneyness, k_grid, expiries = NULL, iv = NULL, 
                          rates = NULL, dividends = NULL, spot = NULL, ...){
  # 
  # The function perform... 
  #
  # ARGUMENTS: 
  # * 
  # 
  # VALUE: 
  # 
  #
  
  
}

Fengler09_LV <- function (l, moneyness, k_grid, expiries = NULL, iv = NULL, 
                           rates = NULL, dividends = NULL, spot = NULL, ...){
  # 
  # The function perform... 
  #
  # ARGUMENTS: 
  # * 
  # 
  # VALUE: 
  # 
  #
  
  
  
}

