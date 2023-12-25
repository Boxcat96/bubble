rm(list = ls())
graphics.off()

library(tidyverse)
library(readxl)
library(openxlsx)
library(mFilter)

# folder
setwd("C:/Users/tkero/OneDrive/経済ファイル/02_分析/04_Bubble/rawdata")

# read data
dat = read_xlsx("DATA_master.xlsx")
head(dat)

dat<-dat %>% 
  filter(year >= ymd("1994-01-01")) %>% 
  filter(year <= ymd("2019-12-31"))
head(dat)

# rename data
y = dat$GDP_yoy
credit = dat$credit_to_GDP
stock = dat$stock_to_GDP

# hp-filter
hp_param = 1600
hp_credit = hpfilter(credit,freq = hp_param)
hp_stock = hpfilter(stock,freq = hp_param)

# calculate percentage of deviation from trend
hp_credit_percent = (hp_credit$cycle/hp_credit$trend)*100
hp_stock_percent = (hp_stock$cycle/hp_stock$trend)*100

# output
data = cbind(y, hp_credit_percent, hp_stock_percent)

# export csv
write.csv(data, "HP.csv")
plot(hp_stock)
