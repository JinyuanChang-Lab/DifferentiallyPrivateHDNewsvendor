# =============================================================================
# File: prepare_real_data.R
#
# Purpose:
#   Construct the analysis data used in the empirical application from the
#   original Corporación Favorita Grocery Sales Forecasting files.
#
# Required input files:
#   data/raw/stores.csv
#   data/raw/items.csv
#   data/raw/transactions.csv
#   data/raw/oil.csv
#   data/raw/holidays_events.csv
#   data/raw/train.csv
#
# Output:
#   data/processed/processed_data.rds
#
# Notes:
#   - The original input files must be downloaded from the official Kaggle
#     competition page.
#   - No manual modification of the original CSV files is required.
#   - All filtering, merging, transformation, and feature construction steps
#     used in the paper are implemented in this script.
# =============================================================================


#-------------------------------------------------------------------------------
#                        Required packages                         
#------------------------------------------------------------------------------- 
library(here)
library(dplyr)
library(lubridate)
library(purrr)



#-------------------------------------------------------------------------------
#                         1 Original data                         
#------------------------------------------------------------------------------- 
# Original data can be downloaded online from https://www.kaggle.com/competitions/favorita-grocery-sales-forecasting/data
# The raw data files should be placed in:
# data/raw/

here::i_am("scripts/real_data/prepare_real_data.R")
raw_dir <- here::here("data", "raw")
raw_files <- c(
  store        = "stores.csv",
  item         = "items.csv",
  transactions = "transactions.csv",
  oil          = "oil.csv",
  holiday      = "holidays_events.csv",
  train        = "train.csv"
)
raw_paths <- setNames(
  file.path(raw_dir, unname(raw_files)),
  names(raw_files)
)

raw_data <- lapply(
  raw_paths,
  read.csv,
  stringsAsFactors = FALSE
)
 

## Original feature data  
# 1.stores
Store_raw <- raw_data[["store"]]  
# 2.items & select the seafood family
Item_raw <- raw_data[["item"]] 
Item_raw1 = Item_raw[Item_raw$perishable==1,]
Item_raw2 = Item_raw1[Item_raw1$family=='SEAFOOD',]
SEAFOOD_item_index <- as.character(Item_raw2$item_nbr)
# 3. transactions
Trans_raw <-  raw_data[["transactions"]]  
Trans_raw <- Trans_raw %>% mutate(date = as.Date(date))  
# 4. oil prices
Oil_raw <- raw_data[["oil"]]  
# 5.  holidays_events 
Holiday_raw <- raw_data[["holiday"]]  

## Original demand data  
Train_raw <- raw_data[["train"]]  
Train_raw <- Train_raw[Train_raw$item_nbr %in% SEAFOOD_item_index,]                        # select the seafood family 
Train_raw <- Train_raw[Train_raw$onpromotion != "", ]                                      # delete the rows while missing onpromotion values
earthquake_index <- which(Train_raw$date == '2016-04-16')[1]                               
Train_raw <- Train_raw[1:(earthquake_index-1),]                                            # select the date before earthquake
startdate_index <- which(Train_raw$date == '2014-04-07')[1] 
Train_raw <- Train_raw[startdate_index:nrow(Train_raw),]                                   # select the date after 2014-04-07 (Monday)
Train_raw <- Train_raw %>% mutate(date = as.Date(date)) %>% 
  mutate( 
    week_start = floor_date(date, unit = "week", week_start = 1),    
    day = wday(date, label = TRUE, abbr = FALSE, locale = "C"),
    month = month(week_start, label = TRUE, abbr = FALSE, locale = "C") 
  )                                                                                        # add day\month
rm(raw_data) #delete raw data

#-------------------------------------------------------------------------------
#                         2 Demand data processing                       
#------------------------------------------------------------------------------- 
## construct "week_id"
Train_raw1 <- Train_raw %>%  mutate( date = as.Date(date),
                                     week_start = floor_date(date, unit = "week", week_start = 1),  # week start from Monday
                                     week_id = as.integer(factor(week_start)) )
## add transactions
Train_raw1 <- Train_raw1 %>% left_join(Trans_raw %>% select(date, store_nbr,  transactions),  by = c("date", "store_nbr"))

Sale_week_summary <- Train_raw1 %>% group_by(store_nbr, week_id) %>%  summarise(  total_sold = sum(unit_sales, na.rm = TRUE)  ) # Weekly total sales per store
Sale_week_summary <- Sale_week_summary %>%  group_by(store_nbr) %>% mutate(store_count = n()) %>% ungroup()                     # Number of observed weeks per store
Sale_week_summary1 <- Sale_week_summary %>% left_join(Store_raw %>% select(  store_nbr,  city),  by = c( "store_nbr"))          # Add city info to weekly sales

## Keep cities that have exactly one store 
Store_raw_new <- Store_raw %>%  group_by(city) %>% mutate(city_store_count = n()) %>% ungroup()                                 # Number of stores in each city
save_city <- Store_raw_new[Store_raw_new$city_store_count==1,]$city                                                             # Cities with a single store
Sale_week_summary_select <- Sale_week_summary1[Sale_week_summary1$city %in% save_city,]                                         # Keep only single-store cities

## Exclude stores with fewer than 95 observed weeks (i.e., ≥10% missing/no-sales weeks)
Sale_week_summary_stockout <- Sale_week_summary_select[Sale_week_summary_select$store_count >= 95,] 

## Past-month average demand (4-week rolling mean, requires 4 prior weeks) 
past_month_ave_demand <- Sale_week_summary_stockout %>% group_by(store_nbr) %>% arrange(week_id, .by_group = TRUE) %>%
  mutate(
    past_month_ave_demand = rowMeans(cbind( lag(total_sold, 1), lag(total_sold, 2), lag(total_sold, 3), lag(total_sold, 4) ), na.rm = FALSE)  
  ) %>% ungroup() %>%   filter(!is.na(past_month_ave_demand))

## Demand data, n = 801
demand_data_pre = past_month_ave_demand[,c("store_nbr","week_id","total_sold")]
demand_data = demand_data_pre[,3][[1]]  

## unique index set 
Index <- past_month_ave_demand[,c("store_nbr","week_id")]


#-------------------------------------------------------------------------------
#                         3 Feature data processing                       
#------------------------------------------------------------------------------- 
##  1. past_month_ave_demand  1 variable
past_month_ave_demand_data <- past_month_ave_demand[,c("store_nbr","week_id","past_month_ave_demand")]

##  2. transactions 1 variable
trans_data_pre <- Train_raw1 %>% group_by(store_nbr, week_id) %>%  summarise(total_trans = sum(transactions, na.rm = TRUE))
trans_data <- Index %>% left_join(trans_data_pre %>% select(store_nbr, week_id,  total_trans),  by = c(  "store_nbr", "week_id"))

##  3. month_dummy  11 variables
month_pre <-   unique(Train_raw1[,c("store_nbr","week_id",'month')]  )
month_pre$month <- factor(month_pre$month,  levels = month.name)   
month_dummies <- model.matrix(~ month - 1, data = month_pre)
month_data_pre <- cbind(month_pre, month_dummies) 
month_dummy_data <- Index %>% left_join(month_data_pre,  by = c(  "store_nbr", "week_id"))   %>% select(-c("month","monthJanuary"))

##  4. wage dummy 2 variables
wage_pre <- Train_raw1 %>% mutate(  wage_15 = if_else(day(date) == 15, 1, 0), wage_30 = if_else(date == ceiling_date(date, "month") - days(1), 1, 0) )
wage_data_pre <- wage_pre %>% group_by(store_nbr, week_id) %>%  summarise(across(  c(wage_15, wage_30), ~ as.integer(sum(.x, na.rm = TRUE) > 0), .names = "any_{.col}" ))
wage_dummy_data <- Index %>% left_join(wage_data_pre,  by = c(  "store_nbr", "week_id"))

##  5. proportion of sales from promotional items  1 variable
Prop_prom_data_pre <- Train_raw1 %>%  mutate(promo_flag = if_else(onpromotion == "False", 0, 1))
Prop_prom_data_pre <- Prop_prom_data_pre %>% group_by(store_nbr, week_id) %>%  summarise(average_promotion = mean(promo_flag, na.rm = TRUE))
Prop_prom_data <- Index %>% left_join(Prop_prom_data_pre,  by = c(  "store_nbr", "week_id"))

##  6. store-level based on Store_raw  10 variables 
# store_type_dummy 4 variables
store_type_dummies <- model.matrix(~ type - 1, data = Store_raw)
store_type_dummies_data_pre <- bind_cols(Store_raw[, "store_nbr", drop = FALSE], store_type_dummies[, 1:5])
store_type_dummy_data <- Index %>% left_join(store_type_dummies_data_pre,  by = c(  "store_nbr" )) %>% select(-c(3)) # no type "A"
# store_cluster_dummy 6 variables
store_cluster_dummies <- model.matrix(~ factor(cluster) - 1, data = Store_raw)
colnames(store_cluster_dummies) <- gsub("^factor\\(cluster\\)", "cluster_", colnames(store_cluster_dummies))
store_cluster_dummies_data_pre <- bind_cols(Store_raw[, "store_nbr", drop = FALSE], store_cluster_dummies[, 1:17])
store_cluster_dummy_data <- Index %>% left_join(store_cluster_dummies_data_pre,  by = c(  "store_nbr" ))  %>% select(-c(4,7,9,10,11,13,14,15,16,18,19)) # all values are "0" 

##  7. Oil_price 1 variable 
Oil_price <- Oil_raw  %>%
  mutate(price = ifelse(is.na(dcoilwtico),
                        zoo::rollapply(dcoilwtico, width = 5, FUN = function(x) mean(x[-3], na.rm = TRUE), fill = NA, align = "center"),
                        dcoilwtico)) %>% mutate(date = as.Date(date)) # fit NA values
Oil_data_pre <- Train_raw1 %>%  left_join(Oil_price %>% select(date, price),   by = c("date" ))
Oil_data_pre <- unique(Oil_data_pre[,c("date","week_id",'price')])
Oil_data_pre[is.na(Oil_data_pre)] <- 0
Oil_data_pre <- Oil_data_pre %>% group_by(week_id) %>% summarise(mean_price = mean(price, na.rm = TRUE))
Oil_data <-  Index %>% left_join(Oil_data_pre,  by = c(  "week_id" ))

## 8. holiday 
Holiday_raw <- Holiday_raw %>%  mutate(date = as.Date(date))
# National holiday dummy 5 variables
Holiday_National <- Holiday_raw[Holiday_raw$locale=='National',]
Holiday_National <- Holiday_National %>% arrange(date) %>%  distinct(date, .keep_all = TRUE)
Holiday_National <- Holiday_National %>%
  mutate(
    is_National_Additional = if_else(locale == "National" & type == "Additional", 1, 0),
    is_National_Bridge     = if_else(locale == "National" & type == "Bridge",     1, 0),
    is_National_Event      = if_else(locale == "National" & type == "Event",      1, 0),
    is_National_Holiday    = if_else(locale == "National" & type == "Holiday",    1, 0),
    is_National_Transfer   = if_else(locale == "National" & type == "Transfer",   1, 0),
    is_National_Workday   = if_else(locale == "National" & type == "Work Day",   1, 0)
  )
National_data_pre <- left_join(Train_raw1, Holiday_National %>% select(date,is_National_Additional, is_National_Bridge,is_National_Event,is_National_Holiday,is_National_Transfer,is_National_Workday), by = "date")
National_data_pre[is.na(National_data_pre)] <- 0
National_data_pre <- National_data_pre %>%
  group_by(store_nbr, week_id) %>%
  summarise(across(
    c(is_National_Additional, is_National_Bridge, is_National_Event, 
      is_National_Holiday,  is_National_Workday ),
    ~ as.integer(sum(.x, na.rm = TRUE) > 0),
    .names = "any_{.col}"
  ))
National_dummy_data <-  Index %>% left_join(National_data_pre,  by = c( "store_nbr", "week_id" )) 

## 9. city dummy  8  variables
city_dummies <- model.matrix(~ city - 1, data = Store_raw)
city_dummies_data_pre <- bind_cols(Store_raw[, "store_nbr", drop = FALSE], city_dummies[, 1:22])
city_dummy_data <- Index %>% left_join(city_dummies_data_pre,  by = c(  "store_nbr" ))  %>% select(-c(3,6,9,10,11,13,16,17,18,19,20,21,22 ,24))

# first level
First_feature_pre <- list(past_month_ave_demand_data, trans_data, month_dummy_data, wage_dummy_data,Prop_prom_data,store_type_dummy_data,store_cluster_dummy_data,Oil_data,National_dummy_data ,city_dummy_data)
First_feature_data  <- reduce(First_feature_pre, inner_join, by = c("store_nbr", "week_id")) %>% select(-c("store_nbr", "week_id"))  

### interactions 
keys <- c("store_nbr", "week_id")

## National_dummy_data:Prop_prom_data
NP1 <- as.matrix(National_dummy_data[ , setdiff(names(National_dummy_data), keys), drop = FALSE])
NP2 <- as.matrix(Prop_prom_data[ , setdiff(names(Prop_prom_data), keys), drop = FALSE])
p1 <- ncol(NP1); q1 <- ncol(NP2)
NP <- NP1[, rep(1:p1, each = q1), drop = FALSE] *  NP2[, rep(1:q1, times = p1), drop = FALSE]

## wage_dummy_data:trans_data
WT1 <- as.matrix(wage_dummy_data[ , setdiff(names(wage_dummy_data), keys), drop = FALSE])
WT2 <- as.matrix(trans_data[ , setdiff(names(trans_data), keys), drop = FALSE])
p3 <- ncol(WT1); q3 <- ncol(WT2)
WT <- WT1[, rep(1:p3, each = q3), drop = FALSE] *  WT2[, rep(1:q3, times = p3), drop = FALSE]

## store_type_dummy_data:past_month_ave_demand_data
SP1 <- as.matrix(store_type_dummy_data[ , setdiff(names(store_type_dummy_data), keys), drop = FALSE])
SP2 <- as.matrix(past_month_ave_demand_data[ , setdiff(names(past_month_ave_demand_data), keys), drop = FALSE])
p4 <- ncol(SP1); q4 <- ncol(SP2)
SP <- SP1[, rep(1:p4, each = q4), drop = FALSE] *  SP2[, rep(1:q4, times = p4), drop = FALSE]

## Oil_data:city_dummy_data
OC1 <- as.matrix(Oil_data[ , setdiff(names(Oil_data), keys), drop = FALSE])
OC2 <- as.matrix(city_dummy_data[ , setdiff(names(city_dummy_data), keys), drop = FALSE])
p5 <- ncol(OC1); q5 <- ncol(OC2)
OC <- OC1[, rep(1:p5, each = q5), drop = FALSE] *  OC2[, rep(1:q5, times = p5), drop = FALSE]

## city_dummy_data:month_dummy_data
CM1 <- as.matrix(city_dummy_data[ , setdiff(names(city_dummy_data), keys), drop = FALSE])
CM2 <- as.matrix(month_dummy_data[ , setdiff(names(month_dummy_data), keys), drop = FALSE])
p6 <- ncol(CM1); q6 <- ncol(CM2)
CM <- CM1[, rep(1:p6, each = q6), drop = FALSE] *  CM2[, rep(1:q6, times = p6), drop = FALSE]

## store_cluster_dummy_data:month_dummy_data
SM1 <- as.matrix(store_cluster_dummy_data[ , setdiff(names(store_cluster_dummy_data), keys), drop = FALSE])
SM2 <- as.matrix(month_dummy_data[ , setdiff(names(month_dummy_data), keys), drop = FALSE])
p7 <- ncol(SM1); q7 <- ncol(SM2)
SM <- SM1[, rep(1:p7, each = q7), drop = FALSE] *  SM2[, rep(1:q7, times = p7), drop = FALSE]

## store_cluster_dummy_data:National_dummy_data
SN1 <- as.matrix(store_cluster_dummy_data[ , setdiff(names(store_cluster_dummy_data), keys), drop = FALSE])
SN2 <- as.matrix(National_dummy_data[ , setdiff(names(National_dummy_data), keys), drop = FALSE])
p8 <- ncol(SN1); q8 <- ncol(SN2)
SN <- SN1[, rep(1:p8, each = q8), drop = FALSE] *  SN2[, rep(1:q8, times = p8), drop = FALSE]

## final standardized feature data
Second_feature_data <- cbind(First_feature_data,NP,WT,SP,OC,CM,SM,SN)   
keep <- vapply(Second_feature_data, function(z) dplyr::n_distinct(z, na.rm = TRUE) > 1, logical(1))
Second_feature_data <- Second_feature_data[, keep, drop = FALSE] 

is_dummy01 <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 2 && all(ux %in% c(0, 1))
}

dummy_cols <- vapply(Second_feature_data, is_dummy01, logical(1))

X_cont  <- Second_feature_data[, !dummy_cols, drop = FALSE]
X_dummy <- Second_feature_data[,  dummy_cols, drop = FALSE] 
X_cont_scaled <- scale(X_cont, center = TRUE, scale = TRUE) 
X_dummy_centered <- scale(X_dummy, center = TRUE, scale = FALSE)

Second_feature_data_scaled <- cbind(
  as.data.frame(X_cont_scaled),
  as.data.frame(X_dummy_centered)
)

Second_feature_data_scaled <- as.matrix(sapply(Second_feature_data_scaled, as.numeric))


 


#-------------------------------------------------------------------------------
#                         4 Final data                         
#------------------------------------------------------------------------------- 
# demand: demand_data, n = 801
# feature: Second_feature_data_scaled, dim = (801,242), without intercept
RD_Y <- demand_data
RD_X <- Second_feature_data_scaled 
 
processed_dir <- here::here("data", "processed")
dir.create( processed_dir, recursive = TRUE, showWarnings = FALSE )
processed_data <- list( RD_Y = RD_Y, RD_X = RD_X )
saveRDS( processed_data, file = file.path(processed_dir, "processed_data.rds"), compress = "gzip")
