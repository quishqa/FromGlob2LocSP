library(qualR)

user_name <- "XXXXXXXXXX"
user_pass <- "XXXXXXXXXX"

start_date <- "01/01/1990"
end_date <- "01/09/2020"

start_date2 <- "30/08/2020"

pin_code <- 99
pdp_code <- 72
ibi_code <- 83

o3_code <- 63
co_code <- 16
nox_code <- 18
so2_code <- 13
pm10_code <- 12
pm25_code <- 57
  
  

ibi_o3_30 <- CetesbRetrieve(user_name, user_pass, o3_code,
                            ibi_code, start_date, end_date)
save(ibi_o3_30, file = "ibi_o3_30_year.Rda")

ibi_nox_30 <- CetesbRetrieve(user_name, user_pass, nox_code,
                             ibi_code, start_date, end_date)
save(ibi_nox_30, file = "ibi_nox_30_year.Rda")
  
ibi_co_30 <- CetesbRetrieve(user_name, user_pass, co_code, 
                            ibi_code, start_date, end_date)
save(ibi_co_30, file = "ibi_co_30_year.Rda")

ibi_pm10_30 <- CetesbRetrieve(user_name, user_pass, pm10_code, 
                             ibi_code, start_date, end_date)
save(ibi_pm10_30, file = "ibi_pm10_30_year.Rda")

ibi_pm25_30 <- CetesbRetrieve(user_name, user_pass, pm25_code, 
                             ibi_code, start_date, end_date)
save(ibi_pm25_30, file = "ibi_pm25_30_year.Rda")

ibi_so2_30 <- CetesbRetrieve(user_name, user_pass, so2_code, 
                             ibi_code, start_date, end_date)
save(ibi_so2_30, file = "ibi_so2_30_year.Rda")
  
  
save(ibi_co_30, 
     ibi_nox_30,
     ibi_o3_30,
     ibi_pm10_30,
     ibi_pm25_30,
     ibi_so2_30,
     file = "ibi_30years.Rda")

# pin_o3_30 <- CetesbRetrieve(user_name, user_pass, o3_code, 
#                             pin_code, start_date, end_date)
# 
# pin_co_30 <- CetesbRetrieve(user_name, user_pass, co_code, 
#                             pin_code, start_date, end_date)
# 
# pin_nox_30 <- CetesbRetrieve(user_name, user_pass, nox_code, 
#                              pin_code, start_date, end_date)
# 
# pdp_o3_30 <- CetesbRetrieve(user_name, user_pass, o3_code, 
#                             pdp_code, start_date, end_date)
# 
# pdp_co_30 <- CetesbRetrieve(user_name, user_pass, co_code, 
#                             pdp_code, start_date, end_date)
