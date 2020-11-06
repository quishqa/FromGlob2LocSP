library(qualR)


# CETESB QUALAR Credentials
user_name <- "XXXXXXXXXX"
user_pass <- "XXXXXXXXXX"

# Start and end date to download
start_date <- "01/01/1990"
end_date <- "01/09/2020"

# Pinheiros and Ibirapuera AQS codes
pin_code <- 99
ibi_code <- 83

# Polluntas codes
o3_code <- 63
co_code <- 16
nox_code <- 18
so2_code <- 13
pm10_code <- 12
pm25_code <- 57
  
  
# Downloading Ibirapuera data

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

# Downloading data for Pinheiros AQS
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
