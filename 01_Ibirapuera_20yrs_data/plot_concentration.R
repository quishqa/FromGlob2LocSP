load("ibi_30years.Rda")

library(openair)
library(scales)

# Calulating mean value from 2000 to 2020
SelectFrom20DayAvg <- function(df, stats = "mean"){
  df_2000 <- selectByDate(df, year = 2000:2020)
  df_m <- timeAverage(df_2000, avg.time = "month", statistic = stats)
  return(df_m)
}

# Selecting pandemic months
SelectCov <- function(df){
  df_m <- SelectFrom20DayAvg(df)
  df_cov <- selectByDate(df_m, year = 2020, month = 3:9)
  return(df_cov)
}


# 20 years trens
PlotTemporalTrend <- function(df, ylabel, col){
  df_m <- SelectFrom20DayAvg(df)
  df_m_sd <- SelectFrom20DayAvg(df, stats = "sd")
  df_cov <- SelectCov(df)
  trend <- lm(pol ~ date, df_m)
  years <- seq(as.POSIXct("2000-01-01", tz = "UTC"),
               as.POSIXct("2020-01-01", tz = "UTC"),
               by = "5 year")
  y_min <- min(df_m$pol - df_m_sd$pol, na.rm = T)
  y_max <- max(df_m$pol + df_m_sd$pol, na.rm = T)
  
  plot(df_m$date, df_m$pol, t = "n", 
       ylim = c(y_min, y_max),
       xlab = "", ylab = "",
       axes = F)
  # abline(model, col = col, lwd = 1.3)
  lines(df_m$date, df_m$pol, col = col, lwd = 1.35)
  points(df_cov$date, df_cov$pol, pch = 19, cex = 1.25)
  points(df_m$date, df_m$pol, pch = 19, col = col)
  arrows(df_m_sd$date, (df_m$pol - df_m_sd$pol),
         df_m_sd$date, (df_m$pol + df_m_sd$pol),
         angle=90, length = 0, col = alpha(col, 0.5))
  axis(2, cex.axis = 1.8)
  axis(1, at = years, labels = seq(2000, 2020, 5),
       cex.axis = 1.8)
  title(ylab = ylabel,
        line = 2.3, cex.lab = 1.8)
  box()
}


pdf("ibirapuera_20years.pdf", width = 14, height = 10)
par(mfrow=c(3, 2), mar=c(2, 4.6, 1.5, 0.7), 
    oma = c(0, 0, 1, 0))
PlotTemporalTrend(ibi_o3_30, 
                  expression("O"[3] * " (" * mu * "g m" ^-3 * ")"),
                  "blue")
PlotTemporalTrend(ibi_nox_30, 
                  expression("NO"[X] * " (ppb)"),
                  "green")
PlotTemporalTrend(ibi_co_30, 
                  expression("CO" * " (ppm)"),
                  "red")
PlotTemporalTrend(ibi_so2_30, 
                  expression("SO"[2] *  " (" * mu * "g m" ^-3 * ")"),
                  "chocolate")
PlotTemporalTrend(ibi_pm10_30, 
                  expression("PM"[10] * " (" * mu * "g m" ^-3 * ")"),
                  "darkorange")
PlotTemporalTrend(ibi_pm25_30, 
                  expression("PM"[2.5] * " (" * mu * "g m" ^-3 * ")"),
                  "cadetblue")
title("Ibirapuera", line = -0.5, outer = TRUE, cex.main = 2)
dev.off()



load("./pin_pm10_pm25_so2.Rda")
load("./pin_tol_30_year.Rda")

pdf("pinheiros_20years.pdf", width = 14, height = 10)
par(mfrow=c(3, 2), mar=c(2, 4.6, 1.5, 0.7), 
    oma = c(0, 0, 1, 0))
PlotTemporalTrend(pin_o3_30, 
                  expression("O"[3] * " (" * mu * "g m" ^-3 * ")"),
                  "blue")
PlotTemporalTrend(pin_nox_30, 
                  expression("NO"[X] * " (ppb)"),
                  "green")
PlotTemporalTrend(pin_co_30, 
                  expression("CO" * " (ppm)"),
                  "red")
PlotTemporalTrend(pin_tol_30, 
                  expression("TOL" *  " (" * mu * "g m" ^-3 * ")"),
                  "chocolate")
PlotTemporalTrend(pin_pm10_30, 
                  expression("PM"[10] * " (" * mu * "g m" ^-3 * ")"),
                  "darkorange")
PlotTemporalTrend(pin_pm25_30, 
                  expression("PM"[2.5] * " (" * mu * "g m" ^-3 * ")"),
                  "cadetblue")
title("Pinheiros", line = -0.5, outer = TRUE, cex.main = 2)
dev.off()





