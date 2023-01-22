transform_csv <- function(data){
  data = select(data, "Open","Date")
  data = rename(data, c("Prix"="Open"))
  data$Prix = as.numeric(data$Prix)
  data$Date = as.Date(data$Date)
  
  data = mutate(data, rendement = log(data$Prix/lag(data$Prix)))
  data = mutate(data, rendement2 = rendement**2)
  return(data[-1,])
}

plot_series_temp <-function(data){
  plot_prix_temps = ggplot(data = data) + geom_line(aes(x = Date,y = Prix))
  plot(plot_prix_temps)
  
  plot_rendement_temps = ggplot(data = data) + geom_line(aes(x = Date,y = rendement))
  plot(plot_rendement_temps)
  
  plot_rendement2_temps = ggplot(data = data) + geom_line(aes(x = Date,y = rendement2))
  plot(plot_rendement2_temps)
}

autocorrelations <-function(data){
  acf(data$rendement,type='correlation', na.action=na.pass, plot=TRUE)
  acf(data$rendement2,type='correlation', na.action=na.pass, plot=TRUE)
}