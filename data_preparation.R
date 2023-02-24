library (tidyverse)

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

# New

transform_csv_with_discount <- function(data,r){ #similar as transform_csv
  data = select(data, "Open","Date")
  data = rename(data, c("Prix"="Open"))
  
  data$Prix = as.numeric(data$Prix)
  data$Date = as.Date(data$Date)
  
  # Discounting module
  n = length(data$Prix)
  discount_vector = double(n)
  discount_vector[1]=1
  for(i in 2:n){
    discount_vector[i] = discount_vector[i-1]*(1+r)
  }
  data$Prix = data$Prix / discount_vector
  
  # End of transform_csv
  data = mutate(data, rendement = log(data$Prix/lag(data$Prix)))
  data = mutate(data, rendement2 = rendement**2)
  return(data[-1,])
  
  
  
}

