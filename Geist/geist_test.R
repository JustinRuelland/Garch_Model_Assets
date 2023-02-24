source(file = "./data_preparation.R",local= TRUE)
data = read.csv("./CAC40_15_19.csv")
data <- transform_csv_with_discount(data,0.001)
