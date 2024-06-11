library(dplyr)
library(tidyr)

UpperLowerEstimates <- function(point,SE) {
  
  alpha <- ((point * (1 - point)) / SE^2 - 1) * point
  beta <- (alpha / point) - alpha
  
  lower_estimate <- qbeta(0.025, alpha, beta)
  upper_estimate <- qbeta(0.975, alpha, beta)
  
  return(data.frame(lower_estimate,upper_estimate))
}
  
data <- data.frame(variable = c("a","b","m","w","e","s","o","l","d","z"), 
                   point = c(0.225, 0.6, 0.0005,0.002,0.009,0.128,0.6,0.009,0.05,0.05),
                   SE = c(0.001,0.1,0.001,0.01,0.05,0.05,0.1,0.05,0.1,0.1))

# I am far more confident on the values  I have chosen for parameters like mortality
# in Kenya than I am for things like mosquito birth rate in that area, so the 
# values for SE reflect this - as will the tornado plot.


results <- t(apply(data[,c("point","SE")], 1, function(row) UpperLowerEstimates(row["point"], row["SE"])))
results

combined_df_by_row <- do.call(rbind, results)
newdata <- cbind(variables = data$variable,combined_df_by_row)
newdata

# I use the values obtained in model to get upper and lower estimates for endemic rate.

endemic <- data.frame(variable = c("a","b","m","w","e","s","o","l","d","z"), 
                      point = c(1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8),
                      Upper = c(1.8,1.81,3.44,13.4,2.49,0.98,1.81,1.02,1.84,1.84),
                      Lower = c(1.8,1.81,1.44,0.41,1.74,4.3,1.78,1.93,0,0))

tornado_plot <- ggplot(endemic, aes(x = variable, y = point)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper), color = "purple", size =0.5) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Tornado Plot",
    x = "Variable",
    y = "Estimate"
  )

# Print the plot
print(tornado_plot)


# We can see there is a large discrepancy between how much altering individual values
# for parameters changes the model. This difference reflects a) the variation in
# the uncertainty I had for each parameter, and b) the true difference in influence
# each variable has on the endemic level of infection within the population.

