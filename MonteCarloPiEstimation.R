# Longer, clearer version 11 lines
radius <- 1
corner_pos <- 1
sample_sizes <- seq(10^2, 10^6, length.out = 1000)
pi_estimates <- vector()
for(iteration in 1:length(sample_sizes)){
  x <- sample(runif(10^5, -corner_pos, corner_pos), sample_sizes[iteration], replace = T)
  y <- sample(runif(10^5, -corner_pos, corner_pos), sample_sizes[iteration], replace = T)
  fraction_samples_in_circle <- (sum(x^2 + y^2 <= radius^2)) / sample_sizes[iteration]
  pi_estimates[iteration] <- ((fraction_samples_in_circle) * ((2*corner_pos)^2)) / radius^2
}
plot(y=pi_estimates, x=sample_sizes, type="l")


# Condensed version (lower readibility) 7 lines
sample_sizes <- seq(10^2, 10^6, length.out = 1000)# A sequence of sample sizes
pi_estimates <- vector()
for(iteration in 1:length(sample_sizes)){
  pos <- list(sample(runif(10^5, -1, 1), sample_sizes[iteration], replace = T), sample(runif(10^5, -1, 1), sample_sizes[iteration], replace = T)) # Taking a sample of random coordinates
  pi_estimates[iteration] <- (sum(pos[[1]]^2 + pos[[2]]^2 <= 1) / sample_sizes[iteration]) * 4  # This is the pi estimate
}
plot(y=pi_estimates, x=sample_sizes, type="l")





for(N in 1:6){ # Number of N values
  in_circle_counter <- 0
  sample_size <- 10^N
  for(position in 1:sample_size){
    x <- runif(1, -1, 1)
    y <- runif(1, -1, 1)
    if(x^2 + y^2 <= 1){
      in_circle_counter <- in_circle_counter + 1
    }
  }
  print((in_circle_counter / sample_size) * 4)
}

# x is same as ((Nb / (Nb+Na)) * AT) / R^2, but R^2 = 1 here so can be ignored
# The 4 is because the total area is a square of sides 2 (around a circle of radius 1)

# It is still taking long to run for the higher sample sizes...



