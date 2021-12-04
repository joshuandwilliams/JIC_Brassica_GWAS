Step <- 0.1
X <- 1
Y <- 2 # Initial values (arbitrary)
for(i in 1:100){
  M <- (2*X) - Y # Differential equation (arbitrary)
  X <- X + Step
  Y <- Y + (Step * M)
  cat("M: ", M, " X: ", X, " Y: ", Y, "\n")
}

# Adaptive step size function to estimate how complicated curve is at given point and adapt step size

# Plot against known function 

X <- c(0)
Y <- c(0)
Step <- 1
for(i in 2:100){
  M <- (X[i-1]^2) # Differential equation (arbitrary)
  X[i] <- X[i-1] + Step
  Y[i] <- Y[i-1] + (Step * M)
  cat("M: ", M, " X: ", X[i], " Y: ", Y[i], "\n")
}

curve((x^3)/3, from=0, to=10, n=100, col="grey")
lines(X, Y, col="red")

install.packages("plot3D")
library("plot3D")
#http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization
# Try with higher dimensions
t <- c(0)
Y <- c(0)
Z <- c(0)
Step <- 0.01
for(i in 2:100){
  # Think of it as X, Y, Z over time instead
  TY <- 2*(t[i-1]^2) -5
  TZ <- t[i-1]^2 + (2*t[i-1]) #My understanding is that you can define as many differential equations as you like - doesn't have to be one for each 2D plane
  t[i] <- t[i-1] + Step
  Y[i] <- Y[i-1] + (Step * TY)
  Z[i] <- Z[i-1] + (Step * TZ)
}
# Plot this in 3d
#par(mfrow=c(1,1))
scatter3D(x=t, y=Y, z=Z, type="l", bty="g", theta=135, phi=0, lwd=4, main="3D Plot of Euler's Method Example",
          xlab = "t", ylab = "Y", zlab = "Z", ticktype="detailed")
par(mfrow=c(1,1))


# Runge Kutta (RK4)
# y(x + h) = y(x) + 1/6(F1 + 2F2 + 2F3 + F4)
# F1 = hf(x,y)
# F2 = hf(x+h/2, y+F1/2)
# F3 = hf(x+h/2, y+F2/2)
# F4 = hf(x+h, y+F3)

# Example: Use Runge-Kutta Method of Order 4 to solve the following, using a stepsize of h=0.1 for 0<=x<=1
# dy/dx = ((5x^2)-y) / e^(x+y), where y(0)=0

# https://www.intmath.com/differential-equations/12-runge-kutta-rk4-des.php



H = 0.1
Y = 1
for(X in seq(from = 0, to = 1, by=H)){
  cat("X", X, "Y", Y, "\n")
  
  # F1
  F1 <- 0.1 * ( ((5*(X^2))-Y) / exp(X + Y) )
  X2 <- X + (H / 2)
  Y2 <- Y + (F1 / 2)

  # F2
  F2 <- 0.1 * ( ((5*(X2^2))-Y2)) / exp(X2 + Y2)
  X3 <- X + (H / 2)
  Y3 <- Y + (F2 / 2)
  
  # F3
  F3 <- 0.1 * ( ((5*(X3^2))-Y3)) / exp(X3 + Y3)
  X4 <- X + H
  Y4 <- Y + F3
  # F4
  F4 <- 0.1 * ( ((5*(X4^2))-Y4)) / exp(X4 + Y4)
  
  # Y
  Y <- Y + ((1/6)*(F1 + (2 * F2) + (2 * F3) + F4))
}

# Numerical recipes textbook and download scripts


# Storing variables so that I can plot

H <- 0.0001 # Stepsize
Y <- 0 # Initial Y
Yvals <- vector() # Vector to contain Y values
Xvals <- seq(from=0, to=1, by=H) # Vector of X values
for(index in 1:length(Xvals)){
  X <- Xvals[index]
  #cat("X", X, "Y", Y, "\n")
  Yvals[index] <- Y
  
  # F1
  F1 <- 0.1 * ( ((5*(X^2))-Y) / exp(X + Y))
  X2 <- X + (H / 2)
  Y2 <- Y + (F1 / 2)
  
  # F2
  F2 <- 0.1 * ( ((5*(X2^2))-Y2)) / exp(X2 + Y2)
  X3 <- X + (H / 2)
  Y3 <- Y + (F2 / 2)
  
  # F3
  F3 <- 0.1 * ( ((5*(X3^2))-Y3)) / exp(X3 + Y3)
  X4 <- X + H
  Y4 <- Y + F3
  # F4
  F4 <- 0.1 * ( ((5*(X4^2))-Y4)) / exp(X4 + Y4)
  
  # Y
  Y <- Yvals[index] + ((1/6)*(F1 + (2 * F2) + (2 * F3) + F4))
}

plot(Xvals, Yvals, type="l")

# Mask function


# Generate 2 more vectors with actual integral and plot that

# Perhaps can be done with ggplot



