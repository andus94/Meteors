library("geosphere")
######
# Calculates Haversine Distance
#
# @param x: longitude/latitude of point(s). 
#           Can be a vector of two numbers, a matrix of 2 columns
# @param y: as above
# @param r: Radius; default = 6051.8 km, radius of Venus
# return:   Haversine Distance
######
dis <- function(x, y, r = 6051.8){
  Haversine_d <- distHaversine(x, y, r)
  return(Haversine_d)
}

K_function <- function(spatial_points, R, Dm = seq(1,100) * pi * 6051.8 / 100){
  # plyr parallelizes the respective matrix computations
  require(plyr)
  k_values <- vector()
  n <- nrow(spatial_points)
  distance <- list()
  
  constant <- (4 * pi * R^2) / (n * (n - 1))
  
  # Create list of lists of distances of each respective point
  for(i in 1:n){
    
    distance[[i]] <- dis(spatial_points[i,], spatial_points[-i,], R)
    
  }
  
  # Find where distances are less than Dm
  index <- list()
  index <- lapply(Dm, FUN = function(x) distance[[i]] <= x)
  summation <- lapply(index, sum)
  summation <- unlist(summation)
  
  for (i in 2:n){
    index <- lapply(Dm, FUN = function(x) distance[[i]] <= x)
    temp <- lapply(index, sum)
    x <- data.frame(summation, unlist(temp))
    summation <- apply(x, 1, FUN = sum)
  }
  
  summation <- constant * summation
  return(summation)
}

data <- read.delim(file = "http://www.lpi.usra.edu/venus/craters/rel3main.txt")

# Clean Data
column_names <- as.character(unlist(data[2,]))
final_data <- data[3:nrow(data),]
colnames(final_data) <- column_names
longitude <- as.numeric(as.character(final_data$Lon))
latitude <- as.numeric(as.character(final_data$Lat))
coordinates <- cbind(longitude, latitude)
Dm = seq(1,100) * pi * 6051.8 / 100

Ripleys_K <- K_function(coordinates, 6051.8)
CSR_K <- (2 * pi * 6051.8^2) * (1- cos(Dm / 6051.8))
plot(1:100, Ripleys_K - CSR_K, xlab = "1:100", main = "Ripleys_K - CSR_K for Elevation of Craters on Venus", type = "l")


Uniform_Venus <- function(Ripleys_K, n = 100){
  
  R <- vector()
  lon <- vector()
  lat <- vector()
  Sample_K <- matrix(nrow = 100, ncol = n)
  # Generate Uniform Distributions
  for (i in 1:n){
    R <- runif(942, min = -6051.8, max = 6051.8)
    lon <- runif(942, min = 0, max = 360)
    lat <- asin(R / 6051.8) * 180 / pi
    coordinates <- data.frame(lon,lat)
    Sample_K[,i] <- K_function(coordinates, 6051.8)  
  }
  # Build dataframe
  bounds <- apply(Sample_K, MARGIN = 1,FUN = function(x) quantile(x, c(0.05,0.95)))
  lb <- bounds[1,]
  ub <- bounds[2,]
  Dm = seq(1,100) * pi * 6051.8 / 100
  CSR_K <- (2 * pi * 6051.8^2) * (1- cos(Dm / 6051.8))
  
  df <- data.frame(1:100, Sample_K)
  names(df) <- c("indicator", "Sample_K")
  # Plot lower and upper bounds for each dm
  library("ggplot2")
  ggplot(df, aes(x = indicator, y = Sample_K)) + geom_line() + geom_ribbon(aes(x = 1:100, ymin = lb - CSR_K, ymax = ub - CSR_K), alpha = .2)
}

