library(dplyr)

create_regions <- function(num_parts_axis, min_x, min_y, max_x, max_y, total_length, data) {
  M <- num_parts_axis
  local_df <- list()
  xy_boundary <- list()
  xy_center_coor <- list()
  add <- total_length / M # add segments
  
  for (jj in 1:M) {
    for (ii in 1:M) {
      region_test11 <- data %>% filter(X > min_x + (ii - 1) * add)
      region_test12 <- region_test11 %>% filter(X < min_x + ii * add)
      region_test13 <- region_test12 %>% filter(Y > min_y + (jj - 1) * add)
      region_test14 <- region_test13 %>% filter(Y < min_y + jj * add)
      max_min_vec <- c(min_x + (ii - 1) * add, min_x + ii * add, min_y + (jj - 1) * add, min_y + jj * add)
      x_y_center_coor_vec <- c((min_x + (ii - 1) * add + min_x + ii * add) / 2, (min_y + (jj - 1) * add + min_y + jj * add) / 2)
      local_df[[length(local_df) + 1]] <- region_test14 # store "regions" data
      xy_boundary[[length(xy_boundary) + 1]] <- max_min_vec # store corresponding x and y boundaries
      xy_center_coor[[length(xy_center_coor) + 1]] <- x_y_center_coor_vec
    }
  }
  
  plot(data$X, data$Y, cex = 0.01, xlim=c(min_x,max_x),ylim=c(min_y,max_y))
  plot_regions_from_whole(num_parts_axis = num_parts_axis, min_x = min_x, min_y = min_y, max_x = max_x, max_y = max_y, total_length = total_length)
  
  print(length(local_df))    # length check
  print(length(xy_boundary)) # length check
  print(length(xy_center_coor))
  
  return(list(local_df = local_df, xy_boundary = xy_boundary, xy_center_coor = xy_center_coor))
}


plot_regions_from_whole <- function(num_parts_axis, min_x, min_y, max_x, max_y, total_length) {
  add <- total_length / num_parts_axis
  for (vy in 0:num_parts_axis) {
    abline(h = min_y + vy * add, col = "red")
  }
  for (hx in 0:num_parts_axis) {
    abline(v = min_x + hx * add, col = "red")
  }
}


plot_graph <-function(num_parts_axis, min_x, min_y, max_x, max_y, total_length, data, c_x, c_y){
  plot(data$X, data$Y, cex = 0.01, xlim=c(min_x,max_x),ylim=c(min_y,max_y))
  plot_regions_from_whole(num_parts_axis = num_parts_axis, min_x = min_x, min_y = min_y, max_x = max_x, max_y = max_y, total_length = total_length)
  abline(v = c_x, col = "blue")
  abline(h = c_y, col = "blue")
}



To_xy_coor = function(RA, DEC, RA_gal = 204.2516509, 
                      DEC_gal = -29.86547, PA_gal = 54, 
                      i_gal = 12.50496329, D_gal = 4.89778819){
  RA = RA - RA_gal
  DEC = DEC - DEC_gal
  PA = PA_gal*pi/180
  i  = i_gal*pi/180
  # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
  x = RA*cos(PA) - DEC*sin(PA)
  y = DEC*cos(PA) + RA*sin(PA)
  # 2: Stretch x values to remove inclination effects
  x = x /cos(i)
  # 3: Convert deg to kpc
  x = x*pi/180*D_gal*1000
  y = y*pi/180*D_gal*1000
  xy_coor = cbind(x,y)
  return(xy_coor)
}



deprojected_distances = function(RA1, DEC1, RA2, DEC2){
  PA = deg2rad(54)
  i = deg2rad(15.3)
  # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
  x1 <- RA1 * cos(PA) - DEC1 * sin(PA)
  y1 <- DEC1 * cos(PA) + RA1 * sin(PA)
  x2 <- RA2 * cos(PA) - DEC2 * sin(PA)
  y2 <- DEC2 * cos(PA) + RA2 * sin(PA)
  # 2: Stretch x values to remove inclination effects
  long_x1 = x1/cos(i)
  long_x2 = x2/cos(i)
  # 3: Compute Euclidean Distances between x1,y1 and x2,y2 to get angular offsets (degrees).
  vec1 = cbind(y1, long_x1)
  vec2 = cbind(y2, long_x2)
  deg_dists = as.matrix(dist(rbind(vec1, vec2)))[1:nrow(vec1), (nrow(vec1)+1):(nrow(vec1)+nrow(vec2))]
  rad_dists = deg2rad(deg_dists)
  # 4: Convert angular offsets to kpc distances using D, and the small-angle approximation.
  Mpc_dists = rad_dists * 4.66
  kpc_dists = Mpc_dists * 1000
  return(kpc_dists)
}


deprojected_distances2 = function(RA1, DEC1, RA2, DEC2, PAd, inc){
  PA = deg2rad(PAd)
  i = deg2rad(inc)
  # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
  x1 <- RA1 * cos(PA) - DEC1 * sin(PA)
  y1 <- DEC1 * cos(PA) + RA1 * sin(PA)
  x2 <- RA2 * cos(PA) - DEC2 * sin(PA)
  y2 <- DEC2 * cos(PA) + RA2 * sin(PA)
  # 2: Stretch x values to remove inclination effects
  long_x1 = x1/cos(i)
  long_x2 = x2/cos(i)
  # 3: Compute Euclidean Distances between x1,y1 and x2,y2 to get angular offsets (degrees).
  vec1 = cbind(y1, long_x1)
  vec2 = cbind(y2, long_x2)
  deg_dists = as.matrix(dist(rbind(vec1, vec2)))[1:nrow(vec1), (nrow(vec1)+1):(nrow(vec1)+nrow(vec2))]
  rad_dists = deg2rad(deg_dists)
  # 4: Convert angular offsets to kpc distances using D, and the small-angle approximation.
  Mpc_dists = rad_dists * 4.66
  kpc_dists = Mpc_dists * 1000
  return(kpc_dists)
}
