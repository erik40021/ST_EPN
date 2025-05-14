get_spatial_aspect_ratio = function(sstobj) {
  coord <- GetTissueCoordinates(sstobj@images[[1]]) # grabs either 'tissue_lowres_image' or 'slice1'
  ratio = (max(coord[, 1]) - min(coord[, 1])) / (max(coord[, 2]) - min(coord[, 2]))
  return(ratio)
}

# define scale_factor as the size a spot should have in a full-size sample (spanning the full 5.5 x 5.5 mm area)
get_spatial_point_size = function(sstobj, scale_factor = 2) {
  coord <- GetTissueCoordinates(sstobj@images[[1]])
  area = (max(max(coord[, 1]) - min(coord[, 1]), max(coord[, 2]) - min(coord[, 2])))^2
  pt_size = max(scale_factor, scale_factor * 1/sqrt(area/2000000))
  return(pt_size)
}
