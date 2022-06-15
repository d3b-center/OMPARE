# function to save network to file
network_to_file <- function(dnet_object, filename){
  ## Extract nodes from dtnetplot list object and remove shape parameter
  vertices <- dnet_object$x$nodes
  vertices$shape <- NULL
  
  ## Convert edge and node data frames to graph object
  net <- graph_from_data_frame(d = dnet_object$x$edges, vertices = vertices, directed = F)
  
  ## Assign layout
  l <- layout_with_fr(net)
  l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
  pdf(file = filename, width = 10, height = 10)
  plot(net, layout = l*2, vertex.size = 4, rescale = T, vertex.label.cex = 0.5)
  dev.off()
}
