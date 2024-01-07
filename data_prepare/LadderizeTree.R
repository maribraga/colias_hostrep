# Function to ladderize tree and correct tip label sequence

LadderizeTree <- function(tree, temp_file = "temp", orientation = "left"){
  if(file.exists(paste0("./", temp_file))){
    stop("The chosen temporary file exists! Please choose an other temp_file name")
  }
  if(orientation == "left"){
    right <- FALSE
  }else{
    right <- TRUE
  }
  tree_temp <- ladderize(tree, right = right)
  write.tree(tree_temp, file = paste0("./", temp_file, ".tre"))
  tree_lad <- read.tree(paste0("./", temp_file, ".tre"))
  file.remove(paste0("./", temp_file, ".tre"))
  return(tree_lad)
}
