#' @name venn_diagram4
#' @description by Niels Hanson
#'  Info: A function to create a venn digagram between three (3) overapping sets a, b, and c.
#'  Requies: The VennDiagram R package, the program will try to download and install if you dont have it.
#'  Inputs: a,b,c - vertical vectors or the first column of matrixes as sets of strings to be compared name_a, name_b, 
#'  and name_c - names that you give to sets a, b, and c, respectively these will appear in the final image 
#'  colors - a vector of colors defined in hex format ie. colors <- c("#B80830","#EACC33","#46E2D9")
#'  it might be a good idea that these mix well as the interlapping classes will be a mix of colors
#'  Colors must be specified in hex ie. "#123456", google "hex colors" for more info
#'  name_output - the name of the output image that will be put in R's current working directory
#'  you can set this yourself before running. ie. setwd("~/my-director/")
#'  A: find all classes for set a, a_b, a_c, a_b_c
#'  number in a, we will use this at the end for a sanity check
#'  make sure that a,b,c are all sets, get rid of duplicates
#' @param a unique element of a set
#' @param b unique element of a set
#' @param c unique element of a set
#' @param d unique element of a set
#' @param name_a
#' @param name_b
#' @param name_c
#' @param name_d
#' @param colors The colors for the four set. Default are c("#542BFD","#FCFF4F","#5CFF4A", "#F63C4B").
#' @param out logical. Whether to save to a file. 
#' @param file_name The name of file to save the plot. 
#' @param euler
#' @param special
#' @export
venn_diagram4 <- function(a, b, c, d, name_a = "A", name_b = "B", name_c = "C", name_d = "D", colors =  c("#542BFD","#FCFF4F","#5CFF4A", "#F63C4B"), out = FALSE, file_name = "default.pdf", euler=FALSE, special=TRUE){
  
  a <- intersect(a,a)
  b <- intersect(b,b)
  c <- intersect(c,c)
  d <- intersect(d,d)
  # length of a
  a_tot <- length(a)
  b_tot <- length(b)
  c_tot <- length(c)
  d_tot <- length(d)
  
  # number shared between a and b and c and d
  a_b_c_d <- intersect(intersect(a,b),intersect(c,d))
  a_b_c_d_tot <- length(a_b_c_d)
  a_b_c_d_tot
  
  # a, b, and c, check to make sure not in d
  a_b_c <- setdiff(intersect(intersect(a,b), c), d)
  a_b_c_tot <- length(a_b_c)
  
  # b,c, and d, checking to make sure not in a
  b_c_d <-setdiff(intersect(intersect(b,c), d), a)
  b_c_d_tot <- length(b_c_d)
  
  # a,c, and d, checking to make sure not in b
  a_c_d <- setdiff(intersect(intersect(a,c), d),b)
  a_c_d_tot <- length(a_c_d)
  
  # a,b, and d, checking to make sure not in c
  a_b_d <- setdiff(intersect(intersect(a,b),d),c)
  a_b_d_tot <- length(a_b_d)
  
  # all pair sets only
  a_b <- setdiff(intersect(a,b),union(c,d))
  a_b_tot <- length(a_b)
  a_c <- setdiff(intersect(a,c),union(b,d))
  a_c_tot <- length(a_c)
  a_d <- setdiff(intersect(a,d),union(b,c))
  a_d_tot <- length(a_d)
  b_c <- setdiff(intersect(b,c),union(a,d))
  b_c_tot <- length(b_c)
  b_d <- setdiff(intersect(b,d),union(a,c))
  b_d_tot <- length(b_d)
  c_d <- setdiff(intersect(c,d),union(a,b))
  c_d_tot <- length(c_d)
  
  # a, b, c, d only
  a_only <- setdiff(a,union(union(b,c),d))
  a_only_tot <- length(a_only)
  b_only <- setdiff(b,union(union(a,c),d))
  b_only_tot <- length(b_only)
  c_only <- setdiff(c,union(union(a,b),d))
  c_only_tot <- length(c_only)
  d_only <- setdiff(d,union(union(a,b),c))
  d_only_tot <- length(d_only)
  
  # sanity checks
  a_tot == (a_only_tot + a_b_tot + a_b_c_tot + a_b_c_d_tot + a_b_d_tot + a_d_tot + a_c_d_tot + a_c_tot)
  b_tot == (b_only_tot + b_c_tot + b_c_d_tot + b_d_tot + a_b_d_tot + a_b_c_d_tot + a_b_c_tot + a_b_tot)
  c_tot == (c_only_tot + c_d_tot + b_c_d_tot + a_b_c_d_tot + a_c_d_tot + a_c_tot + a_b_c_tot + b_c_tot)
  d_tot == (d_only_tot + b_d_tot + a_b_d_tot + a_d_tot + a_c_d_tot + a_b_c_d_tot + b_c_d_tot + c_d_tot)
  
  # try to load VennDiagram else try to download it
  # More info: see Hanbo Chen and Paul C Boutros. BMC Bioinformatics 2011, 12:35 doi:10.1186/1471-2105-12-35
  try(library(VennDiagram), install.packages("VennDiagram")) 
  library(VennDiagram)
  
  
  # create the first list
  offset_a <- 1*(10^6)
  offset_b <- 2*(10^6)
  offset_c <- 3*(10^6)
  offset_d <- 4*(10^6)
  offset_a_b <- 5*(10^6)
  offset_a_c <- 6*(10^6)
  offset_a_d <- 7*(10^6)
  offset_b_c <- 8*(10^6)
  offset_b_d <- 9*(10^6)
  offset_c_d <- 10*(10^6)
  offset_a_b_c <- 11*(10^6)
  offset_b_c_d <- 12*(10^6)
  offset_a_c_d <- 13*(10^6)
  offset_a_b_d <- 14*(10^6)
  offset_a_b_c_d <- 15*(10^6)
  
  list_a_only <- NULL
  list_b_only <- NULL
  list_c_only <- NULL
  list_d_only <- NULL
  list_a_b <- NULL
  list_a_c <- NULL
  list_a_d <- NULL
  list_b_c <- NULL
  list_b_d <- NULL
  list_c_d <- NULL
  list_a_b_c <- NULL
  list_b_c_d <- NULL
  list_a_c_d <- NULL
  list_a_b_d <- NULL
  list_a_b_c_d <- NULL
  
  if(length(a_only) > 0){
    list_a_only <- offset_a:(offset_a+length(a_only)-1);
  }
  if(length(b_only) > 0){
    list_b_only <- offset_b:(offset_b+length(b_only)-1);
  }
  if(length(c_only) > 0){
    list_c_only <- offset_c:(offset_c+length(c_only)-1);
  }
  if(length(d_only) > 0){
    list_d_only <- offset_d:(offset_d+length(d_only)-1);
  }			
  if(length(a_b) > 0){
    list_a_b <- offset_a_b:(offset_a_b+length(a_b)-1);
  }	
  if(length(a_c) > 0){
    list_a_c <- offset_a_c:(offset_a_c+length(a_c)-1);
  }
  if(length(a_d) > 0){
    list_a_d <- offset_a_d:(offset_a_d+length(a_d)-1);
  }
  if(length(b_c) > 0){
    list_b_c <- offset_b_c:(offset_b_c+length(b_c)-1);
  }
  if(length(b_d) > 0){
    list_b_d <- offset_b_d:(offset_b_d+length(b_d)-1);
  }	
  if(length(c_d) > 0){
    list_c_d <- offset_c_d:(offset_c_d+length(c_d)-1);
  }
  if(length(a_b_c) > 0){
    list_a_b_c <- offset_a_b_c:(offset_a_b_c+length(a_b_c)-1);
  }
  if(length(b_c_d) > 0){
    list_b_c_d <- offset_b_c_d:(offset_b_c_d+length(b_c_d)-1);
  }
  if(length(a_c_d) > 0){
    list_a_c_d <- offset_a_c_d:(offset_a_c_d+length(a_c_d)-1);
  }
  if(length(a_b_d) > 0){
    list_a_b_d <- offset_a_b_d:(offset_a_b_d+length(a_b_d)-1);
  }
  if(length(a_b_c_d) > 0){
    list_a_b_c_d <- offset_a_b_c_d:(offset_a_b_c_d+length(a_b_c_d)-1);
  }
  
  temp <- list(
    name_a = c(list_a_only, 
               list_a_b,
               list_a_b_c,
               list_a_c,
               list_a_c_d,
               list_a_d,
               list_a_b_c_d,				
               list_a_b_d
    ),
    name_d = c(list_d_only, 
               list_c_d, 
               list_b_c_d, 
               list_b_d,
               list_a_b_c_d,
               list_a_b_d,
               list_a_c_d,
               list_a_d
    ),
    name_b = c(list_a_b, 
               list_a_b_c,
               list_a_b_c_d,
               list_a_b_d,
               list_b_only,
               list_b_c,
               list_b_c_d,
               list_b_d
    ),
    name_c = c(list_c_only, 
               list_b_c,
               list_a_b_c,
               list_a_c,
               list_a_c_d,
               list_a_b_c_d,
               list_b_c_d,
               list_c_d
    )		
  )
  
  
  names(temp) <- c(name_a,name_d,name_b,name_c)
  
  output<-venn.diagram(
    x = temp,
    sp.cases = special,
    filename = NULL,
    col = "black",
    fill = c(colors[1], colors[4], colors[2], colors[3]),
    alpha = 0.5,
    label.col = c("black"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.default.pos = "text",
    cat.col = c("black"),
    cat.cex = 1.25,
    cat.fontfamily = "serif",
    cat.pos = 0,
    euler.d = euler,
    scaled = euler
  );
  
  try(library("grid"), install.packages("grid")) 
  library(grid)
  if (out == TRUE) {
    png(filename = paste(name_output,'.png', sep=''))
    grid.draw(output)
    dev.off()       
    pdf(file = paste(name_output,'.pdf', sep=''), width=8, height=8)
    grid.draw(output)
    dev.off()
  } else {
    # just plot the graph
    grid.draw(output)
  }
  
  # pack all the sets up and return to the user
  out <- NULL;
  out <- list(a,a_only,a_b,a_c,a_d,
              b,b_only,b_c,b_d,
              c,c_only,c_d,
              d,d_only,
              a_b_c,b_c_d,a_c_d,a_b_d,
              a_b_c_d);
  
  names(out) <- c(name_a,paste(name_a,"only",sep="_"),
                  paste(name_a,name_b,sep="_"),
                  paste(name_a,name_c,sep="_"),
                  paste(name_a,name_d,sep="_"),
                  name_b,
                  paste(name_b,"only",sep="_"),
                  paste(name_b,name_c,sep="_"),
                  paste(name_b,name_d,sep="_"),
                  name_c,
                  paste(name_c,"only",sep="_"),
                  paste(name_c,name_d,sep="_"),
                  name_d,
                  paste(name_d,"only",sep="_"),
                  paste(name_a,name_b,name_c,sep="_"),
                  paste(name_b,name_c,name_d,sep="_"),
                  paste(name_a,name_c,name_d,sep="_"),
                  paste(name_a,name_b,name_d,sep="_"),
                  paste(name_a,name_b,name_c,name_d,sep="_"))
}

#' @name venn_diagram3
#' @description by Niels Hanson
#'  Info: A function to create a venn digagram between three (3) overapping sets a, b, and c.
#'  Requies: The VennDiagram R package, the program will try to download and install if you dont have it.
#'  Inputs: a,b,c - vertical vectors or the first column of matrixes as sets of strings to be compared name_a, name_b, 
#'  and name_c - names that you give to sets a, b, and c, respectively these will appear in the final image 
#'  colors - a vector of colors defined in hex format ie. colors <- c("#B80830","#EACC33","#46E2D9")
#'  it might be a good idea that these mix well as the interlapping classes will be a mix of colors
#'  Colors must be specified in hex ie. "#123456", google "hex colors" for more info
#'  name_output - the name of the output image that will be put in R's current working directory
#'  you can set this yourself before running. ie. setwd("~/my-director/")
#'  A: find all classes for set a, a_b, a_c, a_b_c
#'  number in a, we will use this at the end for a sanity check
#'  make sure that a,b,c are all sets, get rid of duplicates
#' @param a unique element of a set
#' @param b unique element of a set
#' @param c unique element of a set
#' @param name_a
#' @param name_b
#' @param name_c
#' @param colors The colors for the four set. Default are c("#542BFD","#FCFF4F","#5CFF4A", "#F63C4B").
#' @param out logical. Whether to save to a file. 
#' @param file_name The name of file to save the plot. 
#' @param euler
#' @param special
#' @export
venn_diagram3 <- function(a, b, c, name_a = "A", name_b = "B", name_c = "C", colors =  c("#B80830","#EACC33","#46E2D9"), out = FALSE, file_name = "default.pdf", euler=FALSE, special=TRUE){
  # Niels Hanson
  # Info: A function to create a venn digagram between three (3) overapping sets a, b, and c.
  # Requies: The VennDiagram R package, the program will try to download and install if you dont
  #          have it.
  # Inputs: a,b,c - vertical vectors or the first column of matrixes as sets of strings to be compared
  #         name_a, name_b, and name_c - names that you give to sets a, b, and c, respectively
  #         these will appear in the final image
  #         colors - a vector of colors defined in hex format ie. colors <- c("#B80830","#EACC33","#46E2D9")
  #         it might be a good idea that these mix well as the interlapping classes will be a mix of colors
  #         Colors must be specified in hex ie. "#123456", google "hex colors" for more info
  #         name_output - the name of the output image that will be put in R's current working directory
  #         you can set this yourself before running. ie. setwd("~/my-director/")
  
  # A: find all classes for set a, a_b, a_c, a_b_c
  # number in a, we will use this at the end for a sanity check
  
  # make sure that a,b,c are all sets, get rid of duplicates
  a <- intersect(a,a)
  b <- intersect(b,b)
  c <- intersect(c,c)
  # length of a
  a_tot <- length(a);
  a_tot
  # number shared between a and b and c
  a_b_c <- intersect(a,intersect(b,c))
  a_b_c_tot <- length(a_b_c)
  a_b_c_tot
  # number shared with b
  a_b <- setdiff(intersect(a,b),a_b_c)
  a_b_tot <- length(a_b)
  a_b_tot
  # number shared with c
  a_c <- setdiff(intersect(a,c),a_b_c)
  a_c_tot <- length(a_c)
  a_c_tot
  # number in a only
  a_only <- setdiff(a, union(union(a_c,a_b), a_b_c))
  a_only_tot <- length(a_only)
  a_only_tot
  # we should have defined all classificaitons for set 'a'
  # sanity check should resolve TRUE if everything adds up
  a_tot == (a_b_c_tot + a_b_tot + a_c_tot + a_only_tot)
  
  # B: do the same for intersecting classes for b
  b_tot <- length(b)
  b_tot
  # number shared between b and c
  b_c <- setdiff(intersect(b,c),a_b_c)
  b_c_tot <- length(b_c)
  b_c_tot
  # b only
  b_only <- setdiff(b, union(union(b_c,a_b), a_b_c))
  b_only_tot <- length(b_only)
  b_only_tot
  
  # sanity check for B
  b_tot == (a_b_tot + a_b_c_tot + b_c_tot + b_only_tot)
  
  # C: do the same for intersecting classes for c
  c_tot <- length(c)
  c_tot
  c_only <- setdiff(c,union(union(b_c,a_c), a_b_c))
  c_only_tot <- length(c_only)
  c_only_tot
  
  # sanity check
  c_tot == (a_c_tot + b_c_tot + a_b_c_tot + c_only_tot)
  # try to load VennDiagram else try to download it
  # More info: see Hanbo Chen and Paul C Boutros. BMC Bioinformatics 2011, 12:35 doi:10.1186/1471-2105-12-35
  try(library(VennDiagram), install.packages("VennDiagram")) 
  library(VennDiagram)
  
  
  # create the first list
  offset_a <- 1*(10^6)
  offset_a_b <- 2*(10^6)
  offset_a_c <- 3*(10^6)
  offset_b <- 4*(10^6)
  offset_b_c <- 5*(10^6)
  offset_c <- 6*(10^6)
  offset_a_b_c <- 7*(10^6)
  
  list_a_only <- NULL
  list_a_b <- NULL
  list_a_c <- NULL
  list_b_only <- NULL
  list_b_c <- NULL
  list_c_only <- NULL
  list_a_b_c <- NULL
  
  if(length(a_only) > 0){
    list_a_only <- offset_a:(offset_a+length(a_only)-1);
  }
  if(length(a_c) > 0){
    list_a_c <- offset_a_c:(offset_a_c+length(a_c)-1);
  }	
  if(length(a_b) > 0){
    list_a_b <- offset_a_b:(offset_a_b+length(a_b)-1);
  }	
  if(length(b_only) > 0){
    list_b_only <- offset_b:(offset_b+length(b_only)-1);
  }	
  if(length(b_c) > 0){
    list_b_c <- offset_b_c:(offset_b_c+length(b_c)-1);
  }	
  if(length(c_only) > 0){
    list_c_only <- offset_c:(offset_c+length(c_only)-1);
  }
  if(length(a_b_c) > 0){
    list_a_b_c <- offset_a_b_c:(offset_a_b_c+length(a_b_c)-1)
  }
  
  
  temp <- list(
    name_a = c(list_a_only, 
               list_a_b,
               list_a_b_c,
               list_a_c),
    name_b = c(list_b_only, 
               list_a_b, 
               list_a_b_c, 
               list_b_c),
    name_c = c(list_c_only, 
               list_a_c,
               list_a_b_c,
               list_b_c
    )
  )
  
  # temp <- list(
  #           name_a = c(1:(length(a_only)), 
  #                     (length(a_only)):(length(a_only)-1+length(a_b)),
  #                     (length(a_only)+length(a_b)+1):(length(a_only)+length(a_b)+length(a_b_c)),
  #                       (length(a_only)+length(a_b)+length(a_b_c)+1):(length(a_only)+length(a_b)+length(a_b_c)+length(a_c))),
  #           name_b = c((a_tot+1):(a_tot+length(b_only)), 
  #                      (length(a_only)+1):(length(a_only)+length(a_b)), 
  #                      (length(a_only)+length(a_b)+1):(length(a_only)+length(a_b)+length(a_b_c)), 
  #                      (a_tot+length(b_only)+1):(a_tot+length(b_only)+length(b_c))),
  #           name_c = c((c_ind+1):(c_ind+1+length(c_only)), 
  #                      (length(a_only)+length(a_b)+1):(length(a_only)+length(a_b)+length(a_b_c)),
  #                      (length(a_only)+length(a_b)+length(a_b_c)+1):(length(a_only)+length(a_b)+length(a_b_c)+length(a_c)),
  #                      (a_tot+length(b_only)+1):(a_tot+length(b_only)+length(b_c)))
  #           )
  names(temp) <- c(name_a,name_b,name_c)
  
  output<-venn.diagram(
    x = temp,
    sp.cases = special,
    filename = NULL,
    col = "black",
    fill = c(colors[1], colors[2], colors[3]),
    alpha = 0.5,
    label.col = c("black", "white", "black", "white", "white", "white", "black"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.default.pos = "text",
    cat.col = c("black", "black", "black"),
    cat.cex = 1.5,
    cat.fontfamily = "serif",
    cat.dist = c(0.06, 0.06, 0.03),
    cat.pos = 0,
    euler.d = euler,
    scaled = euler
  );
  
  try(library("grid"), install.packages("grid")) 
  library(grid)
  if (out == TRUE) {
    # only produce ouptut files upon request
    png(filename = paste(file_name,'.png', sep=''))
    grid.draw(output)
    dev.off()
    pdf(file = file_name, width=8, height=8)
    grid.draw(output)
    dev.off()
  } else {
    # just plot the figure
    grid.draw(output)
  }
  
  
  out <- list(a,
              a_only,
              a_b,
              a_c,
              b,
              b_only,
              b_c,
              c,
              c_only,
              a_b_c)
  
  names(out) <- c(name_a,
                  paste(name_a,"only",sep="_"),
                  paste(name_a,name_b,sep="_"),
                  paste(name_a,name_c,sep="_"),
                  name_b,
                  paste(name_b,"only",sep="_"),
                  paste(name_b,name_c,sep="_"),
                  name_c,
                  paste(name_c,"only",sep="_"),
                  paste(name_a,name_b,name_c,sep="_"))
  
}
