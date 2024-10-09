library(ggplot2)
library('dplyr')
library('tidyr')
library(spdep)

set.seed(2024) 
###
options(StringsAsFactor = FALSE)
########################################################
########################################################
### Function to calculate CLQ (co-localization quotient)
get_CLQ <- function(cell_type_assignment,coords,nb_list,cell_A,cell_B,
                    number_of_neighbor,cell_B_name){
  ### Cell B is the cell type to check, cell A is attracted to B
  ### CLQ in this case is the degree to which A can be attached to B
  ### How much A is around B cell
  cell_A_indices <- which(cell_type_assignment == cell_A)
  cell_B_indices <- which(cell_type_assignment == cell_B)
  
  cell_A_nb_of_B <- find_cell_type_neighbors(cell_B,nb_list,cell_type_assignment,
                                             cell_A_indices,cell_B_name)
  
  total_cell <- length(which(cell_type_assignment!=0 & is.na(cell_type_assignment)==FALSE))
  Cab <- 0
  for(i in 1:length(cell_A_indices)){
    Cab <- Cab + length(cell_A_nb_of_B[i,1][[1]])
  }
  if(cell_A==cell_B){  #ie. CLQ of t cell to t cell
    if(length(cell_B_indices)<=5 | length(cell_A_indices)<=5){ # min is 5 cells per pop'n
      CLQ_ab <- 0
    }else{
      CLQ_ab <- (Cab/length(cell_A_indices))/((length(cell_B_indices)-1)/(total_cell-1))
    }
  }else{ #ie. CLQ of t cell to fibroblast 
    if(length(cell_B_indices)<=5 | length(cell_A_indices)<=5){ 
      CLQ_ab <- 0
    }else{
      CLQ_ab <- (Cab/length(cell_A_indices))/(length(cell_B_indices)/(total_cell-1))
    }
  }
  return(CLQ_ab)
}
########################################################################
########################################################################
### Function to find the cell types for neighboring cells
find_cell_type_neighbors <- function(cell_B,nb_list,cell_type_assignment,
                                     cell_A_indices,cell_B_name){
  if(length(cell_A_indices)==0){
    
  }else{
    cell_type_nb <- matrix(rep(list(),length(cell_A_indices)),nrow=length(cell_A_indices),ncol=1)
    row.names(cell_type_nb) <- cell_A_indices
    colnames(cell_type_nb) <- cell_B_name
    for(j in 1:length(cell_A_indices)){
      #print(j)
      current_cell_ID <- cell_A_indices[j]
      neighbors <- nb_list[[current_cell_ID]]
      neighbor_types <- cell_type_assignment[neighbors]
      
      cell_type_neighbor_IDs <- neighbors[which(neighbor_types == cell_B)]
      cell_type_nb[j,1][[1]] <- cell_type_neighbor_IDs
    }
    return(cell_type_nb)
  }
}
########################################################
### Function to find N-nearest neighboring cells
KNN_neighbors <- function(coords,number_of_neighbors){
  xxx <- knearneigh(coords,k=number_of_neighbors)
  nb_list <- list()
  for(i in 1:dim(xxx$nn)[1]){
    nb_list[[i]] <- xxx$nn[i,]
  }
  return(nb_list)
}



######################################################

##functions_for_PA_analysis_2_permutated_CLQ_values_modified}

#this part is taking the cell assignment information from the table.
#needs to be sampled with replacement=FALSE (ie. N (permutation number) = X times depend on the input from the user. )
#then CLQ function will be recalled with new assignments.


CLQ_permutated_matrix_gen <- function(iternum, filename, df_c){
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  
  N_PA = iternum # permutation iteration number
  
  #### CELL_TYPES 
  prior_info_cell_types_list <- "cell_types_celesta.csv"
  prior_info <- read.csv(prior_info_cell_types_list,header = TRUE,check.names = FALSE)
  cell_type_num <- seq(1,dim(prior_info)[1],by=1)
  cell_types <- prior_info[,1]
  
  ### Read in cell type assignment file
  cell_type_assignment_file <- read.csv(filename,header=TRUE,check.names = FALSE)
  cell_type_assignment <- cell_type_assignment_file$`Cell type number`
  
  ### get the X and Y coordinates
  coords <- cbind(cell_type_assignment_file$X,cell_type_assignment_file$Y)
  colnames(coords) <- c("X","Y")
  
  ### Build a list with each row is a cell and each columns contain the index for the N-nearest neighbors
  n_neighbors <- 20 # Using 20-nearest neighbors
  nb_list <- KNN_neighbors(coords,number_of_neighbors=n_neighbors)
  
  
  ##################################################################
  
  ### Creating the template dataframe with all the cell-cell combinations. 
  CLQ_matrix_R <- data.frame(matrix(0L, nrow = 0, ncol = 2 ))
  colnames(CLQ_matrix_R) <- c('cellA','cellB')
  for (cA in 1:length(cell_types)){
    for (cB in 1:length(cell_types)){
      CLQ_matrix_R[nrow(CLQ_matrix_R) + 1,] = c(cA, cB)
    }}
  
  ### Generating permutated CLQs
  
  for (r in 1:N_PA){
    
    ### Generating randomized data. Spatial coordinates do not change, only the cell type assignments are randomized, with replacement FALSE, so the cell counts are not changing. 
    cell_type_assignment_randomized <- sample(cell_type_assignment, length(cell_type_assignment), replace = FALSE)
    
    CLQ_array = c()
    ### k is Cell type A
    for(k in 1:length(cell_type_num)){
      count_a = df_c[k, 'Freq']
      ### if the current cell type do not exist in the sample
      if (is.na(count_a)){
        CLQ_array = append(CLQ_array, rep(0,18))
        next
      }
      ### j is Cell type A
      for(j in 1:length(cell_type_num)){
        cell_A <- cell_type_num[k]
        cell_B <- cell_type_num[j]
        cell_B_name <- cell_types[j]
        CLQ_result_R <- get_CLQ(cell_type_assignment_randomized,coords,nb_list,cell_A,cell_B,
                                number_of_neighbor=20,cell_B_name)
        ### Adding CLQ to the array, 
        CLQ_array = append(CLQ_array, CLQ_result_R)
      }
    }
    ### CLQ values for all combinations of the current iteration are added to the dataframe under the columnn P"iteration number", ie P1, P2 ... 
    CLQ_matrix_R[,paste0('P',as.character(r))] <- CLQ_array
  }
  
  ### Save CLQ results file with the sample name
  filename <- paste0(sample_to_check,"_CLQ_Permutated.csv")
  write.csv(CLQ_matrix_R,filename,  row.names = FALSE)
  
  filename <- paste0(sample_to_check,"_CLQ_Permutated.rds")
  saveRDS(CLQ_matrix_R, filename)
  
  
  #############
}

######################################################

get_counts = function(filename){
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  
  prior_info_cell_types_list <- "cell_types_celesta.csv"
  prior_info <- read.csv(prior_info_cell_types_list,header = TRUE,check.names = FALSE)
  cell_types <- prior_info[,1]
  cell_types = cell_types[1:18]
  
  ### Read in cell type assignment file
  cell_type_assignment_file <- read.csv(filename,header=TRUE,check.names = FALSE)
  cell_type_assignment <- cell_type_assignment_file$`Cell type number`
  
  ### table function in R gives the occurences of each variable in that column. so it will be cell type 1  : 21, cell type 2 : 432 etc. 
  counts_original = table(cell_type_assignment)
  
  df_a = data.frame( cell_type = cell_types, cell_type_assignment = c(1:18)) ### cell type number and the names
  df_b<- data.frame(counts_original)
  df_b$cell_type_assignment = as.numeric(as.character(df_b$cell_type_assignment))   ### cell type names and the counts
  df_c= df_a %>% dplyr::left_join( df_b, 
                                   by=c('cell_type_assignment'))
  
  
  filename <- paste0(sample_to_check,"_CellCounts.csv")
  write.csv(df_c,filename,  row.names = FALSE)
  
  return(df_c)
}


######################################################

#File names will be found by the extension. 
files <- (Sys.glob("*cell_type_assignment.csv"))

for (f in files){
  print(f)
  filename_c = f
  
  count_file = get_counts(filename=filename_c)
  
  
  CLQ_permutated = CLQ_permutated_matrix_gen(iternum=500, 
                                             filename = filename_c,
                                             df_c = count_file)
}

### Read in CELESTA cell type assignments

#input is f: the full name of the file. 
#this function will read the file, and will generate the original CLQ matrix, and it will save it as paste0(sample_to_check,"_CLQ.csv")

CLQ_matrix_gen <- function(filename) {
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  
  prior_info_cell_types_list <- "cell_types_celesta.csv"
  
  cell_type_assignment_file <- read.csv(filename,header=TRUE,check.names = FALSE)
  
  cell_type_assignment <- cell_type_assignment_file$`Cell type number`
  
  
  ### get the X and Y coordinates
  coords <- cbind(cell_type_assignment_file$X,cell_type_assignment_file$Y)
  colnames(coords) <- c("X","Y")
  
  ### Read in prior cell-type signature matrix
  
  prior_info <- read.csv(prior_info_cell_types_list,header = TRUE,check.names = FALSE)
  
  ### Build a list with each row is a cell and each columns contain the index for the N-nearest neighbors
  n_neighbors <- 20  #10 # Using 10-nearest neighbors
  nb_list <- KNN_neighbors(coords,number_of_neighbors=n_neighbors)
  
  ### Generate a sequence from 1 to the number of cell types. For example, if there are 26 cell types
  ### So it is a sequence from 1,2,3,...,26
  cell_type_num <- seq(1,dim(prior_info)[1],by=1)
  ### Cell type names from the prior matrix with the same order
  cell_types <- prior_info[,1]
  
  ### Build a matrix with CLQ value for one sample. Each row is a cell type A and each column is a cell type B
  ### Please refer to CELESTA paper for more details on the formula
  CLQ_matrix <- matrix(0L, nrow = length(cell_type_num), ncol = length(cell_type_num))
  row.names(CLQ_matrix) <- cell_types
  colnames(CLQ_matrix) <- cell_types
  
  for(k in 1:length(cell_type_num)){
    for(j in 1:length(cell_type_num)){
      cell_A <- cell_type_num[k]
      cell_B <- cell_type_num[j]
      
      cell_B_name <- cell_types[j]
      CLQ_result <- get_CLQ(cell_type_assignment,coords,nb_list,cell_A,cell_B,
                            number_of_neighbor= 20,cell_B_name)
      CLQ_matrix[k,j] <- CLQ_result
      #CLQmax_matrix[k,j] <- CLQ_result[[1]]
    }
  }
  
  ### This is the original CLQ file. (with real, unmodified data)
  ### Save CLQ results file with the sample name
  filename <- paste0(sample_to_check,"_CLQ.csv")
  write.csv(CLQ_matrix,filename)
  
  
  filename <- paste0(sample_to_check,"_CLQ.rds")
  saveRDS(CLQ_matrix, filename)
  
  
  return(CLQ_matrix)
}


#CLQ_matrix_gen("TAFs1_Erlot_cell_type_assignment.csv")



### Reading in the permutated CLQ file into workspace. 

# inputs: file_name

CLQ_permutated_matrix_gen <- function( filename){
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  filename <- paste0(sample_to_check,"_CLQ_Permutated.csv")
  CLQ_Permutated_file <- read.csv(filename,header=TRUE,check.names = FALSE)
  
  return(CLQ_Permutated_file)
}


#input: iternum (iteration number),  filename, list_of_matrices (permutated CLQ values), CLQ_matrix (original CLQ values), df_c is count table. 

significance_matrix_gen = function(iternum,
                                   filename,
                                   list_of_matrices,
                                   CLQ_matrix_original,
                                   df_c){
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  
  prior_info_cell_types_list <- "cell_types_celesta.csv"
  prior_info <- read.csv(prior_info_cell_types_list,header = TRUE,check.names = FALSE)
  ### Generate a sequence from 1 to the number of cell types. In Irene's data, there are 26 cell types
  ### So it is a sequence from 1,2,3,...,26
  cell_type_num <- seq(1,dim(prior_info)[1],by=1)
  ### Cell type names from the prior matrix with the same order
  cell_types <- prior_info[,1]
  
  
  ### Creating a dataframe for the output.
  df_full <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(df_full) <-c('sample','Cell_A' ,'Cell_A_num','Cell_A_count',
                        'Cell_B','Cell_B_num','Cell_B_count ',
                        'original_CLQ',
                        'percentile')
  
  ### col_array: is the column array consisting of P1, P2 up to P500 (for iteration number 500), It is used to read all the permutated CLQ values from the dataframe, 
  ### and these values are put into an array, called "CLQ_array".
  col_array = c()
  for (i in 1:iternum) { 
    col_array = append(col_array, paste0("P",as.character(i)))
  }
  
  for (cA in 1:length(cell_types)){
    for (cB in 1:length(cell_types)){
      CLQ_array = c()
      
      ### current_row is the row number calculated by the cell types number, ie. row number corresponding to the CLQ of cell type 2 and cell type 5 should be: (2-1)*18 + 5 = 23. 
      current_row = (cA-1)*length(cell_types)+cB
      CLQ_array = unname(unlist(list_of_matrices[current_row,col_array]))
      count_a = df_c[cA, 'Freq']
      count_b = df_c[cB, 'Freq']
      cell_type_A = cell_types[cA]
      cell_type_B = cell_types[cB]  
      
      ### Extracting the original value.
      original = CLQ_matrix_original[cell_type_A,cell_type_B]
      
      
      ### Checking percentile. If original value is on the far right of the permutated values, its rank is 1. (if far-left, it is 0)
      if (original> max(CLQ_array)) {
        CLQ_matrix_Rank = 1 
      } else if (original< min(CLQ_array)) {
        CLQ_matrix_Rank = 0 
      } else {
        ### Percentile calculation: Appending original value to the permutated CLQ array, then sorting, taking the index of the original CLQ and dividing by the total sample size.
        CLQ_matrix_Rank = which(sort(append(CLQ_array, original)) == original)[1]/iternum
      }
      df_full[nrow(df_full) + 1,] = c(sample_to_check, cell_type_A, cA ,count_a,  cell_type_B, cB , count_b, original, CLQ_matrix_Rank)
      ## CLQ_matrix_Rank goes to "percentile" column.
    }}
  
  ### significance : only the left or right 5% percentile. 
  df_full = df_full %>% dplyr::mutate(significant = if_else(percentile <= 0.05 | percentile >= 0.95, TRUE, FALSE))
  ### Save CLQ PA analysis results 
  
  filename <- paste0(sample_to_check,"_CLQ_data_full.csv")
  write.csv(df_full,filename)
  
  return(df_full)
}


get_counts = function(filename){
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  filename <- paste0(sample_to_check,"_CellCounts.csv")
  df_c <- read.csv(filename,header=TRUE,check.names = FALSE)
  
  return(df_c)
}





# for loop: every combination.  

plot_gen = function(iternum, filename, list_of_matrices , CLQ_matrix_original){
  
  col_array = c()
  for (i in 1:iternum) 
  { col_array = append(col_array, paste0("P",as.character(i)))}
  
  sample_to_check <- substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv")) 
  
  
  path = file.path(paste0("./PA_figures_",  sample_to_check ))
  dir.create(path)
  
  #### CELL_TYPES 
  prior_info_cell_types_list <- "cell_types_celesta.csv"
  prior_info <- read.csv(prior_info_cell_types_list,header = TRUE,check.names = FALSE)
  cell_types <- prior_info[,1]
  cell_types = cell_types[1:18]
  
  for (cA in 1:18){
    for (cB in 1:18){
      
      CLQ_array = c()
      current_row = (cA-1)*length(cell_types)+cB
      CLQ_array = unname(unlist(list_of_matrices[current_row,col_array]))
      cell_type_A = cell_types[cA]
      cell_type_B = cell_types[cB]  
      
      original = CLQ_matrix_original[cell_type_A,cell_type_B]
      original = CLQ_matrix[cell_type_A,cell_type_B]
      
      if (is.nan(min(original,CLQ_array))) {
        # print(original)
        next
      }
      png(file= file.path(paste0("./PA_figures_",  sample_to_check  ,"/cell_",cA,"_", cB,".png")),width=1000, height=600, res=200 )
      par(cex.axis=1.0)
      plot = hist(na.omit(CLQ_array), breaks=50, col='grey', main= paste0(cell_type_A,"\n", cell_type_B),
                  cex.main=0.7,
                  cex.lab=0.7,
                  cex.sub=0.7,
                  cex.axis=0.7, 
                  las=1, xlab='CLQ', 
                  xlim=c(min(original,CLQ_array)-5,
                         ifelse(original > max(CLQ_array),original+5,max(CLQ_array)+5)))
      abline(v=original, lwd=3, col="red")
      dev.off()
    }
  }
}


#################################################
#################################################
#################################################

files <- (Sys.glob("*cell_type_assignment.csv"))

for (f in files){
  print(f)
  
  filename_c = f
  
  CLQ_matrix = CLQ_matrix_gen(filename = filename_c)
  
  CLQ_permutated = CLQ_permutated_matrix_gen(filename = filename_c)
  
  count_file = get_counts(filename=filename_c)
  
  significance_matrices = significance_matrix_gen(iternum=500,
                                                  filename = filename_c,
                                                  list_of_matrices = CLQ_permutated,
                                                  CLQ_matrix_original= CLQ_matrix,
                                                  df_c = count_file)
  
  plot_gen(iternum=500,
           filename = filename_c,
           list_of_matrices = CLQ_permutated,
           CLQ_matrix_original= CLQ_matrix)
  
  
}



#################################################
#################################################
#################################################

CLQ_normalization_by_sample <- function(CLQ_before_normalization,cell_count,number_of_nearest_neighbors,
                                        threshold_rare_population=5,prior_matrix,
                                        clipping_right_tail=0.05, clipping_left_tail=0,
                                        plot_distribution=FALSE){
  ### This function requires (1) a named vector with the original CLQ values for one sample 
  ### before normalization, each element need to have a name, which is the two cell types 
  ### in the cell pair, connected by "_"
  ### (2) A cell count file with the number of cells for each cell type in that sample
  ### (3) Number of nearest neighbors in the CLQ calculation
  ### (4) A threshold value cell count for rare cell populations, default is 5
  ### (5) CELESTA input prior cell type signature matrix
  ### (6) Clipping parameters
  
  ### Obtain the cell pair names
  cell_pairs <- names(CLQ_before_normalization)
  for(i in 1:length(CLQ_before_normalization)){
    cell_pair_to_check <- strsplit(cell_pairs[i],"_")
    cell_A <- cell_pair_to_check[[1]][1]
    cell_B <- cell_pair_to_check[[1]][2]
    cell_B_count <- cell_count[which(prior_matrix[,1]==cell_B)]
    cell_A_count <- cell_count[which(prior_matrix[,1]==cell_A)]
    if(cell_B_count<=threshold_rare_population | cell_A_count<=threshold_rare_population){
      CLQ_before_normalization[i] <- NA
    }else{}
  }
  
  if(plot_distribution==TRUE){
    plot(density(CLQ_before_normalization[is.na(CLQ_before_normalization)==FALSE]),col="magenta",
         main="Original CLQs",xlab="CLQs before normalizaiton")
  }
  
  CLQ_before_normalization_no_NA <-  CLQ_before_normalization[which(is.na(CLQ_before_normalization)==FALSE)] 
  ### Order the original CLQ distritbution
  
  order_distribution <- CLQ_before_normalization_no_NA[order(CLQ_before_normalization_no_NA)]
  ### Clipping the original distribution
  max_CLQ <- order_distribution[floor(length(CLQ_before_normalization_no_NA)*(1-clipping_right_tail))]
  if(clipping_left_tail==0){
    min_CLQ <- min(order_distribution)
  }else{min_CLQ <- order_distribution[floor(length(CLQ_before_normalization_no_NA)*clipping_left_tail)]}
  
  original_CLQ_clipped <- CLQ_before_normalization
  for(i in 1:length(CLQ_before_normalization)){
    ##II: if original CLQ > max and it is not NA: we set the CLQ same as Max CLQ.
    if(CLQ_before_normalization[i]>=max_CLQ&is.na(CLQ_before_normalization[i])==FALSE){
      original_CLQ_clipped[i] <- max_CLQ
    }##II: if original CLQ > max and it is not NA: we set the CLQ same as Max CLQ.
    else if(CLQ_before_normalization[i]<=min_CLQ&is.na(CLQ_before_normalization[i])==FALSE){
      original_CLQ_clipped[i] <- min_CLQ
    }
  }
  if(plot_distribution==TRUE){
    plot(density(original_CLQ_clipped[is.na(original_CLQ_clipped)==FALSE]),col="blue",
         main="Original CLQs clipped",xlab="CLQs")
  }
  
  mean_distribution <- mean(original_CLQ_clipped[is.na(original_CLQ_clipped)==FALSE])
  if(abs(mean_distribution-number_of_nearest_neighbors)>3){
    print("Please assess the distribution and consider clipping more")
  }
  sd_distribution <- sd(original_CLQ_clipped[is.na(original_CLQ_clipped)==FALSE])
  ### Z score
  normalized_CLQ <- (original_CLQ_clipped-mean_distribution)/sd_distribution
  
  if(plot_distribution==TRUE){
    plot(density(normalized_CLQ[is.na(original_CLQ_clipped)==FALSE]),col="orange",
         main="Normalized CLQ Z score",xlab="Z scores")
  }
  return(normalized_CLQ)
}
#######################################################################
#######################################################################

call_normalization_fnc = function(filename,  righttail , lefttail){
  sample_to_check = substr(filename, start = 1, stop = nchar(filename) - nchar("_cell_type_assignment.csv"))
  cell_type_assignment_file <- read.csv(filename,header=TRUE,check.names = FALSE)
  
  filename <- paste0(sample_to_check,"_CellCounts.csv")
  count_file <- read.csv(filename,header=TRUE,check.names = FALSE)
  cell_count = count_file$Freq
  
  filename <- paste0(sample_to_check,"_CLQ_data_full.csv")
  significance_matrices <- read.csv(filename,header=TRUE,check.names = FALSE)
  
  CLQ_before_normalization= as.numeric(significance_matrices$original_CLQ)
  
  names(CLQ_before_normalization) <- paste0(significance_matrices$Cell_A, "_", significance_matrices$Cell_B )
  
  prior_info_cell_types_list <- "cell_types_celesta.csv"
  prior_info <- read.csv(prior_info_cell_types_list,header = TRUE,check.names = FALSE)
  
  normalized_CLQ = CLQ_normalization_by_sample(CLQ_before_normalization,cell_count,
                                               number_of_nearest_neighbors = 20,
                                               threshold_rare_population=5,prior_matrix= prior_info,
                                               clipping_right_tail=righttail, clipping_left_tail=lefttail,
                                               plot_distribution=TRUE)
  
  df_normalized_CLQ <- data.frame(normalized_CLQ)
  
  ##
  rownames(df_normalized_CLQ) = 1:324
  df_normalized_CLQ$Cell_A =significance_matrices$Cell_A
  df_normalized_CLQ$Cell_B = significance_matrices$Cell_B 
  df_normalized_CLQ = df_normalized_CLQ[,c(2,3,1)]
  
  ##
  df_normalized_CLQ = df_normalized_CLQ %>%
    dplyr::rename({{sample_to_check}} := "normalized_CLQ" )
  
  #save normalized CLQ values
  filename <- paste0(sample_to_check,"_CLQ_Normalized", "_L",lefttail, "_R", righttail)
  ##
  write.csv(df_normalized_CLQ, paste0(filename,".csv"),  row.names = TRUE)
  ##
  saveRDS(df_normalized_CLQ, paste0(filename,".rds"))
  
  return()
  
}
# 
# #calling the function
# #Example 1 : assembloid sample
# filename = "TAFs1_cell_type_assignment.csv"
# righttail = 0.05
# lefttail = 0
# call_normalization_fnc(filename , righttail , lefttail)
# # 
# 
#calling the function
#Example 2 : clinical specimen
# filename = "lepidic_58.6_cell_type_assignment.csv"
# righttail = 0.05
# lefttail = 0
# call_normalization_fnc(filename , righttail , lefttail)

# filename = "lepidic_58.6_cell_type_assignment.csv"
# righttail = 0.15
# lefttail = 0
# call_normalization_fnc(filename , righttail , lefttail)
# 

files <- (Sys.glob("*cell_type_assignment.csv"))

for (f in files){
  print(f)
  filename=f
  righttail = 0.05
  lefttail = 0
  call_normalization_fnc(filename , righttail , lefttail)
  
}

