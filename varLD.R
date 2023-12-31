# start timer
start_time = proc.time()

## loading packages
source("/proj/yunligrp/users/dinelka/scripts/packages/libraries.R")

# input data
args = commandArgs(trailingOnly = TRUE) # CLI args
studies = c("BioMe", "JHS", "WHI")
chrom = args[1] # chromosome number
window_size = 1e6
n_eig_values = c(5, 10, 25, 50, 100)

# Read LD file
read_LD = function(cohort){
    fread(
        paste("/proj/yunligrp/users/dinelka/scripts/LD_calc_v2/", cohort,  "_", chrom, "_no_filter_0.1_1000000_results/", cohort, "_", chrom, "_no_filter_0.1_1000000_LD.txt.gz", sep = ""), header = TRUE, data.table = TRUE,
              select = c("SNP1", "SNP2", "R2", "+/-corr"), col.names = c("SNP1", "SNP2", "R2", "sign_corr")) 
}

# create matrix from DF
create_LD_matrix = function(df, n, df_name_1, df_name_2){
    sparseMatrix(i = df[[V1_ind]], j = df[[V2_ind]], x = (df[[paste('R2', df_name, sep="_")]] * sign_df1), dims = c(n,n), dimnames = list(vec_SNP,vec_SNP)) +
        sparseMatrix(i = df[[V2_ind]], j = df[[V1_ind]], x = (df[[paste('R2', df_name, sep="_")]] * sign_df1), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP)) +
        sparseMatrix(i = 1:n, j = 1:n, x = rep(1,n), dims = c(n,n), dimnames = list(vec_SNP,vec_SNP))
}

BioMe = read_LD("BioMe")
JHS = read_LD("JHS")
WHI = read_LD("WHI")

df_list = list(JHS, BioMe, WHI)

print(pryr::object_size(df_list))

n_plots = choose(length(df_list), 2)

# Main Function - computing eigen decomposition and plotting graphs

for (n_ev in n_eig_values){
  for (i in (1:(length(df_list)-1))){

    for (j in ((i + 1):length(df_list))){

      df_plot = data.frame(chrom = integer(), window = integer(), common_SNPs = integer(), varLD_score = double())
    
      df_name_1 = studies[i] # dataframe 1 name
      df_name_2 = studies[j] # dataframe 2 name

      ldAdmix_orig = inner_join(df_list[[i]], df_list[[j]], by = c("SNP1", "SNP2"), suffix = c(paste("_", df_name_1, sep=""), paste("_", df_name_2, sep=""))) ## need to check if this is valid!! Important

      size = pryr::object_size(ldAdmix_orig)
      
      print('Done inner joining')

        #plot location
      png(paste('/proj/yunligrp/users/dinelka/scripts/LD_calc/varLD/Figures/', df_name_1, '_', df_name_2, '/', chrom, '.png', sep = ""), width = 1100, height = 600)

      start = round_any(min(ldAdmix_orig$SNP1), window_size, f = ceiling) / window_size
      end = round_any(max(ldAdmix_orig$SNP1), window_size, f = ceiling) / window_size

      for (window in start:end) {

        ldAdmix = ldAdmix_orig[SNP1 <= (window*window_size) & SNP1 > ((window - 1)*window_size)] 

        print(head(ldAdmix))

        sign_df1 = ifelse(ldAdmix[[paste('sign_corr', df_name_1, sep="_")]] == "+", 1, -1)
        sign_df2 = ifelse(ldAdmix[[paste('sign_corr', df_name_2, sep="_")]] == "+", 1, -1)

        print('Done assigning signs')

        vec_SNP = c(unique(ldAdmix$SNP1), unique(ldAdmix$SNP2))

        n2 = length(vec_SNP)

        ldAdmix$V1_ind = match(ldAdmix$SNP1,vec_SNP)
        ldAdmix$V2_ind = match(ldAdmix$SNP2,vec_SNP)

        na_ind = apply(ldAdmix,1,anyNA) # perform function on rows
        ldAdmix = ldAdmix[!na_ind,] 

        n_ldAdmix = dim(ldAdmix)[1] #num of rows

        if (n_ldAdmix <= 1) {
          next
        }

        print('Starting to create the 2 LD matrices')

        ldmat_df1 = sparseMatrix(i = ldAdmix$V1_ind, j = ldAdmix$V2_ind, x = (ldAdmix[[paste('R2', df_name_1, sep="_")]] * sign_df1), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP)) +
          sparseMatrix(i = ldAdmix$V2_ind, j = ldAdmix$V1_ind, x = (ldAdmix[[paste('R2', df_name_1, sep="_")]] * sign_df1), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP)) +
          sparseMatrix(i = 1:n2, j = 1:n2, x = rep(1,n2), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP))

        ldmat_df2 = sparseMatrix(i = ldAdmix$V1_ind, j = ldAdmix$V2_ind, x = (ldAdmix[[paste('R2', df_name_2, sep="_")]] * sign_df2), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP)) +
          sparseMatrix(i = ldAdmix$V2_ind, j = ldAdmix$V1_ind, x = (ldAdmix[[paste('R2', df_name_2, sep="_")]] * sign_df2), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP)) +
          sparseMatrix(i = 1:n2, j = 1:n2, x = rep(1,n2), dims = c(n2,n2), dimnames = list(vec_SNP,vec_SNP))

        print('Done creating ldmat matrices')

        print(nrow(ldmat_df1))

        if (nrow(ldmat_df1) <= 2) {
          next
        }
        
        num_eigen = min(round(0.25 * nrow(ldmat_df1), digits = 0), n_ev) #using only minimum of the top 25% or top 1000 of eigen values 
        
        ev_df1 = eigs_sym(ldmat_df1, num_eigen, which = "LA")
        print('Done calculating df1 EVs')

        ev_df2 = eigs_sym(ldmat_df2, num_eigen, which = "LA")
        print('Done calculating df2 EVs')

        diffEV = ev_df1$values - ev_df2$values

        varLD_score = sum(abs(diffEV))

        print('Done calculating varLD score')

        df_plot = rbind(df_plot, data.frame(chrom = chrom, window = window, common_SNPs = n_ldAdmix, varLD_score = varLD_score))


        print(paste('Window ', window, ' of ', df_name_1, '_', df_name_2, sep = ""))

      }

      mean_varLD = mean(df_plot$varLD_score)
      var_varLD = var(df_plot$varLD_score)

      df_plot$std_varLD_score = (df_plot$varLD_score - mean_varLD)/sqrt(var_varLD)

      #dev.off()

      perc_95 = quantile(df_plot$std_varLD_score, probs = 0.95)
      df_plot$ind = ifelse(df_plot$std_varLD_score > perc_95, 1, 0) #an indicator for if the std varld score is greater than the 95% percentile

      print('Starting to write the varLD data to a table')

      write.table(df_plot, paste('/proj/yunligrp/users/dinelka/scripts/LD_calc_v2/3.0_VarLD/', df_name_1, '_', df_name_2, '/', chrom, '/', chrom, '_ev', n_ev, '.txt', sep = ""), row.names = FALSE)

      plot(df_plot$window, df_plot$varLD_score,type = "b", ylab = "varLD Score", xlab = "Window (500kb)", main = paste("varLD Score Distribution: ", df_name_1, "_", df_name_2, sep=""),
              col = "blue", )

    }
  }
}

# end time
end_time = proc.time() - start_time
print(end_time)
