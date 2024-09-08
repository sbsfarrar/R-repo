library(tidyNano)
library(tidyverse)
library(GrpString)
library(DescTools)

saveDir = paste0(getwd(),"/R/NTA Scripts")
path = file.choose()
#write script to change names to template format
#SampleID_Treatment_Dilution_Rep

data <- nanoimport(path, auto_name = TRUE)

##get most common string in samples

names <- colnames(data[2:ncol(data)])
dil_split <- strsplit(names, "_")
dil_names <- lapply(dil_split, function(x) tail(x, n=1))
dil_names <- unlist(dil_names)
#dil_names <- as.character(rep(1,15))
a <- CommonPatt(names, low=50)
b <- a[which(a$Freq_str == max(a$Freq_str)),]
c <- b[which(b$Length == max(b$Length)),]
sample_base = str_trim(str_remove(c$Pattern, "Rep"))
sample_base_vec = rep(sample_base,length(names))

run_append = length(names)/5 #number of runs
rep_append = length(names)/run_append #number of reps
run_append_str = as.character(seq(1,run_append))
rep_append_str = c("01","02","03","04","05")

run_append_str_vec = rep(run_append_str,each = rep_append)
rep_append_str_vec = rep(rep_append_str,run_append)

d = paste(sample_base_vec, dil_names,sep="_")
e = paste(d,run_append_str_vec,sep="_")
f = paste(e,rep_append_str_vec,sep="_")

colnames(data)[2:ncol(data)] <- f

####a <- "WWDUISBURG-HAMBORNS"
####b <- "QQQQQQDUISBURG (-31.7.29)S"

####A <- strsplit(a, "")[[1]]
####B <- strsplit(b, "")[[1]]

####L <- matrix(0, length(A), length(B))
####ones <- which(outer(A, B, "=="), arr.ind = TRUE)
####ones <- ones[order(ones[, 1]), ]
####for(i in 1:nrow(ones)) {
####  v <- ones[i, , drop = FALSE]
####  L[v] <- ifelse(any(v == 1), 1, L[v - 1] + 1)
####}
####paste0(A[(-max(L) + 1):0 + which(L == max(L), arr.ind = TRUE)[1]], collapse = "")

tidy_data <- data %>%
  nanotidy(sep_var = c("Sample", "Dilution", "Injection", "Tech_rep"))

groups <- tidy_data %>%
  nanocount(Sample,Injection,Tech_rep,name="Total",param_var = True_count)

test_under200 <- tidy_data %>%
  group_by(Sample,Injection,Tech_rep) %>%
  group_map(~ AUC(.x$particle_size[0:201],.x$True_count[0:201]))
  
sub200_count <- as.data.frame(unlist(test_under200))
colnames(sub200_count) <- "under_200nm"
groups = cbind(groups, sub200_count)

test_over200 <- tidy_data %>%
  group_by(Sample,Injection,Tech_rep) %>%
  group_map(~ AUC(.x$particle_size[202:1000],.x$True_count[202:1000]))

super200_count <- as.data.frame(unlist(test_over200))
colnames(super200_count) <- "over_200nm"
groups = cbind(groups, super200_count)

total_df <- groups %>%
  nanolyze(Sample,name="Total",param_var=Total_count, sci_not = TRUE)
sub200_df <- groups %>%
  nanolyze(Sample,name="under_200nm",param_var=under_200nm, sci_not = TRUE)
super200_df <- groups %>%
  nanolyze(Sample,name="over_200nm",param_var=over_200nm, sci_not = TRUE)

final_df <- cbind(total_df, sub200_df[,3:6], super200_df[,3:6])

nanosave(tidy_data, name = paste0(saveDir, "/data"))