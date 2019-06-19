##################
## Directories & files management functions. 
#################

#############################
##  Compare_dir_content
##  Function to detect the instances that differ between two given directories. Using one as reference. 
##  Created: 19-02-12
##  Args: ref_dir_path, dir2_path, strings, rec (recursive) logical, default F.
############################

ref_dir_path = 'haplotypes2'
dir2_path = "haplotypes2 (Selective Sync Conflict)"

compare_dir_content = function(ref_dir_path, dir2_path, rec = F){
  if(!dir.exists(ref_dir_path) | !dir.exists(dir2_path)) stop('One or both directory paths do not exist')
  ls_dir1 = list.files(ref_dir_path)
  ls_dir2 = list.files(dir2_path)
  if(!rec){
    return( setdiff(ls_dir2, ls_dir1) )
  }else{
    ##  Easies way - find each of the files in dir2 and then check if each one of them is in the other directory. 
    ls_dir2 = list.files(dir2_path, recursive = T)
    ls_status = sapply(file.path(ref_dir_path, ls_dir2), file.exists)
    ls_dir2 = ls_dir2[!ls_status]
    return(ls_dir2)
  }
  ##  Return the files that are not present in the 
}

compare_dir_content(ref_dir_path, dir2_path, rec = T)

#############################
##  eliminate_rec_dir_content
##  Function eliminate directories in the same org level with a given name and all their content. (fixed dir file in a complex directory tree)
##  Created: 19-02-12
##  Args: path, dir_name, strings, dirnames can be a vector. 
############################
#path = './haplotypes2'; dir_name = c( "allInputReadsForEachCluster", "clusters", "internalSnpInfo")
path = "/Volumes/All_staff/Investigators/Permar/DATA/Former Lab Students & Personnel/Diana_VC/haplotypes2"
path = '/Users/dveracr/Dropbox\ \(Duke\ Bio_Ea\)/rhCMV/haplotypes2\ \(Selective\ Sync\ Conflict\)'
eliminate_rec_dir = function(path, dir_name){
  if(!dir.exists(path)) stop('Directory does not exist')
  if(length(dir_name) > 1) dir_name = paste(dir_name, collapse = "|") ##  If a vector, then collpase to get a single pattern to search 
  flag = F
  dirs = path
  while(!flag){
    dirs = lapply(dirs, list.dirs, recursive = F, full.name = T)
    dirs = unlist(dirs)
    if(length(dirs) == 0){
      return('No directory name found in the directories structure')
    }
    if( sum(grepl(dir_name, dirs) ) > 0 ) flag = T
  }
  ##  Once it stops, we select only those files with the pattern and erase them and their content. 
  dirs = grep(dir_name, dirs, val = T)
  sapply(paste('rm -r', dirs),  system)
  return(T)
}

#############################
##  dir_tree
##  Function creates a recursive series of list (tree) to  make it easier to copy, remove and edit directories. 
##  Created: 19-02-12
##  Args: ls_files, vector of strings. 
############################

##  Pass all the files from one directory to another. 

dir_tree = function(ls_files, d_to, d_from){
  library(magrittr)
  ##  Create a list of lists with the files paths. 
  dirs_ls = strsplit(ls_files, '/')
  n_levels = max( sapply(dirs_ls, length) ) ##  Get the max number of levels. 
  dir_m = matrix(NA, nrow = length(ls_files), ncol = n_levels)
  for(i in 1:length(ls_files)){
    tmp = dirs_ls[[i]]
    dir_m[i, 1:length(tmp)] = tmp
  }
  ##  Detect files from directories. _. is na the  next column entry. 
  is_file = cbind( apply(dir_m, 2, is.na)[,-1], T)
  
  ####  Now, check that the directories exist from length to right. 
  for(i in 2:cols(dir_m)){
    ##  Get all the unique characters up to the i column. 
    if(i > 1){
      tmp = apply(dir_m[,1:i], 1, paste, collapse= '/') 
    }else{
      tmp = dir_m[,1]
    }
    ##  Likely, things that appear at least 2 should be directories, the other ones have to be sorted. 
    ## Create directories if not existent.
    new_dirs = NULL
    new_dirs = tmp[!is_file[,i]] %>% unique
    if(length(new_dirs) > 0) sapply(file.path(d_to, new_dirs), function(j) if(!dir.exists(j)) dir.create(j)  )
    ##  Pass files 
    if(length(tmp[is_file[,i]]) > 0 ) sapply(tmp[is_file[,i]], function(j) file.copy(file.path(d_from, j), file.path(d_to, j), overwrite = T) )
  }
  
}
