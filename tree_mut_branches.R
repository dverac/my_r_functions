#  Find mutations in a tree.
##  General code for broad use in other phylogenies. 
##  Feb 28th 2019

######
## Generate ancestral datastructures ->  ancestral_reconstruction from treetime
##  Phylo tree from file. 
##  Read seq and get the aa sequence too-> aligned. 
##### 

create_ancestral_data = function(tree_file, aln_file, rdata_file = NULL, obj_names = c('tre', 'nt', 'aa')){
  library(tidyverse);library(ape); library(Biostrings)
  ##  Read the data. 
  tre = read.nexus(tree_file)
  nt = readDNAStringSet(aln_file) ##  Make sure the seq are in frame. 
  
  ##  Translate to AA. 
  aa = gsub('-', '', nt) %>% DNAStringSet() %>% translate
  if(length(unique(width(aa))) > 10){
    ##  Make an alignment if the translated sequences are not of the same length. 
    writeXStringSet(aa, 'tmp.fa')
    system("muscle -in tmp.fa -out tmp_aa.fa")
    aa = readAAStringSet("tmp_aa.fa")
    file.remove("tmp.fa")
    file.remove("tmp_aa.fa")
  }
  ##  If rdata-file is null -> return a list 
  if(!is.null(rdata_file)){
    save(tre, nt, aa, file = rdata_file)
    return(paste("Tree and sequence alignment successfully stored in", rdata_file))
  }else{
    assign(obj_names[1], tre, env = globalenv())
    assign(obj_names[2], nt, env = globalenv())
    assign(obj_names[3], aa, env = globalenv())
    return(paste("Tree and sequence alignments sucessfully loaded in global environment: ", obj_names))
  }
}


######
##  Find mutations between a given sequnece and a reference. StringSets
##  Output dataframe with old, new characters, position and label columns. 
##### 
find_mutations = function(x, ref){
  ##  Check that both are X StringSets. 
  if( ! class(x) %in% c('AAStringSet', 'DNAStringSet') ) stop ("Weird x, not XXStringSet")
  if( ! class(ref) %in% c('AAStringSet', 'DNAStringSet') ) stop ("Weird ref, not XXStringSet")
  if(width(x) != width(ref)) stop()
  ##  Search variable positions and find the prev and new character. 
  cm = consensusMatrix(c(ref,x)) 
  var_pos = apply(cm, 2, function(i) sum(i > 0) > 1) %>% which
  ## Get characters. 
  old_char = sapply(var_pos, function(i) substr(ref, i, i))
  new_char = sapply(var_pos, function(i) substr(x, i, i))
  if(length(var_pos) == 0) old_char = new_char = vector()
  mut = data.frame(pos= var_pos, old = old_char, new = new_char, 
                   label = paste(old_char, var_pos, new_char, sep = ''))
  return( mut )
}

######
## Map mutations to tree branches. 
##  Input: tre - phylo; nt and aa DNA or AAStringsets. 
##  Output: tree tibble with branch.mut_aa / nt with mutation labels. 
##### 
mutations_per_node = function(tre, nt, aa = FALSE){
  library(tidyverse);library(tidytree); library(ape)
  #####
  ##Checking input
  #####  
  ##  Check the tree is a phylo or tree tible str.
  if("tbl_tree" %in% class(tre)){
    tree_dt = tre
  }else if(class(tre) == 'phylo'){
    tree_dt = as_tibble(tre)
  }else{
    stop("Tree is not in a phylo or tbl tree data structure")
  }
  ## Check the seq is a DNAStringSet and all samples have same length.
  if(class(nt) != 'DNAStringSet') stop('Nt is not a DNAStringSet')
  if(length(unique(width(nt))) > 1) stop('Sequences of different length are present')
  ##  Check the names match with the seq.
  if( length(setdiff(tree_dt$label, names(nt)) ) > 0 | nrow(tree_dt) < length(nt) ) stop("Missing entries")
  nt = nt[tree_dt$label] ##  Keep same order in tree and in sequences. 
  ## Check if aminoacid level stuff should be checked too. 
  if(class(aa) == 'logical'){
    if(!aa){
      lvs = list(nt = nt)
    }else{
      lvs = list(nt = nt, aa = translate(nt))
    }
  }else if(class(aa) == 'AAStringSet'){
    lvs = list(nt = nt, aa = aa)
  }else{
    stop('Incorrect aa parameter')
  }
  #####
  ## Main code. 
  #####  
  ##  Put sequences in same order as the tips/nodes in the tibble. 
  for(ch in names(lvs))   lvs[[ch]] = lvs[[ch]][tree_dt$label]
  ##  Add the names of the columns we need. 
  tree_dt[,paste('branch.mut', names(lvs), sep = '_')] = ''
  
  for(nn in nrow(tree_dt):1 ){
    ##  Method: Just compare the sample and their parents  unless I am already in the tree root.
    ##  Keep stuff in index value ->  node is index for the dataframe
    nn_dad = parent(tree_dt, nn)$node
    if(length(nn_dad) ==0 ) next
    for(ch in names(lvs)){
      ##  Get sequences for node and their parent. 
      x = lvs[[ch]][nn]
      ref = lvs[[ch]][nn_dad]
      if(x == ref) tmp = data.frame() else tmp = find_mutations(x, ref)
      if(nrow(tmp) > 0){
        tree_dt[nn, paste('branch.mut', ch, sep = '_')] = paste(tmp$label, collapse = "-")
      }
    }
  }
  return(tree_dt)
}




