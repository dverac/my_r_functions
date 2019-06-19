##  Detecting chimeric sequences. 

## Function to find potential quimeric haplotypes sequences. 
#### h_seq -> DNAStringSet with sequences to investigate -> SAME LENGTH.
#### refs -> Names of the sequences in h_seq to use as references. 
## Returns potential references pairs that can explain haplotype as chimeric instead of de novo mutations. 
find_chimeras = function(h_seq, refs){
  d = as.matrix(ape::dist.dna(as.DNAbin(h_seq), 'raw'))*width(h_seq)[1]
  #print(d[refs, refs])
  pairs = data.frame()
  novo = dif_ref = dif_alt = list()
  ##  Found possible recombinants per each of the sequences used for the query. 
  for(h in setdiff(rownames(d),refs)){
    ix = order(d[h, refs])
    ## Distance between haplotype h to the other parental strains = distance of parentakl strain of h  to other parental seqs + distance of h to its parental.
    #d[h, refs[ix[-1]]] == (d[refs[ix[1]], refs[ix[-1]]] + d[h, refs[ix[1]]])   c = a + b  
    ## If equal: No recombination.  # If different + Recombination or new mutation.
    a = round( d[h, refs[ix[1]]] )  # distance from hap h to parental reference.
    b = round( d[refs[ix[1]], refs[ix[-1]]] )  # distance from parental ref to other non parental references (alt).
    c = round( d[h, refs[ix[-1]]] ) # distance from hap h to other non parental references.
    ix2 = ix[-1][c != (a+b)]  ## Number of other references non parental where the distance doesn't follow c = a + b
    if (length(ix2)>0){
      ##Either solving the system of equations or obtaining the variable sites, the results will be the same, so let's use sets XD
      # z = (b+c-a)/2; x = b-z; y = c-z
      ## Applying sets.
      for (alt in refs[ix2]){
        temp = lapply(h_seq[c(h, refs[ix[1]], alt)], function(i) consensusMatrix(DNAMultipleAlignment(i))[1:4,])
        #which(apply(temp[[1]]+temp[[2]]+temp[[3]], 2, function(i) sum(i>0)) > 1) # total var sites
        a = which(apply(temp[[1]]+temp[[2]], 2, function(i) sum(i>0)) > 1) # var sites between h and ref
        c = which(apply(temp[[1]]+temp[[3]], 2, function(i) sum(i>0)) > 1) # var sites between h and alt
        a2 = sapply(a, function(i) paste(rownames(temp[[1]])[temp[[1]][,i]==1], i, rownames(temp[[2]])[temp[[2]][,i]==1], sep=':') )  ## Adding the nt
        c2 = sapply(c, function(i) paste(rownames(temp[[1]])[temp[[1]][,i]==1], i, rownames(temp[[3]])[temp[[3]][,i]==1], sep=':') )  ## Adding the nt
        y2 = intersect(a2,c2) # de novo mutations
        y = as.numeric(gsub('A|T|C|G|:','',y2))
        x = sort(setdiff(a, y))
        z = sort(setdiff(c, y))
        or = order(c(min(x), max(x),min(z), max(z)-1))
        #or = paste(or2, collapse='')  ## minus 1 in max(z) to account for the same position apparently mutated.
        if(paste(or, collapse='') %in% c('1234','3412')){
          novo[[paste(h, refs[ix[1]], alt, sep=':')]] = y2
          dif_ref[[paste(h, refs[ix[1]], alt, sep=':')]] = setdiff(a2, y2)
          dif_alt[[paste(h, refs[ix[1]], alt, sep=':')]] = setdiff(c2, y2)
          bk = mean( c(min(x), max(x),min(z), max(z))[or[2:3]] )
          pairs = rbind(pairs, data.frame(hap=h, ref=refs[ix[1]], d_ref=length(a), alt=alt, d_alt=length(c), d_novo= length(y), dif_ref=paste(min(x), max(x), sep='-'), dif_alt=paste(min(z), max(z), sep='-'), bk=bk))
        }
      }
    }
  }
  if(nrow(pairs) == 0){
    return(NA)
  }else{
    pairs$c_ref = dif_ref
    pairs$c_alt = dif_alt
    pairs$novo = novo
    ##  Filtering -> taking the entries where the minimum number of new mutations is needed to explain the haplotype.
    ix = unlist(tapply(1:nrow(pairs), pairs$hap, function(i) i[pairs$d_novo[i] == min(pairs$d_novo[i])]))
    pairs = pairs[ix,]
    pairs = pairs[ sort(tapply(1:nrow(pairs), as.character(pairs$hap), function(i) i[1])), ] ## Take just one possible recombination event per haplotype.
    return(pairs)
  }
}
