get_trans_function <-
  function(trans, methyl){
    sample_methyl_pos      <- start(ranges(methyl))
    sample_methyl_seqnames <- data.frame(seqnames(methyl))[[1]]
    trans_metadata         <- mcols(trans)
    
    index_trans            <- 1:length(trans)
    index_methyl           <- data.frame(trans)[, 'xHits']
    seqnames               <- sample_methyl_seqnames[data.frame(trans)[, 'xHits']]
    methyl_pos             <- sample_methyl_pos[data.frame(trans)[, 'xHits']]
    strand                 <- data.frame(trans)[, 'strand']
    trans_ID               <- data.frame(trans)[, 'transcriptsHits']
    
    trans_info             <- data.frame(index_trans  = index_trans,
                                         index_methyl = index_methyl,
                                         seqnames     = seqnames, 
                                         methyl_pos   = methyl_pos, 
                                         strand       = strand,
                                         trans_ID     = trans_ID)
    return(trans_info)
  }

