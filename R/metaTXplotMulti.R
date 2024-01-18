#' plot multiple groups of features.
#' @export metaTXplotMulti
metaTXplotMulti <-
function(    remapresults_list,
               num_bin              = 10,
               includeNeighborDNA   = TRUE,
               relativeProportion   = c(1, 1, 1, 1),
               lambda = 2,
               adjust = 0.15,
               title  = ''
  ){
    Iter <- length(remapresults_list) 
    if(includeNeighborDNA){
      num_bin_fiveUTR        <- num_bin
      num_bin_cds            <- num_bin
      num_bin_threeUTR       <- num_bin
      num_bin_promoter       <- num_bin
      num_bin_tail           <- num_bin
      num_bin_sum            <- num_bin * 5
      
      coord_interval         <- relativeProportion / min(relativeProportion)
      
      coord_promoter_interval<- coord_interval[4] * num_bin
      coord_promoter_Start   <- coord_interval[4] / 2
      coord_promoter_end     <- coord_promoter_interval - coord_interval[4] / 2
      coord_promoter         <- seq(coord_promoter_Start, coord_promoter_end,  coord_interval[4])
      
      coord_5utr_interval    <- coord_interval[1] * num_bin
      coord_5utr_Start       <- coord_promoter_interval + coord_interval[1] / 2
      coord_5utr_end         <- coord_promoter_interval + coord_5utr_interval - coord_interval[1] / 2
      coord_5utr             <- seq(coord_5utr_Start, coord_5utr_end, coord_interval[1])
      
      coord_cds_interval     <- coord_interval[2] * num_bin
      coord_cds_Start        <- coord_promoter_interval + coord_5utr_interval + coord_interval[2] / 2
      coord_cds_end          <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval - coord_interval[2] / 2
      coord_cds              <- seq(coord_cds_Start, coord_cds_end, coord_interval[2])
      
      coord_3utr_interval    <- coord_interval[3] * num_bin
      coord_3utr_Start       <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_interval[3] / 2
      coord_3utr_end         <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval - coord_interval[3] / 2
      coord_3utr             <- seq(coord_3utr_Start, coord_3utr_end, coord_interval[3])
      
      coord_tail_interval    <- coord_interval[4] * num_bin
      coord_tail_Start       <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval + coord_interval[4] / 2
      coord_tail_end         <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval + coord_tail_interval - coord_interval[4] / 2
      coord_tail             <- seq(coord_tail_Start, coord_tail_end, coord_interval[4])
      
      xintercept_1           <- c(coord_promoter_interval, 
                                  coord_promoter_interval + coord_5utr_interval, 
                                  coord_promoter_interval + coord_5utr_interval + coord_cds_interval,
                                  coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval)  
      coord <- c(coord_promoter, coord_5utr, coord_cds, coord_3utr, coord_tail)
      
      width_coord            <- coord_promoter_interval + coord_5utr_interval + coord_cds_interval + coord_3utr_interval + coord_tail_interval
      coord                  <- coord / width_coord
      xintercept_1           <- xintercept_1 / width_coord
      coord_promoter_interval<- coord_promoter_interval / width_coord
      coord_5utr_interval    <- coord_5utr_interval / width_coord
      coord_cds_interval     <- coord_cds_interval / width_coord
      coord_3utr_interval    <- coord_3utr_interval / width_coord
      coord_tail_interval    <- coord_tail_interval / width_coord
      
    }else{
      
      num_bin_fiveUTR        <- num_bin
      num_bin_cds            <- num_bin
      num_bin_threeUTR       <- num_bin
      num_bin_sum            <- num_bin * 3
      
      coord_interval         <- relativeProportion / min(relativeProportion)
      
      coord_5utr_interval    <- coord_interval[1] * num_bin
      coord_5utr_Start       <- coord_interval[1] / 2
      coord_5utr_end         <- coord_5utr_interval - coord_interval[1] / 2
      coord_5utr             <- seq(coord_5utr_Start, coord_5utr_end, coord_interval[1])
      
      coord_cds_interval     <- coord_interval[2] * num_bin
      coord_cds_Start        <- coord_5utr_interval + coord_interval[2] / 2
      coord_cds_end          <- coord_5utr_interval + coord_cds_interval - coord_interval[2] / 2
      coord_cds              <- seq(coord_cds_Start, coord_cds_end, coord_interval[2])
      
      coord_3utr_interval    <- coord_interval[3] * num_bin
      coord_3utr_Start       <- coord_5utr_interval + coord_cds_interval + coord_interval[3] / 2
      coord_3utr_end         <- coord_5utr_interval + coord_cds_interval + coord_3utr_interval - coord_interval[3] / 2
      coord_3utr             <- seq(coord_3utr_Start, coord_3utr_end, coord_interval[3])
      
      xintercept_1           <- c(coord_5utr_interval, 
                                  coord_5utr_interval + coord_cds_interval)  
      
      coord <- c(coord_5utr, coord_cds, coord_3utr) 
      
      width_coord            <- coord_5utr_interval + coord_cds_interval + coord_3utr_interval 
      coord                  <- coord / width_coord
      xintercept_1           <- xintercept_1 / width_coord
      coord_5utr_interval    <- coord_5utr_interval / width_coord
      coord_cds_interval     <- coord_cds_interval / width_coord
      coord_3utr_interval    <- coord_3utr_interval / width_coord
    }
    
    densities_list = list()
    for(j in 1:Iter){
      align_mtr          <- remapresults_list[[j]][[1]]
      width_mtr          <- remapresults_list[[j]][[2]]
      trans_info         <- remapresults_list[[j]][[3]]
      num_bin_sum        <- ncol(align_mtr)
      
      if(includeNeighborDNA){
        weight_start     <- num_bin_promoter + 1
        weight_end       <- num_bin_promoter + num_bin_fiveUTR + num_bin_cds + num_bin_threeUTR
        weight_mtr       <- replicate(num_bin_sum, rowSums(width_mtr[, weight_start:weight_end]))    
        alpha            <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[1]]
        densities <- alpha / colSums(width_mtr) / sum(alpha / colSums(width_mtr))
        
        }else{
        weight_mtr       <- replicate(num_bin_sum, rowSums(width_mtr))
        alpha            <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[1]]
        densities <- alpha / colSums(width_mtr) / sum(alpha / colSums(width_mtr))
      } 
      densities_list[[j]] = densities
    }

    data_plot_list <- list()
    if(includeNeighborDNA){

      for(k in 1:Iter){
        value1 <- densities_list[[k]]
        data_plot            <- data.frame(coord  = coord
                                           , value =  value1
                                           , type  = paste0('group_', k))
        colnames(data_plot)  <- c('coord', 'value', 'type')
        row.names(data_plot) <- 1:nrow(data_plot)
        
        
        coord_interval_seq   <- c(replicate(num_bin,coord_promoter_interval / num_bin),
                                  replicate(num_bin,coord_5utr_interval / num_bin),
                                  replicate(num_bin,coord_cds_interval / num_bin),
                                  replicate(num_bin,coord_3utr_interval / num_bin),
                                  replicate(num_bin,coord_tail_interval / num_bin))
        
        max_coord            <- max(data_plot[, 2])
        max_coord_relative   <- max_coord / sum(data_plot[,2] * coord_interval_seq) 
        max_coord_relative   <- max_coord_relative / adjust  
        
        data_plot_1            <- data_plot
        plot_values            <- data_plot[, 2]
        plot_values            <- c(plot_values[1:num_bin] * coord_promoter_interval,
                                    plot_values[1:num_bin + num_bin] * coord_5utr_interval,
                                    plot_values[1:num_bin + num_bin*2] * coord_cds_interval,
                                    plot_values[1:num_bin + num_bin*3] * coord_3utr_interval,
                                    plot_values[1:num_bin + num_bin*4] * coord_tail_interval)
        plot_values            <-  plot_values / sum(plot_values)
        data_plot_1[, 2]       <- plot_values
        data_plot_list[[k]]    <- data_plot_1
      }
      
      data_plot = c()
      for(i in 1:Iter){
        data_plot = rbind(data_plot, data_plot_list[[i]])
      }
      
        p1 <- 
          ggplot(data_plot, aes(x=coord, group=type, weight=value)) +
          ggtitle(title) + 
          theme(panel.background =element_blank(),
                panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
                axis.text.x = element_blank(), axis.ticks = element_blank(),
                line = element_line(colour = "white", size = 0.5, linetype = 9, lineend = "butt"),
                legend.position = 'bottom', 
                axis.text=element_text(size = 8),
                legend.text=element_text(size = 8),
                axis.title.y =element_text(size = 8 ,hjust=0.5),
                title = element_text(size = 8, face='bold')) + 
          geom_density(adjust = adjust, aes(fill=factor(type), colour = factor(type)),alpha = 0.15) +
          xlab("") +
          ylab("Density
               ") +
          annotate("text", x =  xintercept_1[1] / 2, y = -0.007 * max_coord_relative, label = "Promoter", size = 3) +
          annotate("text", x =  xintercept_1[1] + coord_5utr_interval / 2, y = -0.007 * max_coord_relative, label = "5'UTR", size = 3) +
          annotate("text", x =  xintercept_1[2] + coord_cds_interval / 2, y = -0.007 * max_coord_relative, label = "CDS", size = 3) +
          annotate("text", x =  xintercept_1[3] +  coord_3utr_interval / 2, y = -0.007 * max_coord_relative, label = "3'UTR", size = 3) +
          annotate("text", x =  xintercept_1[4] +  coord_tail_interval / 2, y = -0.007* max_coord_relative, label = "Tail", size = 3) + 
          geom_vline(xintercept= xintercept_1, linetype = 9, size = 0.2,colour = 'black') +
          annotate("rect", xmin = xintercept_1[1], xmax = xintercept_1[2], ymin = -0.0032 * max_coord_relative, ymax = -0.0025 * max_coord_relative, alpha = .8, colour = "black")+
          annotate("rect", xmin = xintercept_1[2], xmax = xintercept_1[3], ymin = -0.0042 * max_coord_relative, ymax = -0.0017 * max_coord_relative, alpha = .3, colour = "black")+
          annotate("rect", xmin = xintercept_1[3], xmax = xintercept_1[4], ymin = -0.0032 * max_coord_relative, ymax = -0.0025 * max_coord_relative, alpha = .8, colour = "black")
      
      
    }else{
      
        
        for(k in 1:Iter){
        value1 <- densities_list[[k]]
        data_plot            <- data.frame(coord  = coord
                                           , value = value1
                                           , type  = paste0('group_', k))
        colnames(data_plot)  <- c('coord', 'value', 'type')
        row.names(data_plot) <- 1:nrow(data_plot)
        
        
        coord_interval_seq   <- c(replicate(num_bin,coord_5utr_interval / num_bin),
                                  replicate(num_bin,coord_cds_interval / num_bin),
                                  replicate(num_bin,coord_3utr_interval / num_bin))
        
        max_coord            <- max(data_plot[, 2])
        max_coord_relative   <- max_coord / sum(data_plot[,2] * coord_interval_seq) 
        max_coord_relative   <- max_coord_relative / adjust  
        
        data_plot_1            <- data_plot
        plot_values            <- data_plot[, 2]
        plot_values            <- c(plot_values[1:num_bin] * coord_5utr_interval,
                                    plot_values[1:num_bin + num_bin*1] * coord_cds_interval,
                                    plot_values[1:num_bin + num_bin*2] * coord_3utr_interval)
        plot_values            <-  plot_values / sum(plot_values)
        data_plot_1[, 2]       <- plot_values
        data_plot_list[[k]]    <- data_plot_1
        }
      
      data_plot = c()
      for(i in 1:Iter){
        data_plot = rbind(data_plot, data_plot_list[[i]])
      }
         p1 <- 
          ggplot(data_plot, aes(x=coord, group=type, weight=value)) +
          ggtitle(title) +
          theme(panel.background =element_blank(),
                panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
                axis.text.x = element_blank(), axis.ticks = element_blank(),
                line = element_line(colour = "white", size = 0.2, linetype = 9, lineend = "butt"),
                legend.position = 'bottom', 
                axis.text=element_text(size = 8),
                legend.text=element_text(size = 8),
                axis.title.y =element_text(size = 8 ,hjust=0.5),
                title = element_text(size = 8, face='bold')) + 
          geom_density(adjust = adjust, aes(fill=factor(type), colour = factor(type)),alpha = 0.15) +
          xlab("") +
          ylab("Density
               ") +
          annotate("text", x = xintercept_1[1] / 2, y = -0.007 * max_coord_relative, label = "5'UTR",  size = 3) +
          annotate("text", x = xintercept_1[1] + coord_cds_interval / 2, y = -0.007 * max_coord_relative, label = "CDS",  size = 3) +
          annotate("text", x = xintercept_1[2] +  coord_3utr_interval / 2, y = -0.007 * max_coord_relative, label = "3'UTR",  size = 3) +
          geom_vline(xintercept= xintercept_1, linetype = 9, size = 0.2,colour = 'black') +
          annotate("rect", xmin = num_bin_sum * 0,  xmax = xintercept_1[1], ymin = -0.0032 * max_coord_relative, ymax = -0.0025 * max_coord_relative, alpha = .8, colour = "black")+
          annotate("rect", xmin = xintercept_1[1], xmax = xintercept_1[2], ymin = -0.0042 * max_coord_relative, ymax = -0.0017 * max_coord_relative, alpha = .3, colour = "black")+
          annotate("rect", xmin =xintercept_1[2], xmax = coord_5utr_interval + coord_cds_interval + coord_3utr_interval, ymin = -0.0032 * max_coord_relative, ymax = -0.0025 * max_coord_relative, alpha = .8, colour = "black")

    } 
    return(p1)
  }
