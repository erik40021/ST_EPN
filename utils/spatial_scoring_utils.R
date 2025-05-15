# all utility functions needed to calculate spatial scores

# ------ PART 1: spatial complexity and coherence functions ------

# calculates spatial gradients for a given sample and window size (-> main, most time-intensive scoring function, recommended to run on HPC)
calculate_spatial_complexity = function(s, w, rand_num, inner_cores) {
  message("calculating 'spatial gradients by window' for sample ", s, " using params w = ", w, ", rand_num = ", rand_num, 
  ", inner_cores = ", inner_cores, " (start time: ", format(Sys.time(), "%H:%M"), ")")
  
  all_spots_progs = all_spots_programs_comp[[s]]
  mp_names = colnames(all_spots_progs)[-ncol(all_spots_progs)]
  all_win_neighbors_table = neighbors_table_funcV2(all_spots_positions[[s]], all_spots_progs)
  neighbors_table = neighbors_table_func_comp(all_spots_positions[[s]], all_spots_progs, mp_names)

  # calculate spatial gradients per spot by measuring the average ratio of spot-to-neighbours in the window around a spot and comparing to min and max expectance estimates
  program_gradients_list = mclapply(rownames(all_spots_progs), mc.cores = inner_cores, function(spot) {
    win_spots = c(spot)
    sapply(1:w, function(j) {
      win_spots <<- unique(c(win_spots,unique(na.omit(as.character(all_win_neighbors_table[win_spots,])))))
      win_spots <<- win_spots[win_spots != "NaN"]
    })
    if (length(win_spots) == 1) return(NA) # no neighbours exist in window
    
    win_rand_neighbors_table_comp = lapply(1:rand_num, function(k) { # shuffling spot positions in reference window to calculate expected min coherence
      win_spots_positions = all_spots_positions[[s]][all_spots_positions[[s]]$V1 %in% win_spots,]
      new_pos = sample(win_spots_positions$V1, length(win_spots_positions$V1), replace = FALSE)
      pos_table = win_spots_positions
      row.names(pos_table) = new_pos 
      pos_table$V1 = new_pos
      pos_table$V2 = win_spots_positions[new_pos, "V2"]
      win_neighbors_table = win_prox_neighbors_table_func_comp(pos_table, all_spots_progs[row.names(pos_table), ])
      return(win_neighbors_table)
    })
    
    all_win_program_abund = all_spots_progs[win_spots, ] # all spots within the window and their programs
    
    programs_gradients_in_spot = sapply(1:length(mp_names), function(prog_ind) { # calculate spatial gradients separately per program
      win_spots_program_scores = all_win_program_abund[, prog_ind]; names(win_spots_program_scores) = rownames(all_win_program_abund)
      win_spots_indices = str_detect(rownames(neighbors_table), paste(paste0(win_spots, ".", prog_ind, "\\b"), collapse = "|")) # fastest to find spots indices using collapsed regex pattern
      program_neighbors_table = neighbors_table[win_spots_indices, ] # all spots within the window and their neighbours
      rownames(program_neighbors_table) = substr(rownames(program_neighbors_table), start = 0, stop = 18)
      program_neighbors_table = program_neighbors_table[win_spots, ] # rename and reorder to match win_spots_program_scores
      
      observed = find_program_coherence(program_neighbors_table, win_spots_program_scores)
      expected_min = get_expected_min(win_rand_neighbors_table_comp, win_spots_program_scores, prog_ind)
      if (observed < expected_min) observed = expected_min # edge case, score will equal 0 <- passiert eigentlich nie, sehr unwahrscheinlich
      gradients = (observed - expected_min) / (1 - expected_min)
      return(gradients)
    })
    return(programs_gradients_in_spot)
  })
  message("all ", nrow(all_spots_progs), " spots done in ", s)
  program_gradients = as.data.frame(t(do.call(cbind, program_gradients_list))); rownames(program_gradients) = rownames(all_spots_progs); colnames(program_gradients) = mp_names
  return(program_gradients)
}


# calculates spatial coherence in a sample per spot and program. coherence in a spot for a program := average ratio of a spot's and its neighbours' program scores (fast enough to be run locally. no need for HPC)
calculate_spatial_coherence = function(s, mp_names, smoothing_win_size = 4, use_norm_scores = T) {
  message("calculating 'spatial coherence by program' for sample ", s)
  
  all_spot_programs = if (use_norm_scores) all_spots_programs_comp_norm[[s]] else all_spots_programs_comp[[s]]  
  all_neighbours = neighbors_table_funcV2(all_spots_positions[[s]], all_spot_programs); all_neighbours[all_neighbours == "NaN"] = NA
  all_neighbours_scores = neighbors_table_func_comp(all_spots_positions[[s]], all_spot_programs)
  
  program_coherences = all_spot_programs[, -(length(mp_names)+1)]; program_coherences[, ] = NA # empty df to store scores of each program per spot
  for (spot in rownames(all_spot_programs)) {
    spot_program_coherence = c()
    for (program in mp_names) {
      spot_program_score = all_spot_programs[spot, program]
      neighbours_scores = na.omit(as.double(all_neighbours_scores[paste(paste0(spot, ".", which(program == mp_names))), ]))
      if (length(neighbours_scores) == 0) { spot_program_coherence[program] = NA; next }
      neighbours_ratios = sapply(neighbours_scores, function(neigh_score) { # get every spot-to-neighbour coherence individually
        ratio = spot_program_score / neigh_score
        scaled_ratio = 1 - abs((ratio - 1) / (ratio + 1))
        return(scaled_ratio)
      })
      coherence = mean(neighbours_ratios, na.rm = T) * spot_program_score
      spot_program_coherence[program] = coherence
    }
    program_coherences[spot, ] = spot_program_coherence
  }
  # smoothen program coherences
  program_coherences_smo = program_coherences
  all_spots = rownames(all_neighbours)
  for (spot in all_spots) {
    win_spots = c(spot)
    for (j in 1:smoothing_win_size) win_spots = unique(c(win_spots, unique(na.omit(as.character(all_neighbours[win_spots, ])))))
    for (program in mp_names) program_coherences_smo[spot, program] = mean(program_coherences[win_spots, program], na.rm = T)
  }
  return(program_coherences_smo)
}

find_program_coherence = function(program_neighbors_table, win_spots_program_scores) {
  # calculates coherence of a program in a window. the ratio of a spot to its 6 neighbours is measure for every spot in the window, then averaged
  neighbors_program_coherence = sapply(1:nrow(program_neighbors_table), function(i) {
    spot_program_score = win_spots_program_scores[i]
    neighbours_ratios = sapply(as.double(program_neighbors_table[i, ]), function(neigh_score) { # get every spot-to-neighbour coherence individually
      ratio = spot_program_score / neigh_score
      scaled_ratio = 1 - abs((ratio - 1) / (ratio + 1))
      return(scaled_ratio)
    })
    coherence = mean(neighbours_ratios, na.rm = T)
    return(coherence)
  })
  mean_win_coherence = mean(neighbors_program_coherence, na.rm = T) #[-spot_index], na.rm = T)
  # mean_win_coherence = (4/5) * mean_win_coherence + (1/5) * neighbors_program_coherence[[spot_index]] # (optional) put more weight to the ratio of the window centre spot
  return(mean_win_coherence)
}

get_expected_min = function(all_rand_tables, spots_program_scores, prog_ind) {
  # calculates program-coherence for each randomly shuffled window composition, to obtain a minimum expectance value
  all_zeroval = sapply(all_rand_tables, function(neighbors_rand_table) {
    neighbors_filtered = neighbors_rand_table[grepl(paste0("\\.", prog_ind, "$"), rownames(neighbors_rand_table)), ]
    rownames(neighbors_filtered) = substr(rownames(neighbors_filtered), start = 0, stop = 18)
    rand_spots_program_scores = spots_program_scores[rownames(neighbors_filtered)]
    rand_obs = find_program_coherence(neighbors_filtered, rand_spots_program_scores) #, which(spot == rownames(neighbors_filtered)))
    return(rand_obs)
  })
  zeroval = mean(all_zeroval) # average across all program-coherences of shuffled window configurations
  return(zeroval)
}

# adapted from Greenwald et al. 24
# finds all adjacent/neighbouring spots' metaprograms for each spot
neighbors_table_func_comp = function(spots_positions, spots_programs_comp) {
  neighbors_table = list()
  for (spot in rownames(spots_programs_comp)) {
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    if (spots_col == 0 | spots_row == 0) {
      c1 = rep(NaN, length(mp_names))
    } else {
      n1 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1] # first neighbour
      if (length(n1) == 0 || (spots_positions$V2[spots_positions$V1 == n1] == 0 | !(n1 %in% rownames(spots_programs_comp)))) {
        c1 = rep(NaN, length(mp_names))
      } else {
        c1 = as.character(spots_programs_comp[n1, -ncol(spots_programs_comp)]) # extract all but the last column 'barcodes'
      }
    }
    if (spots_col == 127 | spots_row == 0) {
      c2 = rep(NaN, length(mp_names))
    } else {
      n2 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (length(n2) == 0 || (spots_positions$V2[spots_positions$V1 == n2] == 0 | !(n2 %in% rownames(spots_programs_comp)))) {
        c2 = rep(NaN, length(mp_names))
      } else {
        c2 = as.character(spots_programs_comp[n2, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 0 | spots_col == 1) {
      c3 = rep(NaN, length(mp_names))
    } else {
      n3 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (length(n3) == 0 || (spots_positions$V2[spots_positions$V1 == n3] == 0 | !(n3 %in% rownames(spots_programs_comp)))) {
        c3 = rep(NaN, length(mp_names))
      } else {
        c3 = as.character(spots_programs_comp[n3, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 126 | spots_col == 127) {
      c4 = rep(NaN, length(mp_names))
    } else {
      n4 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (length(n4) == 0 || (spots_positions$V2[spots_positions$V1 == n4] == 0 | !(n4 %in% rownames(spots_programs_comp)))) {
        c4 = rep(NaN, length(mp_names))
      } else {
        c4 = as.character(spots_programs_comp[n4, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 0 | spots_row == 77) {
      c5 = rep(NaN, length(mp_names))
    } else {
      n5 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (length(n5) == 0 || (spots_positions$V2[spots_positions$V1 == n5] == 0 | !(n5 %in% rownames(spots_programs_comp)))) {
        c5 = rep(NaN, length(mp_names))
      } else {
        c5 = as.character(spots_programs_comp[n5, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 127 | spots_row == 77) {
      c6 = rep(NaN, length(mp_names))
    } else {
      n6 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (length(n6) == 0 || (spots_positions$V2[spots_positions$V1 == n6] == 0 | !(n6 %in% rownames(spots_programs_comp)))) {
        c6 = rep(NaN, length(mp_names))
      } else {
        c6 = as.character(spots_programs_comp[n6, -ncol(spots_programs_comp)])
      }
    }
    neighbors_table[[spot]] = data.frame(c1, c2, c3, c4, c5, c6)
  }
  neighbors_table <- do.call(rbind, neighbors_table)
  return(neighbors_table)
}


# adapted from Greenwald et al. 24
# finds all adjacent/neighbouring spots for each spot
neighbors_table_funcV2 = function(spots_positions, spots_programs) {
  neighbors_table = sapply(spots_programs$barcodes, function(spot) {
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    if (spots_col == 0 | spots_row == 0) {
      c1 = NaN
    } else {
      n1 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (length(n1) == 0 || (spots_positions$V2[spots_positions$V1 == n1] == 0 | !(n1 %in% spots_programs$barcodes))) {
        c1 = NaN
      } else {
        c1 = as.character(n1)
      }
    }
    if (spots_col == 127 | spots_row == 0) {
      c2 = NaN
    } else {
      n2 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (length(n2) == 0 || (spots_positions$V2[spots_positions$V1 == n2] == 0 | !(n2 %in% spots_programs$barcodes))) {
        c2 = NaN
      } else {
        c2 = as.character(n2)
      }
    }
    if (spots_col == 0 | spots_col == 1) {
      c3 = NaN
    } else {
      n3 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (length(n3) == 0 || (spots_positions$V2[spots_positions$V1 == n3] == 0 | !(n3 %in% spots_programs$barcodes))) {
        c3 = NaN
      } else {
        c3 = as.character(n3)
      }
    }
    if (spots_col == 126 | spots_col == 127) {
      c4 = NaN
    } else {
      n4 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (length(n4) == 0 || (spots_positions$V2[spots_positions$V1 == n4] == 0 | !(n4 %in% spots_programs$barcodes))) {
        c4 = NaN
      } else {
        c4 = as.character(n4)
      }
    }
    if (spots_col == 0 | spots_row == 77) {
      c5 = NaN
    } else {
      n5 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (length(n5) == 0 || (spots_positions$V2[spots_positions$V1 == n5] == 0 | !(n5 %in% spots_programs$barcodes))) {
        c5 = NaN
      } else {
        c5 = as.character(n5)
      }
    }
    if (spots_col == 127 | spots_row == 77) {
      c6 = NaN
    } else {
      n6 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (length(n6) == 0 || (spots_positions$V2[spots_positions$V1 == n6] == 0 | !(n6 %in% spots_programs$barcodes))) {
        c6 = NaN
      } else {
        c6 = as.character(n6)
      }
    }
    return(c(c1,c2,c3,c4,c5,c6))
  })
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_programs$barcodes
  return(neighbors_table)
}

# adapted from Greenwald et al. 24
win_prox_neighbors_table_func_comp = function(spots_positions, spots_programs_comp) {
  neighbors_table = list()
  for (spot in rownames(spots_programs_comp)) {
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    if (spots_col == 0 | spots_row == 0) {
      n1 = NA
    } else {
      n1_temp = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (length(n1_temp) == 0) {
        n1 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n1_temp] == 0 | !(n1_temp %in% rownames(spots_programs_comp))){
        n1 = NA
      } else {
        n1 = as.character(spots_programs_comp[n1_temp, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 127 | spots_row == 0) {
      n2 = NA
    } else {
      n2_temp = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (length(n2_temp) == 0) {
        n2 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n2_temp] == 0 | !(n2_temp %in% rownames(spots_programs_comp))){
        n2 = NA
      } else {
        n2 = as.character(spots_programs_comp[n2_temp, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 0 | spots_col == 1) {
      n3 = NA
    } else {
      n3_temp = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (length(n3_temp) == 0) {
        n3 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n3_temp] == 0 | !(n3_temp %in% rownames(spots_programs_comp))){
        n3 = NA
      } else {
        n3 = as.character(spots_programs_comp[n3_temp, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 126 | spots_col == 127) {
      n4 = NA
    } else {
      n4_temp = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (length(n4_temp) == 0) {
        n4 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n4_temp] == 0 | !(n4_temp %in% rownames(spots_programs_comp))){
        n4 = NA
      } else {
        n4 = as.character(spots_programs_comp[n4_temp, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 0 | spots_row == 77) {
      n5 = NA
    } else {
      n5_temp = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (length(n5_temp) == 0) {
        n5 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n5_temp] == 0 | !(n5_temp %in% rownames(spots_programs_comp))){
        n5 = NA
      } else {
        n5 = as.character(spots_programs_comp[n5_temp, -ncol(spots_programs_comp)])
      }
    }
    if (spots_col == 127 | spots_row == 77) {
      n6 = NA
    } else {
      n6_temp = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (length(n6_temp) == 0) {
        n6 = NA
      } else if (spots_positions$V2[spots_positions$V1 == n6_temp] == 0 | !(n6_temp %in% rownames(spots_programs_comp))){
        n6 = NA
      } else {
        n6 = as.character(spots_programs_comp[n6_temp, -ncol(spots_programs_comp)])
      }
    }
    neighbors_table[[spot]] = data.frame(n1,n2,n3,n4,n5,n6)
  }
  neighbors_table <- do.call(rbind, neighbors_table)
  return(neighbors_table)
}


# ---- plotting helper functions ----

program_coherence_per_sample_data = function() {
  all_coherence_scores = data.frame(matrix(nrow = length(sample_ids), ncol = ncol(metaprograms)))
  colnames(all_coherence_scores) = colnames(metaprograms); rownames(all_coherence_scores) = sample_ids
  for (i in 1:length(sample_ids)) {
    mean_prog_coh = colMeans(coherence[[i]], na.rm = T) # mean score of all spots of the MP in one sample
    rel_mean_prog_coh = mean_prog_coh / colMeans(all_spots_programs_comp_norm[[i]][, 1:length(mp_names)], na.rm = T) # adjusted for program abundance
    all_coherence_scores[i, ] = rel_mean_prog_coh 
  }
  data = pivot_longer(all_coherence_scores, cols=everything(), names_to="metaprogram", values_to="score")
  data$sample = rep(sample_ids, each = ncol(all_coherence_scores)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
  average_scores = data %>% group_by(sample) %>% dplyr::summarize(average_score = mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE))
  data$average_score = rep(average_scores[match(unique(data$sample), average_scores$sample), ]$average_score, each = ncol(all_coherence_scores))
  data$sd = rep(average_scores$sd, each = ncol(all_coherence_scores))
  data = data %>% arrange(dplyr::desc(average_score)); data$sample = factor(data$sample, levels = unique(data$sample))
  return(data)
}

program_complexity_per_sample_data = function() {
  complexity = data.frame(matrix(nrow = length(sample_ids), ncol = ncol(metaprograms)))
  colnames(complexity) = colnames(metaprograms); rownames(complexity) = sample_ids
  for (i in 1:length(sample_ids)) {
    mean_prog_coh = colMeans(complexity_win_combined[[i]], na.rm = T) # mean score of all spots of the MP in one sample
    rel_mean_prog_coh = mean_prog_coh #* colMeans(all_spots_norm_comp[[i]], na.rm = T) * 10 # adjusted for program abundance
    complexity[i, ] = rel_mean_prog_coh 
  }
  data = pivot_longer(complexity, cols=everything(), names_to="metaprogram", values_to="score")
  data$sample = rep(sample_ids, each = ncol(complexity)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
  average_scores = data %>% group_by(sample) %>% dplyr::summarize(average_score = mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE))
  data$average_score = rep(average_scores[match(unique(data$sample), average_scores$sample), ]$average_score, each = ncol(complexity))
  data$sd = rep(average_scores$sd, each = ncol(complexity))
  data = data %>% arrange(dplyr::desc(average_score)); data$sample = factor(data$sample, levels = unique(data$sample))
  return(data)
}


plot_score_per_sample = function(data, split.by = NULL, anno = NULL, anno_name = "", order_manual = NULL, ylab = "score") {
  if (is.null(split.by)) split.by = rep("all", length(sample_ids)); names(split.by) = sample_ids # return(plot_score_per_sample_simple(data, ylab))
  # mp_cols = my_cols[1:ncol(metaprograms)]; names(mp_cols) = names(metaprograms)
  y_max <- max(data$score, na.rm = TRUE); y_min = min(data$score, na.rm = TRUE)
  group_sizes = table(split.by); group_sizes = group_sizes[order(group_sizes, decreasing = T)]
  data$group = split.by[match(data$sample, names(split.by))]; data$group = factor(data$group, levels = names(group_sizes))
  if (!is.null(anno)) data$anno = anno[match(data$sample, names(anno))]
  if (!is.null(order_manual)) data$sample = factor(data$sample, levels = order_manual)
  plots <- lapply(levels(data$group), function(group) {
    data = data[data$group == group, ]
    p = ggplot(data, aes(x = sample, y = score, fill = metaprogram)) + ylim(y_min, y_max+0.025) +
      geom_line(aes(x = sample, y = average_score, group = 1), size=0.2, color="grey90") +
      geom_errorbar(aes(ymin=average_score-sd, ymax=average_score+sd), size = 0.3, width = 0.1, color = "black") +
      geom_point(aes(x = sample, y = average_score), color = "black", size = 2) + # points for average score of sample
      geom_point(aes(fill = metaprogram), size = 5, shape = 21, color = "black") + # points for metaprograms
      scale_fill_manual(values = st_mp_cols) + labs(title = group, x = NULL, y = ylab, fill = "Program") + theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.border = element_rect(color = "black", fill = NA, size = 0.5))
    if (!is.null(anno)) { p = p + ggnewscale::new_scale_fill() + geom_tile(aes(x = sample, y = y_max+0.02, fill = anno), height = 0.01) +
        scale_fill_gradientn(colors = c("#fff6f2", "#e6489b", "#5c006f"), name = anno_name, limits = c(min(anno), max(anno))) }
    p
  })
  if (length(unique(split.by)) > 1) {
    plots[1:(length(plots)-1)] = lapply(plots[1:(length(plots)-1)], function(p) p + NoLegend())
    plots[2:length(plots)] =  lapply(plots[2:length(plots)], function(p) p + labs(y = NULL) + theme(axis.text.y = element_blank()))
    relative_widths <- group_sizes / sum(group_sizes) # Dynamically calculate subplot widths based on the number of rows in each group
    return(wrap_plots(plots, widths = relative_widths))
  } else {
    return(plots[[1]])
  }
}

plot_mean_score = function(scores, score_name, separate_wins = F, adjust_by_abund = F, split.by = NULL, add_trendline = T) {
  if (is.null(split.by)) split.by = rep("all", length(sample_ids))
  # --- 1. get data and put into ggplot shape ---
  if (separate_wins) { # use scores ='complexity' here
    data = bind_rows(lapply(unique(split.by), function(group) {
      avg_complexity = data.frame(matrix(nrow = length(win_sizes), ncol = ncol(metaprograms)))
      colnames(avg_complexity) = colnames(metaprograms); rownames(avg_complexity) = win_sizes; sd_avg_complexity = avg_complexity
      for (i in 1:length(win_sizes)) {
        mean_win_complexity = lapply(scores[split.by %in% group], function(s) colMeans(s[[i]], na.rm = T))
        avg_complexity[i, ] = sapply(1:length(mean_win_complexity[[1]]), function(j) { mean(sapply(mean_win_complexity, "[", j), na.rm = T) }) # mean across samples
        sd_avg_complexity[i, ] = sapply(1:length(mean_win_complexity[[1]]), function(j) { sd(sapply(mean_win_complexity, "[", j), na.rm = T) }) # standard deviation across samples
      }
      long_scores <- pivot_longer(avg_complexity, cols=everything(), names_to="metaprogram", values_to="score")
      long_sd <- pivot_longer(sd_avg_complexity, cols=everything(), names_to="metaprogram", values_to="sd")
      data = cbind(long_scores, sd = long_sd$sd); data$win_size = as.factor(rep(win_sizes, each = ncol(avg_complexity)))
      average_scores <- data %>% group_by(metaprogram) %>% dplyr::summarize(average_score = mean(score, na.rm = TRUE))
      data$average_score = rep(average_scores[match(unique(data$metaprogram), average_scores$metaprogram), ]$average_score, length(win_sizes))
      # data = data %>% arrange(dplyr::desc(average_score)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
      data$group = group
      return(data)
    }))
  } else { # use scores = 'complexity_win_combined' or 'coherence here
    mean_prog_score = lapply(1:length(sample_ids), function(i) {
      mean_prog_score = colMeans(scores[[i]], na.rm = T)
      return(if (adjust_by_abund) mean_prog_score/colMeans(all_spots_programs_comp_norm[[i]][, 1:length(mp_names)], na.rm = T) else mean_prog_score) # (optional) adjust for average prog abundance in the sample, scale into more handy range (x10)
    })
    names(mean_prog_score) = split.by
    data = bind_rows(lapply(unique(split.by), function(group) {
      mean_coherence = sapply(1:length(mp_names), function(i) { mean(sapply(mean_prog_score[names(mean_prog_score) == group], "[", i), na.rm = T) }) # mean across samples
      sd_coherence = sapply(1:length(mp_names), function(i) { sd(sapply(mean_prog_score[names(mean_prog_score) == group], "[", i), na.rm = T) }) # standard deviation across samples
      data = data.frame(metaprogram = mp_names, score = mean_coherence, sd = sd_coherence, group = group)
      data = data %>% arrange(dplyr::desc(score)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
      return(data)
    }))
  }
  # --- 2. do the plot(s) ---
  y_max <- max(data$score, na.rm = TRUE) + data$sd[data$score == max(data$score, na.rm = TRUE)]
  y_min = max(0, min(data$score, na.rm = TRUE) - data$sd[data$score == min(data$score, na.rm = TRUE)])
  plot_order = if (!is.null(levels(split.by))) levels(split.by) else unique(split.by)
  plots <- lapply(plot_order, function(group) {
    data = data[data$group == group, ]
    data = data %>% arrange(dplyr::desc(score)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
    if (separate_wins) {
      data = data %>% arrange(dplyr::desc(average_score)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
      p = ggplot(data, aes(x = metaprogram, y = score, color = win_size, group = win_size)) + labs(x = NULL, y = paste("mean", score_name), title = group) + 
        geom_line(aes(x = metaprogram, y = average_score, group = 1), color="grey30") + scale_color_manual(values = c("#aebfd1", "#055b96", "#011726"))
    } else {
      p = ggplot(data, aes(x = metaprogram, y = score)) + labs(x = NULL, y = if (adjust_by_abund) paste("mean", score_name, "/ abundance") else paste("mean", score_name), title = group) +
        if (add_trendline) geom_line(aes(x = metaprogram, y = score, group = 1), size = 1, color="#003967")
    }
    p = p + geom_point(position = position_dodge(width = 0.4), size = 2) + theme_minimal() + coord_cartesian(ylim = c(y_min, y_max)) +
      geom_errorbar(aes(ymin=score-sd, ymax=score+sd), size = 0.3, width = 0.2, position = position_dodge(width = 0.4)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.border = element_rect(color = "black", fill = NA, size = 0.5))
  })
  if (length(unique(split.by)) > 1) {
    plots[1:(length(plots)-1)] = lapply(plots[1:(length(plots)-1)], function(p) p + NoLegend())
    plots[2:length(plots)] =  lapply(plots[2:length(plots)], function(p) p + labs(y = NULL) + theme(axis.text.y = element_blank()))
    return(wrap_plots(plots, nrow = 1))
  } else {
    return(plots[[1]])
  }
}

plot_mean_score_comparative = function(scores, split.by, pair = NULL, score_name = "complexity", order_by_difference = T, cols = NULL) {
  if (is.null(pair)) pair = split.by
  # --- 1. get data and put into ggplot shape ---
  mean_prog_score = lapply(1:length(sample_ids), function(i) colMeans(scores[[i]], na.rm = T))
  names(mean_prog_score) = split.by
  data = bind_rows(lapply(unique(split.by), function(group) {
    mean_coherence = sapply(1:length(mp_names), function(i) { mean(sapply(mean_prog_score[names(mean_prog_score) == group], "[", i), na.rm = T) }) # mean across samples
    sd_coherence = sapply(1:length(mp_names), function(i) { sd(sapply(mean_prog_score[names(mean_prog_score) == group], "[", i), na.rm = T) }) # standard deviation across samples
    data = data.frame(metaprogram = mp_names, score = mean_coherence, sd = sd_coherence, group = group)
    data = data %>% arrange(dplyr::desc(score)); data$metaprogram = factor(data$metaprogram, levels = unique(data$metaprogram))
    return(data)
  }))
  data = data[data$group %in% pair, ]
  # --- 2. do the plot(s) ---
  y_max <- max(data$score, na.rm = TRUE) + data$sd[data$score == max(data$score, na.rm = TRUE)]
  y_min = max(0, min(data$score, na.rm = TRUE) - data$sd[data$score == min(data$score, na.rm = TRUE)])
  plot_order = if (!is.null(levels(split.by))) levels(split.by) else unique(split.by)
  if (order_by_difference) {
    diff_order = sapply(unique(data$metaprogram), function(mp) {x = data[data$metaprogram == mp, "score"]; diff = x[1] - x[2]; names(diff) = mp; diff }); diff_order = diff_order[order(diff_order, decreasing = T)]
    data$metaprogram = factor(data$metaprogram, levels = names(diff_order))
  }
  p = ggplot(data, aes(x = metaprogram, y = score, color = group, group = group)) + labs(x = NULL, y = paste("Mean", score_name), title = "pairwise comparison") + # geom_smooth(method = "lm", se = FALSE) +
    geom_line(size = 1, position = position_dodge(width = 0.2)) + # scale_linetype_manual(values = c("solid", "dashed")) +
    geom_point(position = position_dodge(width = 0.2), size = 2) + theme_minimal() + coord_cartesian(ylim = c(y_min, y_max)) + 
    geom_errorbar(aes(ymin=score-sd, ymax=score+sd), size = 0.3, width = 0.2, position = position_dodge(width = 0.2)) +
    scale_color_manual(values = cols) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.border = element_rect(color = "black", fill = NA, size = 0.5))
  return(p)
}


plot_mp_zone_abundance = function(mp_zone_abund, zone_cols, ylim = c(1, 1), split.by = NULL, group.by = NULL, normalise = T) {
  if (is.null(split.by)) split.by = rep("all", length(sample_ids))
  plots = lapply(unique(split.by), function(group) {
    data = mp_zone_abund[split.by == group]
    count_sums = data.frame(n = colSums(Reduce("+", data))); count_sums$MP = rownames(count_sums); count_sums$n = paste0("n=", count_sums$n)
    data = as.data.frame(Reduce("+", data) / length(data))
    if (normalise) data = apply(data, 2, function(mp) mp / sum(abs(mp)))
    data[c(3,4,5), ] = -data[c(3,4,5), ]# ; data = rbind(data, intermediate_neg = -data[3, ]/2); data[3, ] = data[3, ]/2
    # data = Reduce("+", lapply(data, function(x) replace(x, is.na(x), 0))) / Reduce("+", lapply(data, function(x) !is.na(x)))
    data = as.data.frame(data); data$zone = factor(rownames(data), levels = c("structured_high", "structured_low", "disorganised_high", "disorganised_low", "intermediate")) # unique(rownames(data)))
    data <- data %>% pivot_longer(cols = -zone, names_to = "MP", values_to = "mean_expr")
    # data$count_sums = sapply(data$MP, function(mp) count_sums[mp])
    if (!is.null(group.by)) {
      data = bind_rows(lapply(names(group.by), function(name) {
        data %>% filter(MP %in% group.by[[name]]) %>% group_by(zone) %>% dplyr::summarise(mean_expr = mean(mean_expr), .groups = "drop") %>% mutate(MP = name)
      }))
      count_sums = bind_rows(lapply(names(group.by), function(name) {
        count_sums %>% filter(MP %in% group.by[[name]]) %>% dplyr::summarise(n = paste0("n=", sum(as.integer(sub("n=", "", n)))), .groups = "drop") %>% mutate(MP = name)
      }))
    }
    st_sums = sapply(unique(data$MP), function(mp) sum(data[data$zone %in% c("structured_high", "structured_low") & data$MP == mp, ]$mean_expr))
    data$MP = factor(data$MP, levels = names(st_sums[order(st_sums, decreasing = T)]))
    ggplot(data, aes(x = MP, y = mean_expr, fill = zone)) + geom_bar(stat = "identity", position = "stack") + geom_hline(yintercept=0, color="black") +
      labs(x = NULL, y = "relative abundance", title = group) + scale_fill_manual(values = zone_cols) + 
      geom_text(data = count_sums, aes(x = MP, y = st_sums, label = n), inherit.aes = F, color = "black", angle = 90, hjust = -0.2, size = 2) +
      scale_y_continuous(limits = c(-ylim[1], ylim[2]), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1)) + theme_minimal() + #ylim(-ylim[1], ylim[2]) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid.major.x = element_blank())
  })
  if (length(unique(split.by)) > 1) {
    plots[1:(length(plots)-1)] = lapply(plots[1:(length(plots)-1)], function(p) p + theme(legend.position = "none"))
    plots[2:length(plots)] = lapply(plots[2:length(plots)], function(p) p + labs(y = NULL) + theme(axis.text.y = element_blank()))
  }
  return(plots)
}


                                    
# ------ PART 2: spatial association functions ------


# -----> run on server for higher win_size and rand number <------
calculate_association = function(s, mp_names, pairs, pairs_names, spots_positions, spots_programs_comp, out_dir, max_win_size = 15, rand_num = NULL) {
  message("calculating pairwise association across distances: 0-", max_win_size, " in ", s)

  # spots_programs_comp = all_spots_programs_comp[[s]]
  neighbors_table <- neighbors_table_funcV2(spots_positions, spots_programs_comp)
  
  all_cors = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size + 1)); all_pvals = all_cors
  for (w in 0:max_win_size) {
    peri_abund = get_perimeter_abund(spots_programs_comp, neighbors_table, mp_names, w)
    # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
    all_cor_res <- sapply(1:dim(pairs)[2], function(j) {
      cor_res <- cor.test(spots_programs_comp[, pairs[2, j]], peri_abund[, pairs[1, j]], method = "pearson") # <- betrachtet hier nur eine Richtung (nur cor(A,B_peri) und nicht cor(B, A_peri)), obwohl diese nicht gleich sind!
      return(cor_res)
    })
    all_cors[, w+1] = unlist(all_cor_res[4,]); all_pvals[, w+1] = unlist(all_cor_res[3,])
  }
  rownames(all_cors) <- pairs_names; rownames(all_pvals) = pairs_names
  all_association = list(all_cors, all_pvals)
  saveRDS(all_association, paste0(out_dir, "/", s, "_obs_association_w", max_win_size, ".rds"))
  
  if (!is.null(rand_num)) { # run the same for random positions
    message("calculating reference regional composition in ", s, " for ", rand_num, " randomly shuffled spot positions")
    all_rand <- lapply(1:rand_num, function(z) {
      new_pos_all <- sample(spots_positions$V1[spots_positions$V2 != 0], length(spots_positions$V1[spots_positions$V2 != 0]), replace = FALSE)
      new_spots_positions <- spots_positions
      new_spots_positions$V1[new_spots_positions$V2 != 0] <- new_pos_all

      neighbors_table <- neighbors_table_funcV2(new_spots_positions, spots_programs_comp)

      all_windows <- sapply(0:max_win_size, function(w) {
        association = get_perimeter_abund(spots_programs_comp, neighbors_table, mp_names, w)
        # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
        empty_win_spots = is.na(association[, 1])
        all_cor <- sapply(1:dim(pairs)[2], function(j) {
          a = spots_programs_comp[!empty_win_spots, pairs[2, j]]
          b = association[!empty_win_spots, pairs[1, j]]
          return(cor(a, b))
        })
        all_cor <- data.frame(pair_cors = all_cor); row.names(all_cor) <- pairs_names
        return(all_cor)
      })
      sample_association <- as.data.frame(all_windows); row.names(sample_association) <- pairs_names
      return(sample_association)
    })
    sample_mean_rand_prox <- Reduce("+", all_rand) / length(all_rand)
    sample_sd_rand_prox <- round(apply(array(unlist(all_rand), c(length(pairs_names), max_win_size, rand_num)), c(1,2), sd),4)
    message("saving reference association in ", s, " (", nrow(spots_programs_comp), " spots)")
    saveRDS(list(sample_mean_rand_prox, sample_sd_rand_prox), paste0(out_dir, "/", s, "_reference_regional_association_w", max_win_size, "_rand", rand_num, ".rds"))
    message("finished calculating reference association in ", s, " (", nrow(spots_programs_comp), " spots)")
  }
}


calculate_association_by_ratio = function(s) { # -----> run on server for higher win_size and rand number <------
  message("calculating pairwise association in regional composition in ", s, " using 'max_win_size' = ", max_win_size)
  
  all_spot_progs = all_spots_programs_comp[[s]]
  neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
  all_associations = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size + 1))
  for (w in 0:max_win_size) {
    association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
    # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
    all_associations[, w+1] = sapply(1:dim(pairs)[2], function(j) {
      spot_prog_a = all_spot_progs[, pairs[2, j]]; neighbours_prog_b = association[, pairs[1, j]]
      association_by_ratio = mean(sapply(1:length(spot_prog_a), function(i) {
        ratio = spot_prog_a[i] / neighbours_prog_b[i]; association = (1 - abs((ratio - 1) / (ratio + 1)))
        return(association)
      }), na.rm = T)
      return(association_by_ratio)
    })
  }
  rownames(all_associations) <- pairs_names
  saveRDS(all_associations, paste0(out_dir, s, "_obs_association_by_ratio_w", max_win_size, ".rds"))
  
  # run the same for random positions 
  # message("calculating reference regional composition in ", s, " for ", rand_num, " randomly shuffled spot positions")
  # all_rand <- lapply(1:rand_num, function(z) {
  #   new_pos_all <- sample(all_spots_positions[[s]]$V1[all_spots_positions[[s]]$V2 != 0], length(all_spots_positions[[s]]$V1[all_spots_positions[[s]]$V2 != 0]), replace = FALSE)
  #   spots_positions <- all_spots_positions[[s]]
  #   spots_positions$V1[spots_positions$V2 != 0] <- new_pos_all
  #   
  #   neighbors_table <- neighbors_table_funcV2(spots_positions, all_spot_progs)
  #   
  #   all_windows <- sapply(0:max_win_size, function(w) {
  #     association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
  #     # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
  #     empty_win_spots = is.na(association[, 1])
  #     all_associations <- sapply(1:dim(pairs)[2], function(j) {
  #       spot_prog_a = all_spot_progs[, pairs[2, j]]; neighbours_prog_b = association[, pairs[1, j]]
  #       association_by_ratio = mean(sapply(1:length(spot_prog_a), function(i) {
  #         ratio = spot_prog_a[i] / neighbours_prog_b[i]; association = (1 - abs((ratio - 1) / (ratio + 1)))
  #         return(association)
  #       }), na.rm = T)
  #       return(weighted_association)
  #     })
  #     all_associations <- data.frame(pair_cors = all_associations); row.names(all_associations) <- pairs_names
  #     return(all_associations)
  #   })
  #   sample_association <- as.data.frame(all_windows); row.names(sample_association) <- pairs_names
  #   return(sample_association)
  # })
  # sample_mean_rand_prox <- Reduce("+", all_rand) / length(all_rand)
  # sample_sd_rand_prox <- round(apply(array(unlist(all_rand), c(length(pairs_names), max_win_size, rand_num)), c(1,2), sd),4)
  # message("saving reference association in ", s, " (", nrow(all_spot_progs), " spots)")
  # saveRDS(list(sample_mean_rand_prox, sample_sd_rand_prox), paste0(out_dir, s, "_reference_regional_association_w", max_win_size, "_rand", rand_num, ".rds"))
  # message("finished calculating reference association in ", s, " (", nrow(all_spot_progs), " spots)")
}

# calculates (only) the average gradients per pair of programs and per distance
calculate_pairwise_gradients = function(s, gradient = 1) {
  message("calculating ", if (gradient == 1) "first" else "second", " gradients in ", s, " using 'max_win_size' = ", max_win_size)
  all_spot_progs = all_spots_programs_comp[[s]]
  neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
  all_cors = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size+1))
  for (w in 0:max_win_size) {
    # association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
    gradients = if (gradient == 1) prog_gradients_win_combined[[s]] else derivatives[[s]] # chooses between first or second derivative
    neighbour_gradients = find_gradients(all_spot_progs, neighbors_table, mp_names, w, gradients)
    # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
    all_res <- sapply(1:dim(pairs)[2], function(j) {
      # spot_prog_a = all_spot_progs[, pairs[2, j]]; neighbours_prog_b = association[, pairs[1, j]]
      spot_gradient_a = gradients[, pairs[2, j]]; neighbours_gradient_b = neighbour_gradients[, pairs[1, j]]
      gradients = mean(sapply(1:length(spot_gradient_a), function(i) { spot_gradient_a[i] * neighbours_gradient_b[i] }), na.rm = T)
      return(gradients)
    })
    all_cors[, w+1] = all_res
  }
  rownames(all_cors) <- pairs_names
  return(all_cors)
}


# calculates the gradient-weighted association by correlation
calculate_gradient_weighted_association = function(s, gradient = 1) {
  message("calculating ", if (gradient == 1) "first" else "second", " gradient-weighted, cor-based association in ", s, " using 'max_win_size' = ", max_win_size)
  all_spot_progs = all_spots_programs_comp[[s]]
  neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
  all_cors = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size + 1)); all_pvals = all_cors
  for (w in 0:max_win_size) {
    association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
    gradients = if (gradient == 1) prog_gradients_win_combined[[s]] else derivatives[[s]] # chooses between first or second derivative
    # neighbour_gradients = find_gradients(all_spot_progs, neighbors_table, mp_names, w, gradients)
    
    # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
    all_cor_res <- sapply(1:dim(pairs)[2], function(j) {
      spot_abund = all_spot_progs[, pairs[1, j]]; neighbours_abund = association[, pairs[2, j]]; weights = gradients[, pairs[1, j]]
      cor_res1 <- weighted_cor(neighbours_abund, spot_abund, weights)
      
      spot_abund = all_spot_progs[, pairs[2, j]]; neighbours_abund = association[, pairs[1, j]]; weights = gradients[, pairs[2, j]]
      cor_res2 <- weighted_cor(neighbours_abund, spot_abund, weights)
      
      cor_res = (cor_res1 + cor_res2) / 2 # get average to determine the bi-directional correlation
      return(cor_res)
    })
    all_cors[, w+1] = unlist(all_cor_res[1,]); all_pvals[, w+1] = unlist(all_cor_res[2,])
  }
  rownames(all_cors) <- pairs_names; rownames(all_pvals) = pairs_names
  return(list(all_cors, all_pvals))
}


# modified, >> weighted << Pearson's correlation
weighted_cor <- function(a, b, c) {
  valid_indices <- complete.cases(a, b, c)
  a <- a[valid_indices]; b <- b[valid_indices]; c <- c[valid_indices]
  sum_c <- sum(c)
  mean_a <- sum(c * a) / sum_c; mean_b <- sum(c * b) / sum_c # Gewichtete Mittelwerte
  var_a <- sum(c * (a - mean_a)^2) / sum_c; var_b <- sum(c * (b - mean_b)^2) / sum_c # Gewichtete Varianzen
  cov_ab <- sum(c * (a - mean_a) * (b - mean_b)) / sum_c # Gewichtete Kovarianz
  r <- cov_ab / sqrt(var_a * var_b) # Gewichtete Korrelation
  
  n_eff <- (sum(c)^2) / sum(c^2) # Effektive Stichprobengröße
  t_value <- r * sqrt((n_eff - 2) / (1 - r^2)) # Berechnung des t-Werts
  p_value <- 2 * pt(-abs(t_value), df = n_eff - 2) # Berechnung des p-Werts
  
  return(c(correlation = r, p_value = p_value))
}



# calculates the gradient-weighted association by ratio (1st implementation)
# calculate_gradient_weighted_association_by_ratio_1 = function(s, gradient = 1) {
#   message("calculating ", if (gradient == 1) "first" else "second", " gradient-based association in ", s, " using 'max_win_size' = ", max_win_size)
#   all_spot_progs = all_spots_programs_comp[[s]]
#   neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
#   all_cors = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size+1))
#   for (w in 0:max_win_size) {
#     association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
#     gradients = if (gradient == 1) prog_gradients_win_combined[[s]] else derivatives[[s]] # chooses between first or second derivative
#     neighbour_gradients = find_gradients(all_spot_progs, neighbors_table, mp_names, w, gradients)
#     # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
#     all_res <- sapply(1:dim(pairs)[2], function(j) {
#       spot_prog_a = all_spot_progs[, pairs[2, j]]; neighbours_prog_b = association[, pairs[1, j]]
#       spot_gradient_a = gradients[, pairs[2, j]]; neighbours_gradient_b = neighbour_gradients[, pairs[1, j]]
#       weighted_association = mean(sapply(1:length(spot_prog_a), function(i) {
#         ratio = spot_prog_a[i] / neighbours_prog_b[i]
#         association = (1 - abs((ratio - 1) / (ratio + 1))) * spot_gradient_a[i] * neighbours_gradient_b[i]
#         return(association)
#       }), na.rm = T)
#       return(weighted_association)
#     })
#     all_cors[, w+1] = all_res
#   }
#   rownames(all_cors) <- pairs_names
#   return(all_cors)
# }

# calculates the gradient-weighted association by ratio (2nd implementation)
# calculate_gradient_weighted_association_by_ratio_2 = function(s, gradient = 1) {
#   message("calculating gradient weighted association by ratio in ", s, " using 'max_win_size' = ", max_win_size)
#   all_spot_progs = all_spots_programs_comp[[s]]
#   neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
#   all_associations = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size))
#   for (w in 0:max_win_size) {
#     association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
#     gradients = if (gradient == 1) prog_gradients_win_combined[[s]] else derivatives[[s]]
#     
#     # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
#     associations <- sapply(1:dim(pairs)[2], function(j) {
#       spot_prog_a = all_spot_progs[, pairs[1, j]]; neighbours_prog_a = association[, pairs[1, j]]; weights_a = gradients[, pairs[1, j]]
#       spot_prog_b = all_spot_progs[, pairs[2, j]]; neighbours_prog_b = association[, pairs[2, j]]; weights_b = gradients[, pairs[2, j]]
#       
#       ratio = spot_prog_a / neighbours_prog_b
#       obs_associations = 1 - abs((ratio - 1) / (ratio + 1))
#       association = weighted_unadjusted_mean(obs_associations, weights_a)
#       
#       ratio_rev = spot_prog_b / neighbours_prog_a
#       obs_associations_rev = 1 - abs((ratio_rev - 1) / (ratio_rev + 1))
#       association_rev = weighted_unadjusted_mean(obs_associations_rev, weights_b)
# 
#       return((association + association_rev) / 2) # average of forward and reverse association
#     })
#     all_associations[, w+1] = associations
#   }
#   rownames(all_associations) <- pairs_names
#   return(all_associations)
# }


# j = 29

# prog_a = "Cell Cycle"; prog_b = "Undifferentiated"
# 
# prog_a = "Ependymal"; prog_b = "Undifferentiated"
# a = all_spot_progs[, prog_a]; neighbours_prog_a = association[, prog_a]
# spot_prog_b = all_spot_progs[, prog_b]; b = association[, prog_b]
# 
# ratio_a = a / neighbours_prog_a; w_a = abs((ratio_a - 1) / (ratio_a + 1))
# ratio_b = spot_prog_b / b; w_b = abs((ratio_b - 1) / (ratio_b + 1))
# 
# weighted_associations = sapply(1:length(a), function(i) {
#   ratio = a[i] / b[i]
#   association = (1 - abs((ratio - 1) / (ratio + 1))) * (w_a[i] * w_b[i])
#   return(association)
# })
# max(weighted_associations, na.rm = T)

# spot = rownames(all_spot_progs)[13]
# SpatialPlot(sstobj, cells.highlight = "GATCCCTTTATACTGC-1", pt.size.factor = 3, image.alpha = 0)
# 
# SpatialPlot(sstobj, cells.highlight = win_spots, pt.size.factor = 3, image.alpha = 0)

get_perimeter_abund = function(all_spot_progs, neighbors_table, mp_names, w) {
  # determines the average abundance of each program in the perimeter 'band' of radius w, for each spot
  if (w == 0) return(all_spot_progs[, 1:length(mp_names)]) # edge case: co-localisation
  abund <- t(sapply(rownames(all_spot_progs), function(spot) {
    all_win_spots = spot
    sapply(1:w, function(i) {
      if (i == w) win_spots <<- setdiff(unique(na.omit(as.character(neighbors_table[all_win_spots,]))), all_win_spots)
      all_win_spots <<- unique(c(all_win_spots, unique(na.omit(as.character(neighbors_table[all_win_spots,])))))
    })
    win_abund = colMeans(all_spot_progs[win_spots, 1:length(mp_names)], na.rm = T)
    return(win_abund)
  }))
  return(abund)
}

find_gradients = function(all_spot_progs, neighbors_table, mp_names, w, all_gradients) {
  # determines the average gradient of each program in the perimeter 'band' of radius w, for each spot
  # if (w == 0) return(all_spot_progs[, 1:length(mp_names)]) # edge case co-localisation
  association <- t(sapply(rownames(all_spot_progs), function(spot) {
    all_win_spots = spot
    sapply(1:w, function(i) {
      if (i == w) win_spots <<- setdiff(unique(na.omit(as.character(neighbors_table[all_win_spots,]))), all_win_spots)
      all_win_spots <<- unique(c(all_win_spots, unique(na.omit(as.character(neighbors_table[all_win_spots,])))))
    })
    win_abund = colMeans(all_gradients[win_spots, 1:length(mp_names)], na.rm = T)
    return(win_abund)
  }))
  return(association)
}

calculate_derivate = function(s, mp_names) {
  message("calculating derivative of gradients in ", s)
  spots_filt = all_spots_positions[[s]][all_spots_positions[[s]]$V1 %in% rownames(all_spots_programs_comp[[s]]), ]
  prog_derivatives = prog_gradients_win_combined[[s]]
  for (program in mp_names) {
    x_offset = min(spots_filt$V3) - 1; y_offset = min(spots_filt$V4) -1
    mtx = matrix(nrow = max(spots_filt$V3) - x_offset, ncol = max(spots_filt$V4) - y_offset)
    for (spot in rownames(spots_filt)) {
      x = spots_filt[spot, "V3"] - x_offset; y = spots_filt[spot, "V4"] - y_offset
      mtx[x, y] = prog_gradients_win_combined[[s]][spot, program]
    }
    mtx_filled = fill_adjacent_na(mtx)
    dx <- diff(mtx_filled, differences = 1, lag = 1); dy <- t(diff(t(mtx_filled), differences = 1, lag = 1))
    dx <- rbind(dx, rep(0, ncol(dx))); dy <- cbind(dy, rep(0, nrow(dy))) # append columns of zeros to match dimensions
    gradient_magnitude <- sqrt(dx^2 + dy^2)

    existing_spot_gradients = spots_filt[, c(1,3,4)]
    for (spot in rownames(existing_spot_gradients)) {
      x = existing_spot_gradients[spot, 2] - x_offset; y = existing_spot_gradients[spot, 3] - y_offset
      existing_spot_gradients[spot, 1] = gradient_magnitude[x, y]
    }
    # correct for edge cases scoring abnormally high
    derivative = as.double(existing_spot_gradients[order(rownames(existing_spot_gradients)), ]$V1)
    p95 <- quantile(derivative, 0.95)
    derivative[derivative > p95] <- p95 # replace values greater than the 95th percentile with the 95th percentile
    prog_derivatives[, program] = derivative
  }
  return(prog_derivatives)
}

fill_adjacent_na <- function(mtx) {
  mtx_interpolated <- mtx
  for (i in 1:nrow(mtx)) {
    for (j in 1:ncol(mtx)) {
      if (is.na(mtx[i, j])) next
      neighbors <- list(if (i > 1) c(i - 1, j) else NA, if (i < nrow(mtx)) c(i + 1, j) else NA,
                        if (j > 1) c(i, j - 1) else NA, if (j < ncol(mtx)) c(i, j + 1) else NA)
      for (n in neighbors[!is.na(neighbors)]) {
        if (is.na(mtx[n[1], n[2]])) {
          na_neighbors = c(if (n[1] > 1) mtx[n[1] - 1, n[2]] else NA, if (n[1] < nrow(mtx)) mtx[n[1] + 1, n[2]] else NA,
                           if (n[2] > 1) mtx[n[1], n[2] - 1] else NA, if (n[2] < ncol(mtx)) mtx[n[1], n[2] + 1] else NA)
          mtx_interpolated[n[1], n[2]] <- mean(na_neighbors, na.rm = TRUE)
        }
      }
    }
  }
  return(mtx_interpolated)
}

# calculate_exclusive_gradient_association = function(s, gradient = 1) {
#   message("calculating exclusive ", if (gradient == 1) "first" else "second", " gradient-based association in ", s, " using 'max_win_size' = ", max_win_size)
#   all_spot_progs = all_spots_programs_comp[[s]]
#   neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
#   all_cors = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size))
#   for (w in 1:max_win_size) {
#     association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
#     gradients = if (gradient == 1) prog_gradients_win_combined[[s]] else derivatives[[s]] # chooses between first or second derivative
#     neighbour_gradients = find_gradients(all_spot_progs, neighbors_table, mp_names, w, gradients)
#     # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
#     all_cor_res <- sapply(1:dim(pairs)[2], function(j) {
#       spot_prog_a = all_spot_progs[, pairs[2, j]]; neighbours_prog_a = association[, pairs[2, j]]
#       spot_prog_b = all_spot_progs[, pairs[1, j]]; neighbours_prog_b = association[, pairs[1, j]]
#       
#       # penalise non-exclusive associations by multiplying the observed correlations with exclusivity weights for program a and b
#       ratio_a = spot_prog_a / neighbours_prog_a; exclusivity_a = abs((ratio_a - 1) / (ratio_a + 1)) # if (w > 0) abs((ratio_a - 1) / (ratio_a + 1)) else 1
#       ratio_b = spot_prog_b / neighbours_prog_b; exclusivity_b = abs((ratio_b - 1) / (ratio_b + 1)) # if (w > 0) abs((ratio_b - 1) / (ratio_b + 1)) else 1
#       spot_gradient_a = gradients[, pairs[2, j]]
#       neighbours_gradient_b = neighbour_gradients[, pairs[1, j]]
#       
#       weighted_association = mean(sapply(1:length(spot_prog_a), function(i) {
#         ratio = spot_prog_a[i] / neighbours_prog_b[i]
#         association = (1 - abs((ratio - 1) / (ratio + 1))) * exclusivity_a[i] * exclusivity_b[i]  * 1000 * spot_gradient_a[i] * neighbours_gradient_b[i]
#         return(association)
#       }), na.rm = T)
#       return(weighted_association)
#     })
#     all_cors[, w] = all_cor_res
#   }
#   rownames(all_cors) <- pairs_names
#   return(all_cors)
# }
# 
# calculate_exclusivity = function(s) {
#   message("calculating exclusivity in ", s, " using 'max_win_size' = ", max_win_size)
#   all_spot_progs = all_spots_programs_comp[[s]]
#   neighbors_table <- neighbors_table_funcV2(all_spots_positions[[s]], all_spot_progs)
#   all_exclus = as.data.frame(matrix(nrow = length(pairs_names), ncol = max_win_size+1))
#   for (w in 0:max_win_size) {
#     association = get_perimeter_abund(all_spot_progs, neighbors_table, mp_names, w)
#     # calculate for every pair the correlation between all program A scores in centre spots and program B abundance in the window
#     all_exclusivities <- sapply(1:dim(pairs)[2], function(j) {
#       spot_prog_a = all_spot_progs[, pairs[2, j]]; neighbours_prog_a = association[, pairs[2, j]]
#       spot_prog_b = all_spot_progs[, pairs[1, j]]; neighbours_prog_b = association[, pairs[1, j]]
#       
#       # calculate exclusivity of program a as the ratio of a in centre vs. in neighbour spots (same for b)
#       ratio_a = spot_prog_a / neighbours_prog_a; exclusivity_a = abs((ratio_a - 1) / (ratio_a + 1)) # if (w > 0) abs((ratio_a - 1) / (ratio_a + 1)) else 1
#       ratio_b = spot_prog_b / neighbours_prog_b; exclusivity_b = abs((ratio_b - 1) / (ratio_b + 1)) # if (w > 0) abs((ratio_b - 1) / (ratio_b + 1)) else 1
#       exclusivity = mean(exclusivity_a * exclusivity_b, na.rm = T)
#       return(exclusivity)
#     })
#     all_exclus[, w+1] = all_exclusivities
#   }
#   rownames(all_exclus) <- pairs_names
#   return(all_exclus)
# }
# 
# calculate_derivate2 = function(s, mp_names) {
#   message("calculating derivative of gradients in ", s)
#   spots_filt = all_spots_positions[[s]][all_spots_positions[[s]]$V1 %in% rownames(all_spots_programs_comp[[s]]), ]
#   prog_derivatives = prog_gradients_win_combined[[s]]
#   for (program in mp_names) {
#     x_offset = min(spots_filt$V3) - 1; y_offset = min(spots_filt$V4) -1
#     mtx = matrix(nrow = max(spots_filt$V3) - x_offset, ncol = max(spots_filt$V4) - y_offset)
#     for (spot in rownames(spots_filt)) {
#       x = spots_filt[spot, "V3"] - x_offset; y = spots_filt[spot, "V4"] - y_offset
#       mtx[x, y] = prog_gradients_win_combined[[s]][spot, program]
#     }
#     mtx_filled = zoo::na.approx(mtx, rule = 1)
#     dx <- diff(mtx_filled, differences = 1, lag = 1); dy <- t(diff(t(mtx_filled), differences = 1, lag = 1))
#     dx <- rbind(dx, rep(0, ncol(dx))); dy <- cbind(dy, rep(0, nrow(dy))) # append columns of zeros to match dimensions
#     gradient_magnitude <- sqrt(dx^2 + dy^2)
#     gradient_magnitude[!is.na(mtx) & is.na(gradient_magnitude)] = "NaN"
#     gradients_filled = interpolate_matrix(gradient_magnitude)
#     
#     existing_spot_gradients = spots_filt[, c(1,3,4)]
#     for (spot in rownames(existing_spot_gradients)) {
#       x = existing_spot_gradients[spot, 2] - x_offset; y = existing_spot_gradients[spot, 3] - y_offset
#       existing_spot_gradients[spot, 1] = gradients_filled[x, y]
#     }
#     prog_derivatives[, program] = as.double(existing_spot_gradients[order(rownames(existing_spot_gradients)), ]$V1)
#   }
#   return(prog_derivatives)
# }
# 
# interpolate_matrix = function(mtx) {
#   mtx_interpolated <- mtx
#   for (i in 1:nrow(mtx)) {
#     for (j in 1:ncol(mtx)) {
#       if (!is.na(mtx[i, j]) && mtx[i, j] == "NaN") {
#         neighbors = c(if (i > 1) mtx[i - 1, j] else NA, if (i < nrow(mtx)) mtx[i + 1, j] else NA,
#                       if (j > 1) mtx[i, j - 1] else NA, if (j < ncol(mtx)) mtx[i, j + 1] else NA)
#         mtx_interpolated[i, j] = if (length(neighbors[!is.na(neighbors)]) > 0) mean(as.double(neighbors), na.rm = TRUE) else 0
#       }
#     }
#   }
#   return(mtx_interpolated)
# }


# finds all adjacent/neighbouring spots for each spot
neighbors_table_funcV2 = function(spots_positions, spots_programs) {
  neighbors_table = sapply(spots_programs$barcodes, function(spot) {
    spots_row = spots_positions[spots_positions$V1 == spot, 3]
    spots_col = spots_positions[spots_positions$V1 == spot, 4]
    if (spots_col == 0 | spots_row == 0) {
      c1 = NA
    } else {
      n1 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col - 1]
      if (length(n1) == 0 || (spots_positions$V2[spots_positions$V1 == n1] == 0 | !(n1 %in% spots_programs$barcodes))) {
        c1 = NA
      } else {
        c1 = as.character(n1)
      }
    }
    if (spots_col == 127 | spots_row == 0) {
      c2 = NA
    } else {
      n2 = spots_positions$V1[spots_positions$V3 == spots_row - 1 & spots_positions$V4 == spots_col + 1]
      if (length(n2) == 0 || (spots_positions$V2[spots_positions$V1 == n2] == 0 | !(n2 %in% spots_programs$barcodes))) {
        c2 = NA
      } else {
        c2 = as.character(n2)
      }
    }
    if (spots_col == 0 | spots_col == 1) {
      c3 = NA
    } else {
      n3 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col - 2]
      if (length(n3) == 0 || (spots_positions$V2[spots_positions$V1 == n3] == 0 | !(n3 %in% spots_programs$barcodes))) {
        c3 = NA
      } else {
        c3 = as.character(n3)
      }
    }
    if (spots_col == 126 | spots_col == 127) {
      c4 = NA
    } else {
      n4 = spots_positions$V1[spots_positions$V3 == spots_row & spots_positions$V4 == spots_col + 2]
      if (length(n4) == 0 || (spots_positions$V2[spots_positions$V1 == n4] == 0 | !(n4 %in% spots_programs$barcodes))) {
        c4 = NA
      } else {
        c4 = as.character(n4)
      }
    }
    if (spots_col == 0 | spots_row == 77) {
      c5 = NA
    } else {
      n5 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col - 1]
      if (length(n5) == 0 || (spots_positions$V2[spots_positions$V1 == n5] == 0 | !(n5 %in% spots_programs$barcodes))) {
        c5 = NA
      } else {
        c5 = as.character(n5)
      }
    }
    if (spots_col == 127 | spots_row == 77) {
      c6 = NA
    } else {
      n6 = spots_positions$V1[spots_positions$V3 == spots_row + 1 & spots_positions$V4 == spots_col + 1]
      if (length(n6) == 0 || (spots_positions$V2[spots_positions$V1 == n6] == 0 | !(n6 %in% spots_programs$barcodes))) {
        c6 = NA
      } else {
        c6 = as.character(n6)
      }
    }
    return(c(c1,c2,c3,c4,c5,c6))
  })
  neighbors_table = t(neighbors_table)
  row.names(neighbors_table) = spots_programs$barcodes
  return(neighbors_table)
}





# ----------------- plotting helper functions ----------------

plot_association_heatmap = function(mtx, order_mode = "association_type", title = NA, silent = F, asso_type_tresholds = c(0, 0), 
                                    anno = NULL, anno_cols = NULL, gaps_row = NULL, gaps_col = NULL, fixed_scale = F) {
  if (order_mode == "association_type") {
    types = c("coloc", rep("contact", 5), rep("dist_contact", 10))
    abs_max = apply(mtx, 1, function(pair) pair[which.max(abs(pair))])
    sign = ifelse(abs_max > asso_type_tresholds[1], "pos.", ifelse(abs_max < asso_type_tresholds[2], "neg.", "mixed"))   
    type = apply(mtx, 1, function(pair) types[which.max(abs(pair))])
    df = data.frame(A = rowMeans(mtx), B = factor(type, levels = unique(types)), C = factor(sign, levels = c("pos.", "mixed", "neg.")))
    df_ordered = df[order(df$C, df$B, df$A), ]
    for (type in unique(df_ordered[df_ordered$C == "pos.", ]$B)) { # reorder only pos. associations decreasingly for nicer looks
      df_sub = df_ordered[df_ordered$C == "pos." & df_ordered$B == type, ]
      df_ordered[df_ordered$C == "pos." & df_ordered$B == type, ] = df_sub[order(df_sub$A, decreasing = T), ]
      rownames(df_ordered)[df_ordered$C == "pos." & df_ordered$B == type] = rownames(df_sub[order(df_sub$A, decreasing = T), ])
    }
    mtx = mtx[rownames(df_ordered), ]
    anno = df_ordered[, c(2,3)]; rownames(anno) = rownames(df_ordered); colnames(anno) = c("association_type", "sign")
    gaps_row = which(diff(as.integer(anno$association_type)) != 0) # which(diff(as.integer(anno$sign)) != 0)
    anno_cols = list(association_type = c(coloc="#C43D3C", contact="#208cc0", dist_contact="#ffd60a"),
                     sign = c('pos.'="grey95", 'mixed'="grey75", 'neg.'="grey30")) #c('pos.'="#C43D3C", 'mixed'="grey93", 'neg.'="#246BAE"))
    anno_cols = lapply(anno_cols, function(col) col[names(col) %in% c(anno$association_type, anno$sign)])
  } else if (order_mode == "mean_score") {
    mtx = mtx[order(rowMeans(mtx), decreasing = T), ]
  } else if (order_mode == "colocalisation") {
    mtx = mtx[order(mtx[, 1], decreasing = T), ]
  }
  if (fixed_scale) breaks = seq(-1, 1, length.out = 101) else { max_val <- max(abs(mtx)); breaks <- seq(-max_val, max_val, length.out = 101) }
  return(pheatmap::pheatmap(mtx, cluster_cols = F, cluster_rows = F, main = title, angle_col = 0, border_color = NA, annotation_row = anno, silent = silent,
                     legend_labels = "association", show_row_names = T, color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)), breaks = breaks, labels_col = 0:15,
                     gaps_row = gaps_row, gaps_col = gaps_col, annotation_colors = anno_cols, annotation_names_row = F, annotation_legend = T))
}

plot_mp_association_mtx = function(dist, split.by = NULL, custom_order = NULL, nrow = 1) {
  if (is.null(split.by)) split.by = rep("all", length(sample_ids)); names(split.by) = sample_ids
  
  plots = lapply(unique(split.by), function(group) {
    dist_plots = lapply(unique(dist), function(dist) {
      samples = names(split.by)[split.by == group]
      avg_group_association <- Reduce(function(x, y) Map(`+`, x, y), association[samples])
      avg_group_association <- as.data.frame(lapply(avg_group_association, function(x) x / length(samples)))
      rownames(avg_group_association) = pairs_names; colnames(avg_group_association) = 0:15
      data_unscaled <- reshape2::melt(avg_group_association); data_unscaled$pair_name = rep(rownames(avg_group_association), ncol(avg_group_association))
      names(data_unscaled)[1] = "radius"
      df = data_unscaled[data_unscaled$radius == dist, ]
      df$mp1 = sub(" x .*", "", df$pair_name); df$mp2 = sub("^.*? x ", "", df$pair_name); df = df[, c(2,4,5)] # split the pairs into two columns
      df_rev = df; df_rev$mp1 = df$mp2; df_rev$mp2 = df$mp1; df_diag = data.frame(value = rep(1, length(mp_names)), mp1 = mp_names, mp2 = mp_names)
      df = rbind(df, df_rev, df_diag)
      if (!is.null(custom_order)) {
        levels = if (!is.list(custom_order)) custom_order else if (group %in% names(custom_order)) custom_order[[group]] else unique(df$mp1)
        df$mp1 = factor(df$mp1, levels = levels); df$mp2 = factor(df$mp2, levels = levels) 
      }
      ggplot(df, aes(x = mp1, y = mp2, fill = value, size = abs(value))) + geom_point(shape = 21, color = "white", stroke = NA) + geom_tile(color = "grey70", fill = NA, size = 0.25) +
        scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1, 1)) + scale_size(range = c(0, 7.5)) +  
        coord_fixed() + labs(title = paste0(group, " (dist = ", dist, ")"), fill = "Correlation") + guides(size = "none") + theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(), panel.grid = element_blank())
    })
    dist_plots
  })
  plots = purrr::flatten(plots)
  if (length(plots) > 1) {
    plots[1:(length(plots)-1)] = lapply(plots[1:(length(plots)-1)], function(p) p + theme(legend.position = "none"))
    wrap_plots(plots, nrow = nrow)
  } else plots[[1]]
}

scale_min_max_5_95 <- function(x) {
  x_min <- quantile(x, 0.05, na.rm = TRUE); x_max <- quantile(x, 0.95, na.rm = TRUE) # Calculate the 5th and 95th percentiles
  x_scaled <- -1 + 2 * ((x - x_min) / (x_max - x_min))
  x_scaled[x < x_min] <- -1; x_scaled[x > x_max] <- 1 # set over/under passing values to 1/-1
  return(x_scaled)
}
scale_min_max <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  -1 + ((x - x_min) * 2 / (x_max - x_min))
}

reorder_substrings <- function(s, split_symbol = " x ") {
  parts <- strsplit(s, split_symbol, fixed = TRUE)[[1]] # Split the string by split_symbol 
  if (parts[1] > parts[2]) { s <- paste(parts[2], parts[1], sep = split_symbol) # Reorder alphabetically
  } else { s <- paste(parts[1], parts[2], sep = split_symbol) }
  return(s)
}

target_plot = function(data, radii = "all", add_legend = T) {
  library(plotrix)
  colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
  selected_colors = rev(colors)[round((data$score + 1) * 49.5) + 1]; selected_colors[is.na(selected_colors)] = "grey70"
  
  plot(0, 0, type = "n", xlim = c(-10, 10), ylim = c(-10, 10), xlab = "", ylab = "", asp = 1, xaxt='n', yaxt='n', bty='n')  # no axes, no border
  for (i in 1:nrow(data)) { draw.circle(-2, 0, data$radius[i], col=selected_colors[i], border = if (radii == "all") ifelse(i == 1, "black", "grey70") else if (radii == "none") rgb(0,0,0, alpha = 0) else NULL, 
                                        lwd = if (radii == "all") 0.1 else 1.5) } # ifelse(i == 1, 1.5, 0.1)
  if (add_legend) {
    color_bar <- matrix(1:100, ncol=1)
    rasterImage(as.raster(matrix(colors[color_bar], ncol=1)), xleft=7, ybottom=-4, xright=9, ytop=4, interpolate=TRUE)
    tick_marks <- c(1, 0.5, 0, -0.5, -1); tick_positions <- c(4, 2, 0, -2, -4)
    axis(2, at=tick_positions, labels=FALSE, pos= 9, tck=-0.01); axis(4, at=tick_positions, labels=FALSE, pos= 7, tck=-0.01)
    text(9, tick_positions, labels=tick_marks, cex=0.7, pos=4); text(9, 5.2, labels="association", xpd=TRUE, cex = 0.8)
  }
}

spatial_feature_pies = function(data, pie_size = 0.7, fzone_cols = NULL) {
  if (is.null(fzone_cols)) fzone_cols = c(hyp="black", TMEinf = "#F6CF71", Tundiff= "#6eccfa", Tdiff = "grey90")
  data$V3 = -data$V3
  # scale coordinates to have the same range so that plot ratio is not fucked
  data$x = ((data$V4 - min(data$V4)) / (max(data$V4) - min(data$V4))) * 100
  data$y = ((data$V3 - min(data$V3)) / (max(data$V3) - min(data$V3))) * 100
  ggplot() + geom_scatterpie(data = data, aes(x = x, y = y, r = pie_size, alpha = malignancy), color = NA, cols = names(fzone_cols)) +
    scale_alpha_manual(values = c("non_mal" = 0.3, "mal" = 1)) + scale_fill_manual(values = fzone_cols, name = "zone") + theme_void() + 
    coord_fixed()
}


plot_fzone_sig_abundance = function(sig_fzone_abund, zone_cols, ylim = c(1, 1), split.by = NULL, group.by = NULL, normalise = T, add_var = F) {
  if (is.null(split.by)) split.by = rep("all", length(sample_ids))
  plots = lapply(unique(split.by), function(group) {
    by_treshold_mode = is.list(sig_fzone_abund[[2]]) # if second element is another list, assume abundance was calculated by hard treshold (instead of weighted average)
    if (by_treshold_mode) {
      data = lapply(sig_fzone_abund[[1]][split.by == group], function(s) { s[is.na(s)] = 0; s })
      vars = lapply(sig_fzone_abund[[3]][split.by == group], function(s) { s[is.na(s)] = 0; s })
      anno = data.frame(n = Reduce("+", sig_fzone_abund[[2]][split.by == group]), var = colMeans(Reduce("+", vars)))
      anno$sig = rownames(sig_fzone_abund[[2]]); anno$anno = if(add_var) paste0("n=", anno$n, "\nvar=", round(anno$var, digits = 2)) else paste0("n=", anno$n)
    } else {
      data = lapply(sig_fzone_abund[split.by == group], function(s) { s[is.na(s)] = 0; s }) # <- not finished further from here on! Reason: negative 
    }
    data = as.data.frame(Reduce("+", data))
    if (normalise) data = apply(data, 2, function(sig) sig / sum(abs(sig)))
    data[c(3,4,5), ] = -data[c(3,4,5), ]
    data = as.data.frame(data); data$zone = factor(rownames(data), levels = c("hyp", "hyp-adj", "Tdiff", "Tundiff", "angio-adj"))
    data <- data %>% pivot_longer(cols = -zone, names_to = "sig", values_to = "mean_expr")
    if (!is.null(group.by)) {
      data = bind_rows(lapply(names(group.by), function(name) {
        data %>% filter(sig %in% group.by[[name]]) %>% group_by(zone) %>% summarise(mean_expr = mean(mean_expr), .groups = "drop") %>% mutate(sig = name)
      }))
      anno = bind_rows(lapply(names(group.by), function(name) {
        anno %>% filter(sig %in% group.by[[name]]) %>% summarise(anno = paste0("n=", sum(as.integer(sub("n=", "", n))), "\nvar=", sum(var)), .groups = "drop") %>% mutate(sig = name)
      }))
    }
    hyp_sums = sapply(unique(data$sig), function(sig) sum(data[data$zone %in% c("hyp", "hyp-adj") & data$sig == sig, ]$mean_expr))
    data$sig = factor(data$sig, levels = names(hyp_sums[order(hyp_sums, decreasing = T)]))
    ggplot(data, aes(x = sig, y = mean_expr, fill = zone)) + geom_bar(stat = "identity", position = "stack") + geom_hline(yintercept=0, color="black") +
      labs(x = NULL, y = "relative abundance", title = group) + scale_fill_manual(values = zone_cols) + 
      geom_text(data = anno, aes(x = sig, y = hyp_sums, label = anno), inherit.aes = F, color = "black", angle = 90, hjust = -0.5, size = 2) +
      scale_y_continuous(limits = c(-ylim[1], ylim[2]), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1)) + theme_minimal() + #ylim(-ylim[1], ylim[2]) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid.major.x = element_blank())
  })
  if (length(unique(split.by)) > 1) {
    plots[1:(length(plots)-1)] = lapply(plots[1:(length(plots)-1)], function(p) p + theme(legend.position = "none"))
    plots[2:length(plots)] = lapply(plots[2:length(plots)], function(p) p + labs(y = NULL) + theme(axis.text.y = element_blank()))
  }
  return(plots)
}


plot_zone_sig_composition = function(zone_abund, sig_cols, zone_order = NULL, sig_order = NULL, split.by = NULL, normalise = T, split_y_axis = NULL, ylim = c(NA,NA)) {
  if (is.null(split.by)) split.by = rep("all", length(sample_ids))
  if (length(zone_abund) == 2 & is.list(zone_abund[[1]])) { # case for fzone input data, which is based on average zone abund per spot (transform first here)
    zone_abund = lapply(1:length(zone_abund[[1]]), function(i) sweep(zone_abund[[1]][[i]], 2, zone_abund[[2]][, i], `*`))
  }
  plots = lapply(unique(split.by), function(group) {
    data = zone_abund[split.by == group]; data <- lapply(data, function(s) { s[is.na(s)] = 0; s })
    data = as.data.frame(Reduce("+", data))
    if (normalise) data = apply(data, 1, function(zone) zone / sum(abs(zone)))
    if (!is.null(split_y_axis)) data[split_y_axis, ] = -data[split_y_axis, ]
    data = as.data.frame(data); data$sig = factor(rownames(data), levels = if (!is.null(sig_order)) sig_order else unique(rownames(data)))
    data <- data %>% pivot_longer(cols = -sig, names_to = "zone", values_to = "mean_expr")
    data$zone = factor(data$zone, levels = if (!is.null(zone_order)) zone_order else unique(data$zone))
    ggplot(data, aes(x = zone, y = mean_expr, fill = sig)) + geom_bar(stat = "identity", position = "stack") + geom_hline(yintercept=0, color="black") +
      labs(x = NULL, y = "relative abundance", title = group) + scale_fill_manual(values = sig_cols) + 
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), panel.grid.major.x = element_blank()) +
      if (!is.null(split_y_axis)) scale_y_continuous(limits = c(-ylim[1], ylim[2]), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c(1, 0.5, 0, 0.5, 1))
  })
  if (length(unique(split.by)) > 1) {
    plots[1:(length(plots)-1)] = lapply(plots[1:(length(plots)-1)], function(p) p + theme(legend.position = "none"))
    plots[2:length(plots)] = lapply(plots[2:length(plots)], function(p) p + labs(y = NULL) + theme(axis.text.y = element_blank()))
  }
  return(plots)
}

plot_all_sig_zone_abundances = function(szone_abund, fzone_abund, out_dir, add_var = F) {
  # a) abundance of signatures in STRUCTURAL zones 
  plots = plot_mp_zone_abundance(szone_abund, szone_cols, ylim = c(1, 1.1))[[1]] + ggtitle(paste0("all (min_score = ", min_score_threshold, ")"))
  ggsave(paste0("01a_sig_abund_in_szones_", min_score_threshold, ".png"), plot = plots, width = 5, height = 4, path = out_dir)
  plots = plot_mp_zone_abundance(szone_abund, szone_cols, ylim = c(1, 1.1), split.by = subtype)
  ggsave(paste0("01b_sig_abund_in_szones_split-by-subtype_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  plots = plot_mp_zone_abundance(szone_abund, szone_cols, ylim = c(1, 1.1), split.by = zfta_or_not)
  ggsave(paste0("01c_sig_abund_in_szones_split-by-ZFTA-vs-nonZFTA_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 9, height = 4, path = out_dir)
  plots = plot_mp_zone_abundance(szone_abund, szone_cols, ylim = c(1, 1.1), split.by = hyp_class)
  ggsave(paste0("01d_sig_abund_in_szones_split-by-hyp-status_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  
  # b) abundance of signatures in FUNCTIONAL zones 
  plots = plot_fzone_sig_abundance(fzone_abund, fzone_cols, ylim = c(1, 1.1), add_var = add_var)[[1]] + ggtitle(paste0("all (min_score = ", min_score_threshold, ")"))
  ggsave(paste0("02a_sig_abund_in_fzones_", min_score_threshold, ".png"), plot = plots, width = 5, height = 4, path = out_dir)
  plots = plot_fzone_sig_abundance(fzone_abund, fzone_cols, ylim = c(1, 1.1), add_var = add_var, split.by = subtype)
  ggsave(paste0("02b_sig_abund_in_fzones_split-by-subtype_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  plots = plot_fzone_sig_abundance(fzone_abund, fzone_cols, ylim = c(1, 1.1), add_var = add_var, split.by = zfta_or_not)
  ggsave(paste0("02c_sig_abund_in_fzones_split-by-ZFTA-vs-nonZFTA_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 9, height = 4, path = out_dir)
  plots = plot_fzone_sig_abundance(fzone_abund, fzone_cols, ylim = c(1, 1.1), add_var = add_var, split.by = hyp_class)
  ggsave(paste0("02d_sig_abund_in_fzones_split-by-hyp-status_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  
  # c) COMPOSITION of structural and functional zones regarding the signatures
  plots = plot_zone_sig_composition(szone_abund, sig_cols, zone_order = szone_order)[[1]]
  ggsave(paste0("03a_sig_compositions_of_szones_", min_score_threshold, ".png"), plot = plots, width = 7, height = 4, path = out_dir)
  plots = plot_zone_sig_composition(szone_abund, sig_cols, zone_order = szone_order, split.by = subtype)
  ggsave(paste0("03b_sig_compositions_of_szones_split-by-subtype_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  plots = plot_zone_sig_composition(szone_abund, sig_cols, zone_order = szone_order, split.by = zfta_or_not)
  ggsave(paste0("03c_sig_compositions_of_szones_split-by-ZFTA-vs-nonZFTA_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 9, height = 4, path = out_dir)
  plots = plot_zone_sig_composition(szone_abund, sig_cols, zone_order = szone_order, split.by = hyp_class)
  ggsave(paste0("03d_sig_compositions_of_szones_split-by-hyp-status_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  
  plots = plot_zone_sig_composition(fzone_abund, sig_cols, zone_order = fzone_order)[[1]]
  ggsave(paste0("04a_sig_compositions_of_fzones_", min_score_threshold, ".png"), plot = plots, width = 7, height = 4, path = out_dir)
  plots = plot_zone_sig_composition(fzone_abund, sig_cols, zone_order = fzone_order, split.by = subtype)
  ggsave(paste0("04b_sig_compositions_of_fzones_split-by-subtype_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
  plots = plot_zone_sig_composition(fzone_abund, sig_cols, zone_order = fzone_order, split.by = zfta_or_not)
  ggsave(paste0("04c_sig_compositions_of_fzones_split-by-ZFTA-vs-nonZFTA_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 9, height = 4, path = out_dir)
  plots = plot_zone_sig_composition(fzone_abund, sig_cols, zone_order = fzone_order, split.by = hyp_class)
  ggsave(paste0("04d_sig_compositions_of_fzones_split-by-hyp-status_", min_score_threshold, ".png"), plot = wrap_plots(plots, nrow = 1), width = 12, height = 4, path = out_dir)
}


get_sig_scores = function(sigs, plot_sigs = F, out_dir = NULL, nrow = 2) {
  source(file.path(utils_dir, "Single cell/seurat_utils.R"))
  if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, r = T)
  overview_subtitles = c(names(sigs), "[Counts]", "[Features]", "[H&E]")
  all_spots_programs_comp = list(); all_spots_programs_comp_norm = list(); all_spots_programs_comp_spot_norm = list()
  for (s in sample_ids) {
    message(">>> loading sample '", s, "' [", format(Sys.time(), "%d.%m. %X"), "] <<<")
    sstobj = readRDS(paste0("Objects/sstobj_", s, ".rds"))
    ratio = get_spatial_aspect_ratio(sstobj); pt_size = get_spatial_point_size(sstobj, scale_factor = 3.5)
    message("transfering sigs to ", s)
    sstobj = try_add_module_score(sstobj, features = sigs, name = "sig", ctrl = 100, min_ctrl = 50, nbin = 24, min_bin = 18, verbose = T)
    sstobj = transform_modscores_by_cell(sstobj, score_names = paste0("sig", 1:length(sigs)), transformation = "normalisation",
                                         min_score_threshold = 0, transformed_score_names = paste0("sig", 1:length(sigs), "_spot-norm")) # transform module scores across spots
    sstobj = transform_modscores_by_program(sstobj, score_names = paste0("sig", 1:length(sigs)), transformation = "min_max",
                                            transformed_score_names = paste0("sig", 1:length(sigs), "_sig-norm")) # transform module scores across programs
    spots_programs_comp = FetchData(sstobj, vars = paste0("sig", 1:length(sigs))); colnames(spots_programs_comp) = names(sigs)
    spots_programs_comp_spot_norm = FetchData(sstobj, vars = paste0("sig", 1:length(sigs), "_spot-norm")); colnames(spots_programs_comp_spot_norm) = names(sigs)
    spots_programs_comp_norm = FetchData(sstobj, vars = paste0("sig", 1:length(sigs), "_sig-norm")); colnames(spots_programs_comp_norm) = names(sigs)
    all_spots_programs_comp_norm[[s]] = spots_programs_comp_norm; all_spots_programs_comp[[s]] = spots_programs_comp; all_spots_programs_comp_spot_norm[[s]] = spots_programs_comp_spot_norm
    if (plot_sigs) {
      p = c(as.list(SpatialPlot(sstobj, pt.size.factor = pt_size, features = paste0("sig", 1:length(sigs)), image.alpha = 0)), 
            list(SpatialFeaturePlot(sstobj, features = "nCount_Spatial", pt.size.factor = pt_size, image.alpha = 0) + theme(aspect.ratio = ratio), 
                 SpatialFeaturePlot(sstobj, features = "nFeature_Spatial", pt.size.factor = pt_size, image.alpha = 0) + theme(aspect.ratio = ratio), 
                 SpatialFeaturePlot(sstobj, features = NULL, alpha = 0) + NoLegend()))
      p = lapply(1:length(p), function(i) p[[i]] + ggtitle(overview_subtitles[i]) + labs(fill = NULL) + theme(aspect.ratio = ratio, plot.margin = unit(c(0,0,0,0), "mm"), plot.title = element_text(hjust = 0.5)))
      p = (wrap_plots(p[1:length(sigs)], nrow = nrow) | wrap_plots(p[(length(sigs)+1):length(p)], ncol = 1)) + plot_layout(widths = c(12,1))
      ggsave(filename = paste0(s, "_a_overview_raw-scores_counts_features_H&E.png"), plot = p, width = 12, height = 6, path = out_dir)
    
      plots = as.list(SpatialPlot(sstobj, pt.size.factor = pt_size, features = paste0("sig", 1:length(sigs)), image.alpha = 0))
      plots = lapply(1:length(plots), function(i) plots[[i]] + ggtitle(overview_subtitles[i]) + labs(fill = NULL) + theme(aspect.ratio = ratio, plot.margin = unit(c(0,0,0,0), "mm"), plot.title = element_text(hjust = 0.5)))
      global_max <- max(sapply(plots, function(p) max(p$data[, 4], na.rm = TRUE)))
      plots = lapply(plots, function(p) p + scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = c(0, global_max)))
      ggsave(filename = paste0(s, "_b_scores_same-scales.png"), plot = wrap_plots(plots, nrow = nrow), width = 12, height = 6, path = out_dir)
      
      plots = as.list(SpatialPlot(sstobj, pt.size.factor = pt_size, features = paste0("sig", 1:length(sigs), "_spot-norm"), image.alpha = 0))
      plots = lapply(1:length(plots), function(i) plots[[i]] + ggtitle(overview_subtitles[i]) + labs(fill = NULL) + theme(aspect.ratio = ratio, plot.margin = unit(c(0,0,0,0), "mm"), plot.title = element_text(hjust = 0.5)))
      global_max <- max(sapply(plots, function(p) max(p$data[, 4], na.rm = TRUE)))
      plots = lapply(plots, function(p) p + scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = c(0, global_max)))
      ggsave(filename = paste0(s, "_c_scores_same-scales_spot-norm.png"), plot = wrap_plots(plots, nrow = nrow), width = 12, height = 6, path = out_dir)
    }
  }
  return(list(all_spots_programs_comp, all_spots_programs_comp_norm, all_spots_programs_comp_spot_norm))
}

get_sig_szone_abund = function(all_spots_programs_comp, complexity_zones, min_score_threshold) {
  # use raw scores and a hard threshold to consider only spots with significant expression of a signature (emphasises differences much better)
  sig_szone_abund = lapply(sample_ids, function(s) { # simpler variant without significance testing
    sig_scores = all_spots_programs_comp[[s]]
    zone_abund = sapply(colnames(sig_scores), function(sig) {
      tapply(sig_scores[, sig], as.character(complexity_zones[[s]][rownames(sig_scores)]), function(zone_spot_scores) sum(zone_spot_scores > min_score_threshold))
    })
    if (is.null(dim(zone_abund))) { zone_abund = data.frame(t(zone_abund)); rownames(zone_abund) = unique(complexity_zones[[s]][rownames(sig_scores)]); colnames(zone_abund) = colnames(sig_score) }
    for (row_name in c("structured_high", "structured_low", "intermediate", "disorganised_low", "disorganised_high")) { # add missing zones to allow later comparisons!
      if (!row_name %in% rownames(zone_abund)) { zone_abund = rbind(zone_abund, rep(0, ncol(zone_abund))); rownames(zone_abund)[nrow(zone_abund)] = row_name }
    }
    zone_abund = zone_abund[c("structured_high", "structured_low", "intermediate", "disorganised_low", "disorganised_high"), ] # bring in consistent order and remove 'intermediate'
    return(zone_abund)
  }); names(sig_szone_abund) = sample_ids
  return(sig_szone_abund)
}

get_sig_fzone_abund_by_threshold = function(all_spots_programs_comp, functional_zones, min_score_threshold) {
  all_n_spots = list(); all_vars = list()
  sig_fzone_abund = lapply(sample_ids, function(s) {
    # calculate the average zone composition of all spots with significant sig expression
    n_spots = c(); vars = c()
    abund = sapply(colnames(all_spots_programs_comp[[s]]), function(sig) {
      sig_expr = functional_zones[[s]][all_spots_programs_comp[[s]][, sig] > min_score_threshold, ]
      n_spots <<- c(n_spots, nrow(sig_expr))
      vars[[sig]] <<- apply(sig_expr, 2, var, na.rm = T)
      colMeans(sig_expr, na.rm = T)
    })
    names(n_spots) = colnames(all_spots_programs_comp[[s]])
    all_n_spots[[s]] <<- n_spots; all_vars[[s]] <<- data.frame(vars)
    abund
  }); names(sig_fzone_abund) = sample_ids
  return(list(sig_fzone_abund, as.data.frame(all_n_spots), all_vars))
}

get_sig_fzone_abund_by_weighting = function(all_scores, functional_zones) {
  sig_fzone_abund = lapply(sample_ids, function(s) {
    sapply(colnames(all_scores[[s]]), function(sig) {
      zones = functional_zones[[s]]
      sig_expr = all_scores[[s]][, sig]
      colSums(zones * sig_expr, na.rm = T) / sum(abs(sig_expr))
    })
  }); names(sig_fzone_abund) = sample_ids
  return(sig_fzone_abund)
}


reshape_df_pairs_to_mtx = function(df, mp_names) {
  value_matrix <- matrix(NA, nrow = length(mp_names), ncol = length(mp_names), dimnames = list(mp_names, mp_names))
  pvalue_matrix <- matrix(NA, nrow = length(mp_names), ncol = length(mp_names), dimnames = list(mp_names, mp_names))
  for (row in seq_len(nrow(df))) {
    pair <- as.character(df$pair_name[row])
    value <- df$value[row]; pvalue <- df$pvalue[row]
    parts <- unlist(strsplit(pair, " x ")) # Split the pair into its components
    label1 <- parts[1]; label2 <- parts[2]
    idx1 <- match(label1, mp_names); idx2 <- match(label2, mp_names) # Find the indices for the mp_names
    value_matrix[idx1, idx2] <- value; value_matrix[idx2, idx1] <- value
    pvalue_matrix[idx1, idx2] <- pvalue; pvalue_matrix[idx2, idx1] <- pvalue
  }
  value_matrix[is.na(value_matrix)] = 1; pvalue_matrix[is.na(pvalue_matrix)] = 1
  return(list(value_matrix, pvalue_matrix))
}


combine_pvals_dfs = function(pval_dfs) {
  nrow_df = nrow(pval_dfs[[1]]); ncol_df = ncol(pval_dfs[[1]])
  combined_pvals <- matrix(NA, nrow = nrow_df, ncol = ncol_df) # Output matrix
  for (i in 1:nrow_df) {
    for (j in 1:ncol_df) {
      pvals_to_combine <- sapply(pval_dfs, function(df) df[i, j])
      combined_pvals[i, j] = average_indp_pvals_by_fisher(pvals_to_combine) # Extract the p-values at position (i, j) from all dataframes
    }
  }
  combined_pvals_df <- as.data.frame(combined_pvals); colnames(combined_pvals_df) <- colnames(pval_dfs[[1]]); rownames(combined_pvals_df) <- rownames(pval_dfs[[1]])
  # if (adjust) combined_pvals_df = pmin(combined_pvals_df * (nrow_df * ncol_df), 1) # adjust for multiple testing using bonferroni correction (number of tests = n_pairs * n_distances)
  return(combined_pvals_df)
}








                                    


