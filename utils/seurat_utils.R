try_add_module_score <- function(sobj, features, name, ctrl = 100, ctrl_decrement = 5, min_ctrl = 50, nbin = 24, min_bin = 18, verbose = F) {
  current_ctrl <- ctrl  # Startwert für ctrl
  repeat {
    if (current_ctrl < min_ctrl) break
    tryCatch({
      sobj <- Seurat::AddModuleScore(sobj, slot = "data", features = features, name = name, ctrl = current_ctrl)
      break # wenn erfolgreich, beende die Schleife
    }, error = function(e) { current_ctrl <<- current_ctrl - ctrl_decrement; if (verbose) message(e) }) # try with reduced number of ctrl features
  }
  if (current_ctrl >= min_ctrl) {
    if (current_ctrl != ctrl) message("[WARNING] Highest working number of ctrl features lower than desired: ", current_ctrl, " instead of ", ctrl)
    return(sobj)
  } else {
    message("[WARNING] No working number of ctrl features found for min_ctrl = ", min_ctrl, "!")
    current_bin = nbin
    repeat {
      if (current_bin < min_bin) break
      tryCatch({
        sobj <- Seurat::AddModuleScore(sobj, slot = "data", features = features, name = name, ctrl = ctrl, nbin = current_bin)
        break # wenn erfolgreich, beende die Schleife
      }, error = function(e) { current_bin <<- current_bin - 1; if(verbose) message(e) }) # try with reduced number of bins
    }
    if (current_bin >= min_bin) { message("[WARNING] Highest working number of bins lower than desired: ", current_bin, " instead of ", nbin, " (low variance in data?)")
    } else { message("[ERROR] No working number of bins found for min_bin = ", min_bin, "!") }
    return(sobj)
  }
}

add_annotation_from_excelsheet = function(sobj, metadata, col_name, annotation_id, sample_id_col_name = "Study ID") {
  subtype_annotation = data.frame(metadata[[col_name]], "s_ids" = metadata[[sample_id_col_name]])
  ljoined_anno = merge(subtype_annotation, data.frame("s_ids" = sobj@meta.data$orig.ident), by = "s_ids", all.x = T, sort = F)
  sobj@meta.data[[annotation_id]] = ljoined_anno[[2]]
  return(sobj)
}
