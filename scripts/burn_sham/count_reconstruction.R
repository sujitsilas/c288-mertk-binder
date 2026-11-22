# The two denominators in obs allow the nested AnnData normalization to be inverted.
h5_vector <- function(file, key) as.vector(rhdf5::h5read(file, key))

h5_category <- function(file, key) {
  prefix <- paste0("obs/", key, "/")
  categories <- h5_vector(file, paste0(prefix, "categories"))
  codes <- h5_vector(file, paste0(prefix, "codes"))
  stopifnot(all(codes >= 0 & codes < length(categories)))
  categories[1L + codes]
}

reconstruct_counts <- function(file, groups, detect_by=NULL) {
  vec <- function(key) h5_vector(file, key)
  genes <- vec("var/_index")
  lib <- vec("obs/nCount_RNA")
  intermediate_lib <- vec("obs/nCount_originalexp")
  features <- vec("obs/nFeature_RNA")
  ptr <- as.double(vec("X/indptr"))
  n <- length(lib)
  stopifnot(length(groups) > 0L, !is.null(names(groups)),
            all(vapply(groups, is.factor, logical(1))),
            all(lengths(groups) == n),
            is.null(detect_by) || detect_by %in% names(groups),
            !anyDuplicated(genes), length(ptr) == n + 1L, all(diff(ptr) == features),
            all(lib == vec("obs/total_counts")), all(lib > 0), all(intermediate_lib > 0),
            rhdf5::h5readAttributes(file, "X")[["encoding-type"]] == "csr_matrix")
  initialize <- function(group) {
    matrix(0, length(genes), nlevels(group), dimnames=list(genes, levels(group)))
  }
  by_group <- lapply(groups, initialize)
  detected_by_group <- if (!is.null(detect_by)) initialize(groups[[detect_by]]) else NULL
  detected <- numeric(length(genes))
  sum_error <- integer(n)
  max_round_error <- 0
  max_log_error <- 0
  for (start in seq(1, n, by=384)) {
    cells <- start:min(start + 383L, n)
    boundaries <- ptr[start:(tail(cells, 1) + 1L)]
    percell <- as.integer(diff(boundaries))
    offset <- boundaries[1L]
    length_chunk <- sum(percell)
    log_values <- as.double(rhdf5::h5read(file, "X/data", start=offset + 1, count=length_chunk))
    index <- as.integer(rhdf5::h5read(file, "X/indices", start=offset + 1, count=length_chunk))
    inner_factor <- rep(intermediate_lib[cells] / 10000, times=percell)
    raw_factor <- rep(lib[cells] / 10000, times=percell)
    inferred <- expm1(expm1(log_values) * inner_factor) * raw_factor
    counts <- round(inferred)
    stopifnot(all(is.finite(counts)), all(counts >= 1),
              all(counts <= rep(lib[cells], times=percell)))
    max_round_error <- max(max_round_error, abs(inferred - counts))
    roundtrip <- log1p(log1p(counts / raw_factor) / inner_factor)
    max_log_error <- max(max_log_error, abs(roundtrip - log_values))
    mat <- Matrix::sparseMatrix(i=index + 1L, p=as.integer(boundaries - offset), x=counts,
                                dims=c(length(genes), length(cells)))
    stopifnot(length(mat@x) == length_chunk)
    sum_error[cells] <- colSums(mat) - lib[cells]
    for (group in names(groups)) {
      labels <- groups[[group]]
      indicator <- Matrix::sparseMatrix(i=seq_along(cells), j=as.integer(labels[cells]), x=1,
                                        dims=c(length(cells), nlevels(labels)))
      by_group[[group]] <- by_group[[group]] + as.matrix(mat %*% indicator)
      if (identical(group, detect_by)) {
        detected_by_group <- detected_by_group + as.matrix((mat > 0) %*% indicator)
      }
    }
    detected <- detected + tabulate(index + 1L, nbins=length(genes))
  }
  qc <- data.frame(cells=n, nonzero_entries=tail(ptr, 1), total_library_counts=sum(lib),
                   cells_matching_library_totals=sum(sum_error == 0),
                   cells_matching_detected_features=sum(diff(ptr) == features),
                   max_cell_total_error=max(abs(sum_error)),
                   max_pre_round_integer_distance=max_round_error,
                   max_roundtrip_logX_error=max_log_error)
  if (any(sum_error != 0) || max_log_error > 0.0001) {
    stop("Count reconstruction did not pass validation: ", paste(qc, collapse=", "))
  }
  stopifnot(all(vapply(by_group, sum, numeric(1)) == sum(lib)))
  list(by_group=by_group, detected=detected, detected_by_group=detected_by_group, qc=qc)
}
