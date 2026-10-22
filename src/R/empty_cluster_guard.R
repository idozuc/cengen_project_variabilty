# Source after src/R/clustering.R: it wraps estimate_resolution_stability().
#
# Guard for a pipeline edge case: a G-component fit whose hard labels use fewer
# than G clusters (an empty component) made estimate_resolution_stability()
# stop, aborting the whole cell type. Here such a resolution instead gets the
# same "no stability" result the pipeline assigns to structurally unsupported
# resolutions, so it is rejected as unstable and the other resolutions proceed.
# For fits without an empty component this is a no-op.
.pipeline_estimate_resolution_stability <- estimate_resolution_stability

estimate_resolution_stability <- function(scores, full_labels, g, n_repeats = 50L, ...) {
  if (nlevels(factor(full_labels)) != g) {
    return(list(
      ari = rep(NA_real_, n_repeats),
      median_ari = NA_real_,
      successful_repeats = 0L,
      failed_repeats = n_repeats
    ))
  }
  .pipeline_estimate_resolution_stability(scores, full_labels, g, n_repeats, ...)
}
