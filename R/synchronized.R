# Synchronized sample normalization -----------------------------------------

# A synchronized analysis preserves the acquired target values.  It only
# scales DNA using an explicitly exported FlowJo G1 population.  The required
# inputs and the G1 containment check deliberately fail closed: no gate is
# inferred from file order, a nearby sample, or a fitted replacement.

synchronized_row_signatures <- function(data, columns) {
  encoded <- lapply(columns, function(column) {
    value <- data[[column]]
    type <- typeof(value)
    rendered <- if (is.double(value)) {
      ifelse(is.na(value), "<NA>", sprintf("%.17g", value))
    } else {
      ifelse(is.na(value), "<NA>", encodeString(as.character(value), quote = '"'))
    }
    paste(type, rendered, sep = ":")
  })
  do.call(paste, c(encoded, sep = "\u001f"))
}

synchronized_identity_columns <- function() {
  "event_index"
}

synchronized_validate_event_index <- function(data, sample_id, population_label) {
  index <- data$event_index
  if (anyNA(index) || any(!nzchar(trimws(as.character(index)))) ||
      anyDuplicated(index)) {
    stop("Invalid, missing, or duplicate event_index values for ", sample_id,
         " / ", population_label, ".", call. = FALSE)
  }
}

synchronized_validate_g1_containment <- function(
    complete, g1, sample_id, minimum_g1_events, population_label = "G1",
    required_columns = synchronized_identity_columns()
) {
  missing_parent <- setdiff(required_columns, names(complete))
  missing_child <- setdiff(required_columns, names(g1))
  if (length(missing_parent) || length(missing_child)) {
    stop("Required event identity/export columns are missing for ", sample_id,
         " / ", population_label, ": ",
         paste(unique(c(missing_parent, missing_child)), collapse = ", "),
         ".", call. = FALSE)
  }
  synchronized_validate_event_index(complete, sample_id, "Single Cells")
  synchronized_validate_event_index(g1, sample_id, population_label)
  if (nrow(g1) < minimum_g1_events) {
    stop("Too few ", population_label, " events for ", sample_id,
         ": found ", nrow(g1), ", require at least ", minimum_g1_events,
         ".", call. = FALSE)
  }
  parent_counts <- table(synchronized_row_signatures(complete, required_columns))
  child_counts <- table(synchronized_row_signatures(g1, required_columns))
  available <- parent_counts[names(child_counts)]
  available[is.na(available)] <- 0L
  if (any(child_counts - available > 0L)) {
    stop(population_label,
         " rows are not an exact multiset subset of Single Cells for ",
         sample_id, ".", call. = FALSE)
  }
  list(
    status = "validated",
    method = "exact_event_index_containment_fixed_identity_column",
    shared_channels = required_columns,
    single_cells_events = nrow(complete),
    child_population = population_label,
    child_events = nrow(g1),
    g1_events = nrow(g1),
    unmatched_or_excess_g1_events = 0L
  )
}

synchronized_normalize_table <- function(events, dna_channel, target_channel,
                                         anchor, dna_2n_value) {
  events$target_raw <- events[[target_channel]]
  events$dna_norm <- events[[dna_channel]] / anchor * dna_2n_value
  events$baseline <- NA_real_
  events$target_norm <- events$target_raw
  events$target_bgsub <- events$target_raw
  events
}
