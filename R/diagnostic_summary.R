#' Print aggregated diagnostic information for a PanelMatch analysis
#'
#' Prints a diagnostic report that is otherwise obtainable via separate
#' calls to \code{summary.PanelMatch()}, \code{compare_treated_observations()},
#' \code{summary.PanelBalance()}, and \code{placebo_test()}. This function
#' does not perform any new estimation, matching, or refinement -- it only
#' aggregates and formats results the user has already computed and
#' supplied. Sections for which the corresponding object is not provided are
#' simply reported as not available; this function will not run those steps
#' on the user's behalf, consistent with the package's philosophy of
#' requiring users to explicitly perform and inspect each stage of an
#' analysis. Each underlying piece already has its own \code{print}/
#' \code{summary}/\code{plot} methods (\code{PanelMatch}, \code{PanelBalance},
#' \code{PanelEstimate}); this function is a higher-level report layered on
#' top, not a new object type.
#'
#' If \code{pm.object} has \code{qoi = "ate"}, the covariate balance and
#' matched-vs-unmatched-treated sections are automatically computed for both
#' \code{"att"} and \code{"atc"} and reported side by side, since both are
#' always available for an \code{ate} \code{PanelMatch} object.
#'
#' In addition to the printed report, this function raises R warnings for
#' four conditions: too few matched sets, both as a proportion (via
#' \code{empty.set.threshold}) and in absolute terms (via
#' \code{min.matched.sets}); poor covariate balance persisting after
#' refinement (via \code{balance.threshold}); and any placebo test estimate
#' whose confidence interval excludes 0 (if \code{placebo.results} is
#' supplied -- any such estimate is flagged, since there is no tunable
#' threshold for this check).
#'
#' @param pm.object A \code{PanelMatch} object. Required.
#' @param panel.data A \code{PanelData} object corresponding to
#'   \code{pm.object}. Required only if \code{covariates} is supplied.
#' @param covariates Character vector of covariate names to pass to
#'   \code{compare_treated_observations()}, comparing matched and unmatched
#'   treated units. If \code{NULL} (default), this section is omitted. If
#'   \code{pm.object} has \code{qoi = "ate"}, this comparison is run for
#'   both \code{"att"} and \code{"atc"} automatically.
#' @param pb.object A \code{PanelBalance} object, as returned by
#'   \code{get_covariate_balance()}. If \code{NULL} (default), the balance
#'   section is omitted. Only refined balance results are shown, regardless
#'   of whether \code{pb.object} was built with \code{include.unrefined =
#'   TRUE}. If \code{pm.object} has \code{qoi = "ate"}, balance results for
#'   both \code{"att"} and \code{"atc"} are reported automatically.
#'
#'   \code{pb.object} may contain balance results for multiple
#'   \code{PanelMatch} configurations, since \code{get_covariate_balance()}
#'   accepts more than one via \code{...}. Only the configuration matching
#'   \code{pm.object} is used; the others (if any) are ignored. This match is
#'   made by variable name: \code{get_covariate_balance()} names each
#'   configuration using the argument name at its own call site (e.g.
#'   \code{get_covariate_balance(pm.obj, ...)} produces a configuration named
#'   \code{"pm.obj"}), so \code{pm.object} must be passed to
#'   \code{diagnostic_summary()} using that same variable name for the match
#'   to succeed -- passing it through an intermediate variable with a
#'   different name, or as a more complex expression, will fail to match and
#'   raise an error naming the configurations that were found instead.
#' @param placebo.results Output of \code{placebo_test(..., plot = FALSE)}.
#'   If \code{NULL} (default), the placebo section is omitted.
#' @param balance.threshold Numeric. Refined covariate-period balance
#'   statistics (in standard deviations) exceeding this value in absolute
#'   terms are counted and flagged (addressing poor post-refinement
#'   balance). Default 0.2.
#' @param empty.set.threshold Numeric, between 0 and 1. If the proportion of
#'   treated observations with no matched controls exceeds this threshold
#'   (for any QOI present in \code{pm.object}), this is flagged (addressing
#'   too few matched sets, as a proportion). Default 0.25.
#' @param min.matched.sets Numeric. If the number of non-empty matched sets
#'   (for any QOI present in \code{pm.object}) falls below this count, this
#'   is flagged (addressing too few matched sets, in absolute terms).
#'   Default 10.
#' @param digits Number of digits to round printed numeric output to.
#'   Default 3.
#' @param ... Additional arguments passed to \code{summary.PanelMatch()}.
#'
#' @return Invisibly returns a plain list with the following components,
#'   for programmatic use (e.g. exporting a table for a manuscript
#'   appendix):
#'   \describe{
#'     \item{checks}{A \code{data.frame} with one row per threshold check,
#'       giving the QOI, configuration (if applicable), metric, a
#'       human-readable label, the observed value, the threshold, which
#'       direction is considered problematic, and whether the check is
#'       flagged.}
#'     \item{matched.set.summary}{As returned by \code{summary.PanelMatch()}.}
#'     \item{matched.treated.summary}{As returned by
#'       \code{compare_treated_observations()}, if \code{covariates} was
#'       supplied.}
#'     \item{balance.summary}{As returned by \code{summary.PanelBalance()},
#'       if \code{pb.object} was supplied.}
#'     \item{placebo.table}{The placebo test table, if \code{placebo.results}
#'       was supplied.}
#'   }
#'
#' @export
diagnostic_summary <- function(pm.object,
                               panel.data = NULL,
                               covariates = NULL,
                               pb.object = NULL,
                               placebo.results = NULL,
                               balance.threshold = 0.2,
                               empty.set.threshold = 0.25,
                               min.matched.sets = 10,
                               digits = 3,
                               ...) {
  
  # captured before pm.object is otherwise touched -- this is how pb.object
  # (which may contain balance results for several PanelMatch configurations,
  # since get_covariate_balance() accepts multiple objects via ...) gets
  # matched back to this specific pm.object. get_covariate_balance() names
  # each config using the deparsed argument name at ITS call site (e.g.
  # get_covariate_balance(pm.obj, ...) produces a config named "pm.obj"), so
  # pm.object must be passed to diagnostic_summary() using that same variable
  # name for the match to succeed.
  pm.name <- deparse(substitute(pm.object))
  
  if (!inherits(pm.object, "PanelMatch")) {
    stop("pm.object must be a PanelMatch object.")
  }
  if (!is.null(covariates) && is.null(panel.data)) {
    stop("panel.data must be supplied when covariates is specified.")
  }
  if (!is.null(pb.object) && !inherits(pb.object, "PanelBalance")) {
    stop("pb.object must be a PanelBalance object.")
  }
  if (!is.numeric(empty.set.threshold) || empty.set.threshold < 0 || empty.set.threshold > 1) {
    stop("empty.set.threshold must be numeric, between 0 and 1.")
  }
  if (!is.numeric(min.matched.sets) || min.matched.sets < 0) {
    stop("min.matched.sets must be a non-negative number.")
  }
  
  qoi.in <- attr(pm.object, "qoi")
  is.ate <- identical(qoi.in, "ate")
  
  # ============================================================
  # 1. Compute all pieces (no new estimation/matching -- these
  #    functions just reformat what the user already computed)
  # ============================================================
  
  # --- matched set sizes: reuse summary.PanelMatch() as-is ---
  matched.set.summary <- summary(pm.object, ...)
  
  # --- empty matched set proportions and absolute counts ---
  # a low proportion of empty sets can still leave very few usable matched
  # sets in absolute terms (e.g. 3 of 5 treated units matched looks fine as
  # a proportion, but 3 matched sets is thin) -- hence both checks.
  empty.set.props <- lapply(matched.set.summary, function(df) {
    empty.row   <- df[df$quantity == "Number of empty matched sets", "value"]
    treated.row <- df[df$quantity == "Number of treated units", "value"]
    if (length(empty.row) == 1 && length(treated.row) == 1 && treated.row > 0) {
      empty.row / treated.row
    } else {
      NA_real_
    }
  })
  n.matched.sets <- lapply(matched.set.summary, function(df) {
    empty.row   <- df[df$quantity == "Number of empty matched sets", "value"]
    treated.row <- df[df$quantity == "Number of treated units", "value"]
    if (length(empty.row) == 1 && length(treated.row) == 1) {
      treated.row - empty.row
    } else {
      NA_real_
    }
  })
  
  # --- matched vs. unmatched treated comparison: reuse compare_treated_observations() as-is ---
  # compare_treated_observations() reads attr(pm.object, "qoi") to select
  # which matched sets to use (pm.object[[qoi]]). For qoi = "ate" objects
  # there is no pm.object[["ate"]] element (only "att"/"atc"), so for that
  # case we run it once per QOI, each time pointing at a version of
  # pm.object with its qoi attribute overridden -- this doesn't compute
  # anything new, it just directs compare_treated_observations() at the
  # correct existing subset.
  matched.treated.summary <- NULL
  if (!is.null(covariates)) {
    if (is.ate) {
      matched.treated.summary <- lapply(c("att", "atc"), function(q) {
        pm.sub <- pm.object
        attr(pm.sub, "qoi") <- q
        compare_treated_observations(pm.sub, panel.data, covariates)
      })
      names(matched.treated.summary) <- c("att", "atc")
    } else {
      matched.treated.summary <- compare_treated_observations(pm.object, panel.data, covariates)
    }
  }
  
  # --- covariate balance: reuse summary.PanelBalance() as-is, then just count/flag ---
  # only refined balance is considered -- summary.PanelBalance(include.unrefined = FALSE)
  # never returns the "_unrefined" columns, so there's nothing to filter out here.
  flag.balance <- function(mat) {
    list(
      n.exceed = sum(abs(mat) > balance.threshold, na.rm = TRUE),
      n.total  = sum(!is.na(mat))
    )
  }
  
  balance.summary <- NULL
  balance.flags <- NULL
  if (!is.null(pb.object)) {
    # pb.object may contain balance results for several PanelMatch
    # configurations (get_covariate_balance(pm1, pm2, ...)) -- only show the
    # one that corresponds to pm.object. Subsetting by integer position
    # (rather than by name) matters here: [.PanelBalance subsets its parallel
    # unrefined.balance.results attribute using the same index, and that
    # attribute's names carry a "_unrefined" suffix, so a name-based index
    # would fail to line the two up correctly.
    pos <- match(pm.name, names(pb.object))
    if (is.na(pos)) {
      stop(sprintf(
        paste0(
          "pb.object does not contain a configuration named '%s'. diagnostic_summary() ",
          "matches entries in pb.object to pm.object by variable name, which must match ",
          "the name used when pm.object was passed to get_covariate_balance() (e.g. ",
          "get_covariate_balance(%s, panel.data = ..., covariates = ...)). Configurations ",
          "found in pb.object: %s."
        ),
        pm.name, pm.name, paste(names(pb.object), collapse = ", ")
      ))
    }
    pb.sub <- pb.object[pos]
    
    if (is.ate) {
      balance.summary <- lapply(c("att", "atc"), function(q) summary(pb.sub, qoi = q, include.unrefined = FALSE))
      names(balance.summary) <- c("att", "atc")
      balance.flags <- lapply(balance.summary, function(qoi.list) lapply(qoi.list, flag.balance))
    } else {
      balance.summary <- summary(pb.sub, include.unrefined = FALSE)
      balance.flags <- lapply(balance.summary, flag.balance)
    }
  }
  
  # --- placebo test: pull pass/fail directly from the CI already computed by placebo_test() ---
  placebo.summary <- NULL
  if (!is.null(placebo.results)) {
    est <- placebo.results$estimates
    ci  <- placebo.results$conf.intervals
    if (is.null(est) || is.null(ci)) {
      stop("placebo.results must contain 'estimates' and 'conf.intervals'. Update placebo_test() to retain conf.intervals, then re-run placebo_test(..., plot = FALSE).")
    }
    
    # a placebo estimate "fails" if its confidence interval excludes 0
    significant <- !(ci[, 1] <= 0 & ci[, 2] >= 0)
    
    placebo.table <- data.frame(
      period      = rownames(ci),
      estimate    = as.numeric(est),
      ci.lower    = ci[, 1],
      ci.upper    = ci[, 2],
      significant = ifelse(significant, "Yes", "No"),
      row.names   = NULL
    )
    
    placebo.summary <- list(
      table            = placebo.table,
      ci.colnames      = colnames(ci),
      n.significant    = sum(significant, na.rm = TRUE),
      n.total          = sum(!is.na(significant)),
      prop.significant = mean(significant, na.rm = TRUE)
    )
  }
  
  # ============================================================
  # 2. Build the checks table (single source of truth for both
  #    the warnings below and the printed NOTE lines)
  # ============================================================
  
  rows <- list()
  
  for (nm in names(matched.set.summary)) {
    prop   <- empty.set.props[[nm]]
    n.sets <- n.matched.sets[[nm]]
    
    rows[[length(rows) + 1]] <- data.frame(
      qoi = nm, config = NA_character_, metric = "empty_set_proportion",
      label = sprintf("%s: %.1f%% of treated units unmatched (threshold: %.1f%%)",
                      toupper(nm), 100 * prop, 100 * empty.set.threshold),
      value = prop, threshold = empty.set.threshold, direction = "above",
      flagged = isTRUE(!is.na(prop) && prop > empty.set.threshold),
      stringsAsFactors = FALSE
    )
    
    rows[[length(rows) + 1]] <- data.frame(
      qoi = nm, config = NA_character_, metric = "n_matched_sets",
      label = sprintf("%s: %d non-empty matched sets (threshold: %d)",
                      toupper(nm), as.integer(n.sets), as.integer(min.matched.sets)),
      value = n.sets, threshold = min.matched.sets, direction = "below",
      flagged = isTRUE(!is.na(n.sets) && n.sets < min.matched.sets),
      stringsAsFactors = FALSE
    )
  }
  
  add_balance_rows <- function(bal.flags, qoi.label) {
    out <- list()
    for (cfg in names(bal.flags)) {
      fl <- bal.flags[[cfg]]
      if (is.null(fl)) next
      out[[length(out) + 1]] <- data.frame(
        qoi = qoi.label, config = cfg, metric = "balance",
        label = sprintf("%s: %d of %d balance stats exceed |%.2f| SD (%s)",
                        toupper(qoi.label), fl$n.exceed, fl$n.total, balance.threshold, cfg),
        value = fl$n.exceed, threshold = 0, direction = "above",
        flagged = fl$n.exceed > 0,
        stringsAsFactors = FALSE
      )
    }
    out
  }
  if (!is.null(balance.flags)) {
    if (is.ate) {
      for (q in names(balance.flags)) {
        rows <- c(rows, add_balance_rows(balance.flags[[q]], q))
      }
    } else {
      rows <- c(rows, add_balance_rows(balance.flags, qoi.in))
    }
  }
  
  if (!is.null(placebo.summary)) {
    ps <- placebo.summary
    rows[[length(rows) + 1]] <- data.frame(
      qoi = NA_character_, config = NA_character_, metric = "placebo",
      label = sprintf("Placebo: %d of %d estimates have CIs excluding 0",
                      ps$n.significant, ps$n.total),
      value = ps$n.significant, threshold = 0, direction = "above",
      flagged = ps$n.significant > 0,
      stringsAsFactors = FALSE
    )
  }
  
  checks <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(
    qoi = character(0), config = character(0), metric = character(0),
    label = character(0), value = numeric(0), threshold = numeric(0),
    direction = character(0), flagged = logical(0), stringsAsFactors = FALSE
  )
  rownames(checks) <- NULL
  
  # ============================================================
  # 3. Warnings (reviewer comments 2 and 3): derived directly
  #    from the checks table above, so this can't drift out of
  #    sync with what gets printed
  # ============================================================
  
  warn_group <- function(metric, prefix) {
    sub <- checks[checks$metric == metric & checks$flagged, , drop = FALSE]
    if (nrow(sub) > 0) {
      warning(paste0(prefix, ": ", paste(sub$label, collapse = "; ")), call. = FALSE)
    }
  }
  warn_group("empty_set_proportion", "Empty matched set check flagged")
  warn_group("n_matched_sets",       "Minimum matched sets check flagged")
  warn_group("balance",              "Covariate balance check flagged")
  warn_group("placebo",              "Placebo test check flagged")
  
  # ============================================================
  # 4. Print the report
  # ============================================================
  
  section.rule <- strrep("=", 60)
  sub.rule     <- strrep("-", 60)
  
  check.note <- function(metric.name, qoi = NA_character_, config = NA_character_) {
    sub <- checks[checks$metric == metric.name & checks$flagged, , drop = FALSE]
    if (!is.na(qoi))    sub <- sub[sub$qoi == qoi, , drop = FALSE]
    if (!is.na(config)) sub <- sub[sub$config == config, , drop = FALSE]
    for (lbl in sub$label) cat(sprintf("NOTE: %s\n", lbl))
  }
  
  print_balance_block <- function(bal.list, qoi.label) {
    for (nm in names(bal.list)) {
      cat(sprintf("\nConfiguration: %s\n", nm))
      print(round(bal.list[[nm]], digits))
      check.note("balance", qoi = qoi.label, config = nm)
    }
  }
  
  # formats the plain list returned by compare_treated_observations() --
  # that function returns unclassed data, so this formatting lives here
  # rather than as a dispatched print method or standalone function.
  print_matched_treated <- function(x) {
    cat("Treated unit-times with matched controls:              ", x$n_has_match, "\n")
    cat("Treated unit-times with no matched controls:           ", x$n_no_match,
        sprintf(" (%s%% of treated observations)", round(x$pct_no_match, digits)), "\n")
    cat("Unmatched treated unit-times with no viable comparison:", x$n_no_viable,
        sprintf(" (%s%% of unmatched treated observations)", round(x$pct_no_viable, digits)), "\n")
    
    if (is.null(x$covariate_diffs)) {
      cat("\nNo covariate comparison available (no cohort had both matched and unmatched treated units).\n")
      return(invisible(NULL))
    }
    
    cat("\nAverage covariate differences (matched - unmatched), by treatment cohort:\n\n")
    
    print_df                    <- x$covariate_diffs
    print_df$mean_diff          <- round(print_df$mean_diff, digits)
    print_df$weighted_mean_diff <- round(print_df$weighted_mean_diff, digits)
    print(print_df, row.names = FALSE)
  }
  
  cat(section.rule, "\n")
  cat("PanelMatch Diagnostic Summary\n")
  cat(section.rule, "\n")
  
  section.num <- 0
  
  # --- 1. Matched Set Sizes -------------------------------------------------
  section.num <- section.num + 1
  cat("\n", sub.rule, "\n", sep = "")
  cat(sprintf("[%d] MATCHED SET SIZES\n", section.num))
  cat(sub.rule, "\n")
  for (nm in names(matched.set.summary)) {
    cat(sprintf("\nQOI: %s\n", toupper(nm)))
    print(matched.set.summary[[nm]], row.names = FALSE)
    if (is.null(matched.treated.summary)) {
      # only show these notes here if section [2] below (which reports the
      # same figures in more detail) isn't present
      check.note("empty_set_proportion", qoi = nm)
      check.note("n_matched_sets", qoi = nm)
    }
  }
  
  # --- 2. Matched vs. Unmatched Treated Units --------------------------------
  if (!is.null(matched.treated.summary)) {
    section.num <- section.num + 1
    cat("\n", sub.rule, "\n", sep = "")
    cat(sprintf("[%d] MATCHED VS. UNMATCHED TREATED UNITS\n", section.num))
    cat(sub.rule, "\n")
    
    if (is.ate) {
      for (q in names(matched.treated.summary)) {
        cat(sprintf("\nQOI: %s\n", toupper(q)))
        print_matched_treated(matched.treated.summary[[q]])
        check.note("empty_set_proportion", qoi = q)
        check.note("n_matched_sets", qoi = q)
      }
    } else {
      cat("\n")
      print_matched_treated(matched.treated.summary)
      check.note("empty_set_proportion", qoi = qoi.in)
      check.note("n_matched_sets", qoi = qoi.in)
    }
  }
  
  # --- 3. Covariate Balance --------------------------------------------------
  section.num <- section.num + 1
  cat("\n", sub.rule, "\n", sep = "")
  cat(sprintf("[%d] COVARIATE BALANCE\n", section.num))
  cat(sub.rule, "\n")
  if (!is.null(balance.summary)) {
    if (is.ate) {
      for (q in names(balance.summary)) {
        cat(sprintf("\nQOI: %s\n", toupper(q)))
        print_balance_block(balance.summary[[q]], q)
      }
    } else {
      print_balance_block(balance.summary, qoi.in)
    }
  } else {
    cat("\nNot provided. Pass a PanelBalance object via pb.object to include this section.\n")
  }
  
  # --- 4. Placebo Test --------------------------------------------------------
  section.num <- section.num + 1
  cat("\n", sub.rule, "\n", sep = "")
  cat(sprintf("[%d] PLACEBO TEST\n", section.num))
  cat(sub.rule, "\n")
  if (!is.null(placebo.summary)) {
    tbl <- placebo.summary$table
    tbl$estimate <- round(tbl$estimate, digits)
    tbl$ci.lower <- round(tbl$ci.lower, digits)
    tbl$ci.upper <- round(tbl$ci.upper, digits)
    names(tbl)[names(tbl) == "ci.lower"] <- placebo.summary$ci.colnames[1]
    names(tbl)[names(tbl) == "ci.upper"] <- placebo.summary$ci.colnames[2]
    cat("\n")
    print(tbl, row.names = FALSE)
    cat("\n")
    check.note("placebo")
  } else {
    cat("\nNot provided. Pass the output of placebo_test(..., plot = FALSE) via placebo.results to include this section.\n")
  }
  
  cat("\n", section.rule, "\n", sep = "")
  
  # ============================================================
  # 5. Return the underlying pieces invisibly, as a plain list
  #    (no new class -- each piece already has its own methods
  #    via its original object type)
  # ============================================================
  
  invisible(list(
    checks                   = checks,
    matched.set.summary      = matched.set.summary,
    matched.treated.summary  = matched.treated.summary,
    balance.summary          = balance.summary,
    placebo.table            = if (!is.null(placebo.summary)) placebo.summary$table else NULL
  ))
}
