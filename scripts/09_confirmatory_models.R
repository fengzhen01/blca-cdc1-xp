options(stringsAsFactors = FALSE)

root <- normalizePath(".", winslash = "/")
input <- file.path(root, "results", "GSE176307_scores_and_outcomes.csv")
output <- file.path(root, "results", "GSE176307_regression_results.csv")

d <- read.csv(input, check.names = FALSE)
names(d)[1] <- "sample_id"
d <- d[!duplicated(d$sample_id), ]
d$response <- ifelse(d[["io.response"]] %in% c("CR", "PR"), 1,
                     ifelse(d[["io.response"]] %in% c("SD", "PD"), 0, NA))
d$cross_z <- as.numeric(scale(d$cross_presentation))
d$tmb_z <- as.numeric(scale(as.numeric(d$tmb)))
d$ecog_num <- suppressWarnings(as.numeric(d$ecog))
d$fgfr <- factor(d[["fgfr.mutation.yn"]])
d$progressed_event <- ifelse(d$progressed == "Yes", 1,
                             ifelse(d$progressed == "No", 0, NA))
d$death_event <- ifelse(d$alive == "No", 1,
                        ifelse(d$alive == "Yes", 0, NA))

tidy_glm <- function(model, label) {
  s <- summary(model)$coefficients
  terms <- rownames(s)
  data.frame(
    model = label,
    endpoint = "objective_response_CR_PR",
    term = terms,
    estimate = s[, "Estimate"],
    effect = exp(s[, "Estimate"]),
    conf_low = exp(s[, "Estimate"] - 1.96 * s[, "Std. Error"]),
    conf_high = exp(s[, "Estimate"] + 1.96 * s[, "Std. Error"]),
    p_value = s[, "Pr(>|z|)"],
    n = nobs(model),
    row.names = NULL
  )
}

results <- list()
fit1 <- glm(response ~ cross_z, data = d, family = binomial(), na.action = na.omit)
results[[length(results) + 1]] <- tidy_glm(fit1, "univariable")

# Prespecified, parsimonious adjustment. With only 16 responders, no larger model is justified.
fit2 <- glm(response ~ cross_z + tmb_z, data = d, family = binomial(), na.action = na.omit)
results[[length(results) + 1]] <- tidy_glm(fit2, "adjusted_for_TMB")

if (requireNamespace("survival", quietly = TRUE)) {
  tidy_cox <- function(model, label, endpoint) {
    s <- summary(model)$coefficients
    ci <- summary(model)$conf.int
    data.frame(
      model = label,
      endpoint = endpoint,
      term = rownames(s),
      estimate = s[, "coef"],
      effect = s[, "exp(coef)"],
      conf_low = ci[, "lower .95"],
      conf_high = ci[, "upper .95"],
      p_value = s[, "Pr(>|z|)"],
      n = model$n,
      row.names = NULL
    )
  }
  pfs_fit <- survival::coxph(survival::Surv(pfs, progressed_event) ~ cross_z,
                             data = d, na.action = na.omit)
  os_fit <- survival::coxph(survival::Surv(overall_survival, death_event) ~ cross_z,
                            data = d, na.action = na.omit)
  results[[length(results) + 1]] <- tidy_cox(pfs_fit, "univariable", "PFS")
  results[[length(results) + 1]] <- tidy_cox(os_fit, "univariable", "OS")
}

out <- do.call(rbind, results)
write.csv(out, output, row.names = FALSE, fileEncoding = "UTF-8")
print(out)
