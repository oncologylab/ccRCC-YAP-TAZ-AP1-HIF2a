# Figure 3A: treatment effects on log2 tumor-volume fold changes.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript figures/figure_3/tumor_growth.R data/figure3a_tumor_growth.csv output_directory")
}
if (!requireNamespace("nlme", quietly=TRUE)) stop("Install the R package nlme")
if (file.exists(args[2])) stop("Choose a new output directory")
dat <- read.csv(args[1], stringsAsFactors=FALSE)
required <- c("group", "label", "mouse", "day", "log2_fold_change")
if (!all(required %in% names(dat))) stop("Required columns: ", paste(required, collapse=", "))
if (anyNA(dat[required])) stop("Missing input values")
for (column in c("group", "day", "log2_fold_change")) {
    if (!is.numeric(dat[[column]]) || any(!is.finite(dat[[column]]))) {
        stop("Finite numeric values required: ", column)
    }
}
if (!setequal(unique(dat$group), 0:5)) stop("Groups 0 through 5 are required")
if (any(dat$mouse == "") || any(dat$label == "")) stop("Missing mouse or group identifiers")
if (anyDuplicated(dat[c("mouse", "day")])) stop("Duplicate mouse/day observations")
if (any(vapply(split(dat$group, dat$mouse), function(x) length(unique(x)), integer(1)) != 1)) {
    stop("Each mouse must belong to one group")
}
if (any(vapply(split(dat$label, dat$group), function(x) length(unique(x)), integer(1)) != 1)) {
    stop("Each group must have one label")
}
if (any(dat$day < 0 | dat$day > 51)) stop("Figure 3A uses days 0 through 51")

# Source group order: vehicle, TEAD-i, HIF2a-i, combination, shYT+vehicle, shYT+HIF2a-i.
pairs <- list(c(0, 1), c(0, 4), c(2, 3), c(2, 5))
fits <- list()
results <- list()
for (pair in pairs) {
    name <- paste(pair, collapse="_")
    repeated_anova <- pair[1] == 0
    last_day <- if (repeated_anova) 37 else 51
    d <- dat[dat$group %in% pair & dat$day <= last_day, ]
    group_a <- unique(d$label[d$group == pair[1]])
    group_b <- unique(d$label[d$group == pair[2]])
    d$group <- factor(d$group, levels=pair)
    d$time <- factor(d$day)
    d$mouse <- factor(d$mouse)
    contrasts(d$group) <- contr.sum(2)
    contrasts(d$time) <- contr.sum(nlevels(d$time))
    if (repeated_anova) {
        if (!all(table(d$mouse, d$time) == 1)) stop("RM ANOVA requires complete repeated measurements")
        model <- stats::aov(log2_fold_change ~ group * time + Error(mouse / time), data=d)
        treatment <- summary(model)[["Error: mouse"]][[1]]
        f_value <- treatment["group", "F value"]
        df_num <- treatment["group", "Df"]
        df_den <- treatment["Residuals", "Df"]
        p_value <- treatment["group", "Pr(>F)"]
        method <- "Two-way repeated-measures ANOVA"
    } else {
        # Include every available value, including day 37; no imputation or exclusions.
        model <- nlme::lme(log2_fold_change ~ group + time, random=~1|mouse,
                           data=d, method="REML", na.action=na.fail,
                           control=nlme::lmeControl(msMaxIter=200, tolerance=1e-9, msTol=1e-9))
        treatment <- summary(model)$tTable["group1", ]
        f_value <- unname(treatment["t-value"]^2)
        df_num <- 1
        df_den <- unname(treatment["DF"])
        p_value <- unname(treatment["p-value"])
        method <- "Additive REML mixed-effects model"
    }
    annotation <- if (p_value < 0.0001) "****" else if (p_value < 0.001) "***" else
                  if (p_value < 0.01) "**" else if (p_value < 0.05) "*" else "ns"
    results[[name]] <- data.frame(comparison=name, group_a=group_a, group_b=group_b,
        method=method, effect="Treatment main effect", day_min=min(d$day), day_max=max(d$day),
        n_observations=nrow(d), n_mice=nlevels(d$mouse), df_num=df_num, df_den=df_den,
        F=f_value, p_value=p_value, adjustment="None (separate pairwise models)", annotation=annotation)
    fits[[name]] <- model
}
result <- do.call(rbind, results)
if (!dir.create(args[2], recursive=TRUE)) stop("Could not create output directory")
write.table(result, file.path(args[2], "treatment_effects.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(fits, file.path(args[2], "models.rds"))
capture.output(sessionInfo(), file=file.path(args[2], "sessionInfo.txt"))
print(result, row.names=FALSE)
