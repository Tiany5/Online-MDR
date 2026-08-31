project <- Sys.getenv("OMDRC_ROOT", unset = getwd())
setwd(project)

r_files <- list.files(".", pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
syntax_errors <- character()

for (file in r_files) {
  parsed <- tryCatch(parse(file), error = identity)
  if (inherits(parsed, "error")) {
    syntax_errors[file] <- conditionMessage(parsed)
  }
}

required <- sort(c(
  "data.table", "doParallel", "dplyr", "foreach", "ggplot2", "ggpubr",
  "glmnet", "grid", "isotree", "iterators", "kedd", "Matrix", "mgcv",
  "onlineFDR", "parallel", "patchwork", "RColorBrewer", "REBayes",
  "rstudioapi", "scales", "stats", "tidyr", "utils"
))

status <- vapply(required, requireNamespace, logical(1L), quietly = TRUE)
versions <- vapply(
  required,
  function(pkg) if (status[[pkg]]) as.character(packageVersion(pkg)) else NA_character_,
  character(1L)
)

cat("R_VERSION\t", R.version.string, "\n", sep = "")
cat("R_FILES\t", length(r_files), "\n", sep = "")
cat("SYNTAX_ERRORS\t", length(syntax_errors), "\n", sep = "")
if (length(syntax_errors)) {
  for (file in names(syntax_errors)) cat("SYNTAX_ERROR\t", file, "\t", syntax_errors[[file]], "\n", sep = "")
}
cat("PACKAGES\t", length(required), "\n", sep = "")
for (pkg in required) {
  cat("PACKAGE\t", pkg, "\t", if (status[[pkg]]) versions[[pkg]] else "MISSING", "\n", sep = "")
}

smoke <- list(
  onlineFDR = function() {
    result <- onlineFDR::SAFFRON(c(0.001, 0.2, 0.03, 0.8, 0.01), alpha = 0.1)
    stopifnot(is.data.frame(result), "R" %in% names(result), nrow(result) == 5L)
  },
  isotree = function() {
    set.seed(42)
    x <- matrix(rnorm(200), ncol = 2L)
    model <- isotree::isolation.forest(x, ntrees = 10L, seed = 42L)
    scores <- predict(model, x)
    stopifnot(length(scores) == nrow(x), all(is.finite(scores)))
  },
  glmnet = function() {
    set.seed(42)
    x <- matrix(rnorm(200), ncol = 4L)
    y <- rbinom(nrow(x), 1L, 0.4)
    model <- glmnet::glmnet(x, y, family = "binomial")
    stopifnot(inherits(model, "glmnet"))
  },
  mgcv = function() {
    set.seed(42)
    dat <- data.frame(x = runif(100), y = rbinom(100, 1L, 0.5))
    model <- mgcv::gam(y ~ s(x, k = 5L), data = dat, family = binomial())
    stopifnot(inherits(model, "gam"))
  },
  project_omdrc = function() {
    env1 <- new.env(parent = globalenv())
    env2 <- new.env(parent = globalenv())
    sys.source("code-semi/OMDRC.R", envir = env1)
    sys.source("application/OMDRC.R", envir = env2)
    stopifnot(exists("OMDRC_DD", envir = env1), exists("OMDRC_DD", envir = env2))
  }
)

smoke_ok <- logical(length(smoke))
names(smoke_ok) <- names(smoke)
for (name in names(smoke)) {
  result <- tryCatch({
    smoke[[name]]()
    NULL
  }, error = identity)
  smoke_ok[[name]] <- is.null(result)
  cat(
    "SMOKE\t", name, "\t",
    if (is.null(result)) "OK" else paste("FAILED", conditionMessage(result)),
    "\n", sep = ""
  )
}

if (length(syntax_errors) || any(!status) || any(!smoke_ok)) quit(status = 1L)
