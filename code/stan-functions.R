# functions for fitting models with censored preds in stan ----------------

# modify_stancode
modify_stancode_censored <- function(scode_raw, var_xcens, lower_bound = NULL) {
  
  var_xcens <- make.names(str_remove_all(var_xcens, "[\\._]"))
  n_var_xcens <- length(var_xcens)
  scode <- scode_raw
  
  if (length(lower_bound) == 1 && is.numeric(lower_bound)) lower_bound <- rep(lower_bound, n_var_xcens)
  
  for (i in seq_len(n_var_xcens)) {
    # modify data block:
    n_cens <- glue("int<lower=0> Ncens_{var_xcens[i]};  // number of left-censored")
    j_cens <- glue("array[Ncens_{var_xcens[i]}] int<lower=1> Jcens_{var_xcens[i]};  // positions of left-censored")
    u <- glue("vector[Ncens_{var_xcens[i]}] U_{var_xcens[i]};  // left-censoring limits")
    
    # modifications to parameters block:
    y_cens <- glue("vector<upper=U_{var_xcens[i]}>[Ncens_{var_xcens[i]}] Ycens_{var_xcens[i]};  // estimated left-censored")
    if (!is.null(lower_bound)) {
      if (!is.na(lower_bound[i])) y_cens <- str_replace(y_cens, "^vector<", glue("vector<lower={lower_bound[i]}, "))
    }
    
    # modifications to model block:
    yl <- glue("Yl_{var_xcens[i]}[Jcens_{var_xcens[i]}] = Ycens_{var_xcens[i]}; // add imputed left-censored values")
    
    # modifications to generated quantities block
    gq <- glue("// vector combining observed and missing responses
        vector[N_{var_xcens[i]}] Yl_{var_xcens[i]}_gq = Y_{var_xcens[i]};
        Yl_{var_xcens[i]}_gq[Jmi_{var_xcens[i]}] = Ymi_{var_xcens[i]};
        Yl_{var_xcens[i]}_gq[Jcens_{var_xcens[i]}] = Ycens_{var_xcens[i]}; // add imputed left-censored values
        ")
    
    scode <- scode %>%
      # modifications to data block:
      str_replace(
        "(data \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1", n_cens, j_cens, u), collapse = "\n  ")
      ) %>%
      # modifications to parameters block:
      str_replace(
        "(parameters \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1\n  ", y_cens), collapse = "")
      ) %>%
      # modifications to model block:
      str_replace(
        "(model \\{\n(.|\n)*?)(?=\n    mu_)",
        paste(c("\\1\n    ", yl), collapse = "")
      ) %>%
      str_replace(
        "(generated quantities \\{\n(.|\n)*?)(?=\n\\})",
        paste(c("\\1\n    ", gq), collapse = "")
      )
  }
  
  class(scode) <- "brmsmodel"
  
  return(scode)
  
}


# modify_standata
modify_standata <- function(sdata, model_data, var_xcens, cens_ind) {
  
  varstan <- make.names(str_remove_all(var_xcens, "_"))
  
  for(i in seq_along(var_xcens)) {
    # make logical censoring indicator:
    cens_logical <- model_data[[cens_ind[i]]] == "left"
    
    if (sum(cens_logical) > 0) { # append to standata list only if there are left-censored values
      # number of left-censored
      sdata[[paste0("Ncens_", varstan[i])]] <- sum(cens_logical) 
      # positions of left-censored:
      sdata[[paste0("Jcens_", varstan[i])]] <- as.array(seq_len(nrow(model_data))[cens_logical])
      # left-censoring limits
      # N.B. this assumes that censored values have been replaced by their censoring-limit in the data
      sdata[[paste0("U_", varstan[i])]] <- model_data[[var_xcens[i]]][cens_logical] 
      
    } else message(glue("No left-censored {varstan[i]} values."))
  }
  
  sdata
}

# fit_model --------------------------------------------------------------

fit_stan_model <- function(file,
                           seed,
                           bform,
                           model_data,
                           bpriors = NULL,
                           sample_prior = "no",
                           overwrite = FALSE,
                           var_xcens = NULL,
                           cens_ind = NULL,
                           lower_bound = NULL,
                           knots = NULL,
                           stancode = NULL,
                           ...) {
  
  # load model if it exists
  model_saved <- get_model(file)
  
  # prepare stancode
  code <- brms::make_stancode(
    bform,
    data = model_data,
    prior = bpriors,
    sample_prior = sample_prior
  )
  
  # prepare standata
  data <- brms::make_standata(
    bform,
    data = model_data,
    prior = bpriors,
    sample_prior = sample_prior,
  )
  
  # modify stancode 
  scode <- modify_stancode_censored(
    scode_raw = code,
    var_xcens = var_xcens,
    lower_bound = lower_bound)
  
  # modify standata 
  data <- modify_standata(
    sdata = data,
    model_data = model_data,
    var_xcens = var_xcens, 
    cens_ind = cens_ind)
  
  # fit model:
  
  stanmod <- if (length(model_saved$csvs) > 0 && !overwrite) {
    # read the previously generated cdmdstan files
    brms::read_csv_as_stanfit(model_saved$csvs)
  } else {
    fit_cmdstan_model(
      scode, data, seed, model_saved$path, model_saved$basename, file, ...
    )
  } 
  
  # feed back into brms:
  
  if (length(model_saved$rds) > 0 && !overwrite) {
    brmsmod <- brm(
      bform,
      data = model_data,
      prior = bpriors,
      family = family,
      knots = knots,
      file = file,
      file_refit = "never"
    )
  } else {
    brmsmod <- brm(
      formula = bform,
      data = model_data,
      empty = TRUE
    )
    # save empty fit:
    if (!is.null(file)) brms:::write_brmsfit(brmsmod, file = file)
  }
  
  # add stan model to fit slot:
  brmsmod$fit <- stanmod
  brmsmod <- rename_pars(brmsmod)
  brmsmod$model <- code # replace the original Stan code with the modified code
  
  return(brmsmod)
}


get_model <- function(file) {
  path <- str_remove(file, "\\/[^\\/]+$")# remove base filename
  bname <- str_extract(file, "[^\\/]+$") # extract base filename
  # list csv and rds files matching file path/basename combo:
  csvfiles <- list.files(path = path,
                         pattern = paste0("^", paste0(bname, "[-_]\\d\\.csv")),
                         full.names = TRUE)
  rdsfiles <- list.files(path = path,
                         pattern = paste0("^", paste0(bname, "\\.rds")),
                         full.names = TRUE)
  list(csvs = csvfiles, path = path, basename = bname, rds = rdsfiles)
}

fit_cmdstan_model <- function(code, data, seed, path, basename, file, ...) {
  model_setup <- cmdstanr::cmdstan_model(stan_file = cmdstanr::write_stan_file(code), compile = FALSE)
  model_setup$format(overwrite_file = TRUE, canonicalize = TRUE, backup = FALSE)
  model_setup$compile()
  model <- model_setup$sample(data = data, seed = seed, ...)
  model$save_output_files(
    dir = path,
    basename = basename,
    random = FALSE,
    timestamp = FALSE
  )
  # rstan::read_stan_csv(model$output_files())
  brms::read_csv_as_stanfit(model$output_files())
}


