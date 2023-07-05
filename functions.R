# Functions for things

### Annotate methylation sites
# Use genome_sites file, which has all genome annotations together for K-12 and
#   columns 'Type' (category of feature like promoter), 'Site' (name of feature),
#   'Left' and 'Right' (start and end of feature)
# Input a file like AvWT which has a column 'start' (position of methylated site)
# Other metadata columns do not matter
annotateMethylSites <- function(methyl_df, meta_df, location) {
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[meta_df$Left <= position &
                     meta_df$Right >= position,]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'No_Feature'] <- '1'
      next
    }
    sites_at_position <- meta_df[meta_df$Left <= position &
                                   meta_df$Right >= position,]
    for (i in 1:nrow(sites_at_position)) {
      methyl_df[methyl_df[[location]] == position, toString(sites_at_position[i,'Type'])] <- sites_at_position[i,'Site']
    }
  }
  return(methyl_df)
}

### Annotate methylation sites relative to TSS
# Similar to the above, but only using Transcription-Units
# Use genome_sites file, which has all genome annotations together for K-12 and
#   columns 'Type' (category of feature like promoter), 'Site' (name of feature),
#   'Left' and 'Right' (start and end of feature)
# Input a file like AvWT which has a column 'start' (position of methylated site)
# Other metadata columns do not matter
annotateTSS <- function(methyl_df, meta_df, location, size) {
  meta_df <- meta_df %>% 
    filter(Type == 'Transcription-Units')
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[(meta_df$Strand == '+' &
                      meta_df$Left - size <= position &
                      meta_df$Left + size >= position) |
                     (meta_df$Strand == '-' &
                      meta_df$Right - size <= position &
                      meta_df$Right + size >= position), ]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'NoTSS'] <- 'X'
      next
    }
    SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
                                     meta_df$Left - size <= position &
                                     meta_df$Left + size >= position, ]
    AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
                                         meta_df$Right - size <= position &
                                         meta_df$Right + size >= position, ]
    for (i in 1:nrow(SenseTU_at_position)) {
      if (nrow(SenseTU_at_position) == 0) {
        next
      }
      methyl_df[methyl_df[[location]] == position, paste0('RelPos_+', i)] <-
        position - SenseTU_at_position[i, 'Left']
    }
    for (i in 1:nrow(AntisenseTU_at_position)) {
      if (nrow(AntisenseTU_at_position) == 0) {
        next
      }
      methyl_df[methyl_df[[location]] == position, paste0('RelPos_-', i)] <-
        AntisenseTU_at_position[i, 'Right'] - position
    }
  }
  return(methyl_df)
}

### Annotate methylation sites relative to TTS
# Similar to the above, but only using Transcription-Units
# Use genome_sites file, which has all genome annotations together for K-12 and
#   columns 'Type' (category of feature like promoter), 'Site' (name of feature),
#   'Left' and 'Right' (start and end of feature)
# Input a file like AvWT which has a column 'start' (position of methylated site)
# Other metadata columns do not matter
annotateTTS <- function(methyl_df, meta_df, location, size) {
  meta_df <- meta_df %>% 
    filter(Type == 'Transcription-Units')
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[(meta_df$Strand == '+' &
                      meta_df$Right-size <= position &
                      meta_df$Right+size >= position) |
                     (meta_df$Strand == '-' &
                      meta_df$Left-size <= position &
                      meta_df$Left+size >= position),]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'NoTTS'] <- 'X'
      next
    }
    SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
                                     meta_df$Right-size <= position &
                                     meta_df$Right+size >= position,]
    AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
                                         meta_df$Left-size <= position &
                                         meta_df$Left+size >= position,]
    for (i in 1:nrow(SenseTU_at_position)) {
      if(nrow(SenseTU_at_position) == 0) {
        next
      }
      methyl_df[methyl_df[[location]] == position, paste0('RelPos_+',i)] <- position-SenseTU_at_position[i,'Right']
    }
    for (i in 1:nrow(AntisenseTU_at_position)) {
      if(nrow(AntisenseTU_at_position) == 0) {
        next
      }
      methyl_df[methyl_df[[location]] == position, paste0('RelPos_-',i)] <- AntisenseTU_at_position[i,'Left']-position
    }
  }
  return(methyl_df)
}

# Define function to calculate variance by coverage
# Input df, 
# Outputs dataframe with 2 columns, coverage and variance
# Calculate variance of change in evolved methylation sites by coverage
# Generate df with coverage (min to max) and variance at that coverage
# Determine which coverages are (and are not) present in dataset
var_by_coverage <- function(dataset) {
  f_all_coverage <- rep(NA, max(dataset$Coverage_Sample))
  for (cov in 1:max(dataset$Coverage_Sample)) {
    if (cov %in% dataset$Coverage_Sample) {
      f_all_coverage[cov] <- TRUE
    } else {
      f_all_coverage[cov] <- FALSE
    }
  }
  #Create vector containing all coverage levels present
  f_coverage_levels <- rep(NA, sum(f_all_coverage))
  cov_i <- 1
  for (cov in 1:length(f_all_coverage)) {
    if (f_all_coverage[cov] == TRUE) {
      f_coverage_levels[cov_i] <- cov
      cov_i <- cov_i + 1
    }
  }
  #Calculate variance of all coverage levels
  f_coverage_dist <- as.data.frame(cbind(f_coverage_levels, rep(NA, length(f_coverage_levels))))
  names(f_coverage_dist)[1] <- "coverage"
  names(f_coverage_dist)[2] <- "variance"
  for (cov in 1:nrow(f_coverage_dist)) {
    temp <- dataset[dataset$Coverage_Sample == f_coverage_dist[cov,1],]
    f_coverage_dist[cov,2] <- var(temp["Percent_Methyl_Sample"]-temp["Ancestor_Mean"])
  }
  return(f_coverage_dist)
}

# Function to calculate sequencing depth across non-overlapping windows of arbitrary size
# Input dataframe, column with position, column with coverages, and window size
# Returns df with right position of each window and average depth within window
calculate_methyl_site_depth <- function(df, position_col, cov_col, w_size, calc_log2 = FALSE) {
  out_df <- tibble(position = integer(),
                   coverage = double())
  for (i in 1:ceiling(max(df[[position_col]])/w_size)) {
    out_df[i,1] <- i*w_size
    out_df[i,2] <- mean(df[between(df[[position_col]], (i-1)*w_size-1, i*w_size), cov_col][[1]])
  }
  if (calc_log2 == TRUE) {
    out_df$log2_coverage <- log2(out_df$coverage)
  }
  return(out_df)
}

# Calculate rolling mean of methylation values to show whole chromosome at once
# Default using MG1655 genome which is 4641652 long, but can put in other genome size
methylRollingMedian <- function(df, position_col, methyl_col, w_size, genome_size = 4641652, method = "exact") {
  tstart <- Sys.time()
  # Combine input into one df
  df <- bind_cols(position = df[[position_col]], methyl = df[[methyl_col]])
  
  if (method == "exact") {
  # Create empty df of every genomic position + extra for the end to wrap around
    all_pos <- data.frame(position = seq(1, genome_size + w_size, 1),
                          methyl = NA)
    all_pos[all_pos$position %in% df$position,]$methyl <- df$methyl
    # Take the first w_size of sites from the beginning of the chromosome and add them to the end
    all_pos[all_pos$position > genome_size, "methyl"] <- all_pos[all_pos$position <= w_size, "methyl"]
    # Calculate median sliding window 
    all_pos$med_methyl <- zoo::rollapply(all_pos$methyl, w_size + 1, 
                                         median, na.rm = TRUE, partial = TRUE,
                                         align = "center")
    out_df <- all_pos[1:genome_size,]
  }
  
  if (method == "fast") {
    nsites <- nrow(df)
    beginning_sites <- df %>% 
      filter(position <= w_size) %>% 
      mutate(position = position + genome_size)
    df <- bind_rows(df, beginning_sites)
    df <- data.matrix(df)
    # TODO sort df by position
    out_df <- tibble(position = as.numeric(rep(NA, nsites)),
                     mean_methyl = as.double(rep(NA, nsites)))
    out_df <- data.matrix(out_df)
    for (i in 1:nsites) {
      if (i %% 1000 == 0) {
        print(i)
      }
      out_df[[i,'position']] <- df[[i,'position']]
      out_df[[i,'mean_methyl']] <- median(df[df[,1] >= df[[i,1]] & df[,1] < (df[[i,1]] + w_size), 2])
    }
    out_df <- as_tibble(out_df)
  }
  print(Sys.time() - tstart)
  return(out_df)
}

# Calculate M value from B value
methylBtoM <- function(B, alpha = 0.001) {
  log2((B + alpha) / (1 - B + alpha))
  }

# Calculate B value from M value
methylMtoB <- function(M) {2^M / (1 + 2^M)}

# Function to remove positional bias in coverage
detrendCoverage <- function(coverage, position) {
  # Coverage and position are both numeric vectors
  # Output vector of coverages
  suppressMessages(
    cov_pos <- bind_cols(position, coverage)
  )
  colnames(cov_pos) <- c("position", "coverage")
  fit <- cov_pos %>% 
    MASS::rlm(coverage ~ ns(position, 3), data = .)
  cov_pos <- cov_pos %>% 
    add_residuals(fit)
  cov_pos$new_cov <- cov_pos$resid + (median(coverage) - median(cov_pos$resid))
  return(cov_pos$new_cov)
}

# Quantile normalization on coverage and then equivalent transformation on #methylated reads
# Optional detrend coverage and position
methylationNormalizeQuantiles <- function(df,
                                          alpha = 0.001,
                                          shift_pos = TRUE,
                                          normalize_position = TRUE,
                                          rescale = FALSE,
                                          plots = TRUE) {
  # df must have columns 'Position' and coverage starting with 'cov' and betas
  # starting with 'methyl'
  df_noNA <- df %>% tidyr::drop_na()
  coverage_df <- df_noNA %>% 
    dplyr::select(starts_with('cov'))
  if (normalize_position == TRUE) {
    for (i in 1:ncol(coverage_df)) {
      coverage_df[[i]] <- detrendCoverage(coverage = coverage_df[[i]],
                                          position = df_noNA[[2]])
    }
  }
  methyl_colnames <- df_noNA %>% 
    dplyr::select(starts_with('methyl')) %>% 
    colnames()
  methyl_df <- df_noNA %>% 
    dplyr::select(starts_with('methyl')) %>% 
    as.matrix()
  # Convert methylation beta values to number of methylated reads
  methyl_df_reads <- round(as.matrix(coverage_df)*methyl_df*100, digits = -2)/100
  # Normalize coverage by quantile for each sample
  coverage_df_normalized <- coverage_df %>% 
    as.matrix() %>% 
    preprocessCore::normalize.quantiles(x = .)
  colnames(coverage_df_normalized) <- colnames(coverage_df)
  colnames(methyl_df_reads) <- methyl_colnames
  # For each sample, fit lm to normalized coverage ~ original coverage,
  # then extract the intercept and slope from the lm object,
  # then use formula y = a + bx to normalize each sample's # methylated reads
  transformabx <- function(a, b, x) {a + b*x}
  methyl_df_normalized <- matrix(data = numeric(length = 1),
                                 nrow = nrow(methyl_df),
                                 ncol = ncol(methyl_df))
  colnames(methyl_df_normalized) <- methyl_colnames
  models <- vector('list', ncol(methyl_df))
  if (plots == TRUE) {
    model_plots <- vector('list', ncol(methyl_df))
    cov_plots <- vector('list', ncol(methyl_df))
    methyl_plots <- vector('list', ncol(methyl_df))
    methyl_read_plots <- vector('list', ncol(methyl_df))
  }
  r_squared_cov <- vector('list', ncol(methyl_df))
  r_squared_methyl <- vector('list', ncol(methyl_df))
  # Shows RMSE between beta==norm_beta (theoretical "perfect" transformation)
  # and actual best fit line
  rmse_list <- vector('list', ncol(methyl_df))

  for (i in 1:ncol(methyl_df)) {
    fit <- lm(coverage_df_normalized[,i] ~ as.matrix(coverage_df)[,i])
    models[[i]] <- fit
    r_squared_cov[[i]] <- summary(fit)$r.squared
    a <- fit$coefficients[1]
    b <- fit$coefficients[2]
    methyl_df_normalized[,i] <- 
      transformabx(a, b, methyl_df_reads[,i]) / 
      (coverage_df_normalized[,i] + alpha)
    rmse_fit <- lm(methyl_df_normalized[,i] ~ methyl_df[,i])
    rmse_pred <- seq(0, 1, 1 / (length(methyl_df[,i]) - 1))
    rmse_list[[i]] <- rmse(rmse_fit$coefficients[2]*rmse_pred + rmse_fit$coefficients[1],
                           rmse_pred)
    r_squared_methyl[[i]] <- summary(rmse_fit)$r.squared
    if (rescale == TRUE) {
      # Reassign all values below 0 to 0 and all above max(1, 99th percentile) to max(1, 99th percentile)
      # to the 99th percentile number 
      methyl_df_normalized[methyl_df_normalized[,i] < 0, i] <- 0
      methyl_df_normalized[methyl_df_normalized[,i] > 
                             max(1, quantile(methyl_df_normalized[,i], 0.99)), i] <- max(1, quantile(methyl_df_normalized[,i], 0.99))
      methyl_df_normalized[,i] <- methyl_df_normalized[,i] / max(methyl_df_normalized[,i])
    }
    if (plots == TRUE) {
      model_plots[[i]] <- autoplot(fit)
      cov_plots[[i]] <- ggplot(data = data.frame(coverage = as.matrix(coverage_df)[,i],
                                                 norm_coverage = coverage_df_normalized[,i]),
                               aes(x = coverage,
                                   y = norm_coverage)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
      methyl_plots[[i]] <- ggplot(data = data.frame(beta = methyl_df[,i],
                                                    norm_beta = methyl_df_normalized[,i]),
                                  aes(x = beta,
                                      y = norm_beta)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
      methyl_read_plots[[i]] <- ggplot(data = data.frame(methyl_reads = methyl_df_reads[,i],
                                                         norm_methyl_reads = transformabx(a, b, methyl_df_reads[,i])),
                                       aes(x = methyl_reads,
                                           y = norm_methyl_reads)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
    }
  }
  if (shift_pos == TRUE) {
    methyl_df_normalized = methyl_df_normalized + abs(min(methyl_df_normalized,
                                                          na.rm = T))
  }
  suppressMessages(
    out_df <- bind_cols(df_noNA$Position, as.data.frame(coverage_df_normalized), as.data.frame(methyl_df_normalized))
  )
  colnames(out_df)[1] <- 'Position'
  if (plots == TRUE) {
    return(list(data = out_df,
                models = models,
                model_plots = model_plots,
                cov_plots = cov_plots,
                r_squared_cov = r_squared_cov,
                r_squared_methyl = r_squared_methyl,
                methyl_plots = methyl_plots,
                methyl_read_plots = methyl_read_plots,
                RMSE = rmse_list))
  } else {
    return(list(data = out_df,
                models = models,
                r_squared_cov = r_squared_cov,
                r_squared_methyl = r_squared_methyl,
                RMSE = rmse_list))
  }

}

# Define custom theme for ggplots
theme_stone <- function() {
  theme_bw() %+replace%
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(linewidth = rel(2),
                                          fill = NA),
          panel.ontop = TRUE,
          strip.background = element_rect(linewidth = rel(2)))
}


# Build methylRawList object from dataframe
# Must have position column "Position" and coverage and beta values columns
# that start with "cov" and "methyl"
buildMethylRawList <- function(df,
                               assembly = "NC_000913.3",
                               context = "GATC",
                               resolution = "base") {
  coverage_df <- df %>% 
    dplyr::select(starts_with("cov"))
  methyl_df <- df %>% 
    dplyr::select(starts_with("methyl"))
  n_samples <- ncol(coverage_df)
  sample_ids <- strsplit(colnames(coverage_df), split = "_")
  sample_ids <- matrix(unlist(sample_ids), ncol = n_samples)[2,]
  methylRawList_obj <- new("methylRawList")
  
  for (i in 1:n_samples) {
    methylRaw_obj <- new("methylRaw")
    methylRaw_obj@sample.id = sample_ids[i]
    methylRaw_obj@assembly = assembly
    methylRaw_obj@context = context
    methylRaw_obj@resolution = resolution
    methylRaw_obj@.Data = list(rep(assembly, nrow(df)),
                               df$Position,
                               df$Position,
                               c(rep(c("+", "-"), floor(nrow(df / 2))), "+"),
                               coverage_df[,i],
                               methyl_df[,i] * coverage_df[,i],
                               (1 - methyl_df[,i]) * coverage_df[,i])
    methylRaw_obj@names = c("chr", "start", "end", "strand", "coverage", "numCs", "numTs")
    methylRaw_obj@row.names = rownames(df)
    methylRawList_obj[[i]] <- methylRaw_obj
  }
  
  return(methylRawList_obj)
}

methylationNormalize <- function(df,
                                 alpha = 0.001,
                                 normalize_position = TRUE,
                                 rescale = FALSE,
                                 plots = TRUE) {
  # df must have columns 'Position' and coverage starting with 'cov' and betas
  # starting with 'methyl'
  df_noNA <- df %>% tidyr::drop_na()
  coverage_df <- df_noNA %>% 
    dplyr::select(starts_with('cov'))
  if (normalize_position == TRUE) {
    for (i in 1:ncol(coverage_df)) {
      coverage_df[[i]] <- detrendCoverage(coverage = coverage_df[[i]],
                                          position = df_noNA[[2]])
    }
  }
  methyl_colnames <- df_noNA %>% 
    dplyr::select(starts_with('methyl')) %>% 
    colnames()
  methyl_df <- df_noNA %>% 
    dplyr::select(starts_with('methyl')) %>% 
    as.matrix()
  # Convert methylation beta values to number of methylated reads
  methyl_df_reads <- round(as.matrix(coverage_df)*methyl_df*100, digits = -2)/100
  # Normalize coverage by quantile for each sample
  coverage_df_quantnormalized <- coverage_df %>% 
    as.matrix() %>% 
    preprocessCore::normalize.quantiles(x = .)
  colnames(coverage_df_quantnormalized) <- colnames(coverage_df)
  colnames(methyl_df_reads) <- methyl_colnames
  # For each sample, fit lm to normalized coverage ~ original coverage,
  # then extract the intercept and slope from the lm object,
  # then use formula y = a + bx to normalize each sample's # methylated reads
  transformabx <- function(a, b, x) {a + b*x}
  coverage_df_normalized <- matrix(data = numeric(length = 1),
                                   nrow = nrow(coverage_df),
                                   ncol = ncol(coverage_df))
  colnames(coverage_df_normalized) <- colnames(coverage_df)
  methyl_df_normalized <- matrix(data = numeric(length = 1),
                                 nrow = nrow(methyl_df),
                                 ncol = ncol(methyl_df))
  colnames(methyl_df_normalized) <- methyl_colnames
  models <- vector('list', ncol(methyl_df))
  if (plots == TRUE) {
    model_plots <- vector('list', ncol(methyl_df))
    cov_plots <- vector('list', ncol(methyl_df))
    methyl_plots <- vector('list', ncol(methyl_df))
    methyl_read_plots <- vector('list', ncol(methyl_df))
  }
  r_squared_cov <- vector('list', ncol(methyl_df))
  r_squared_methyl <- vector('list', ncol(methyl_df))
  # Shows RMSE between beta==norm_beta (theoretical "perfect" transformation)
  # and actual best fit line
  rmse_list <- vector('list', ncol(methyl_df))
  
  for (i in 1:ncol(methyl_df)) {
    fit <- lm(coverage_df_quantnormalized[,i] ~ as.matrix(coverage_df)[,i])
    models[[i]] <- fit
    r_squared_cov[[i]] <- summary(fit)$r.squared
    a <- fit$coefficients[1]
    b <- fit$coefficients[2]
    coverage_df_normalized[,i] <- transformabx(a, b, as.matrix(coverage_df)[,i])
    methyl_df_normalized[,i] <- 
      transformabx(a, b, methyl_df_reads[,i]) / 
      (coverage_df_normalized[,i] + alpha)
    rmse_fit <- lm(methyl_df_normalized[,i] ~ methyl_df[,i])
    rmse_pred <- seq(0, 1, 1 / (length(methyl_df[,i]) - 1))
    rmse_list[[i]] <- rmse(rmse_fit$coefficients[2]*rmse_pred + rmse_fit$coefficients[1],
                           rmse_pred)
    r_squared_methyl[[i]] <- summary(rmse_fit)$r.squared
    if (rescale == TRUE) {
      methyl_df_normalized[,i] <- (methyl_df_normalized[,i] - rmse_fit$coefficients[1]) / rmse_fit$coefficients[2]
      methyl_df_normalized[methyl_df_normalized[,i] < 0, i] <- 0
      methyl_df_normalized[methyl_df_normalized[,i] > 1, i] <- 1
    }
    if (plots == TRUE) {
      model_plots[[i]] <- autoplot(fit)
      cov_plots[[i]] <- ggplot(data = data.frame(coverage = as.matrix(coverage_df)[,i],
                                                 norm_coverage = coverage_df_normalized[,i]),
                               aes(x = coverage,
                                   y = norm_coverage)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
      methyl_plots[[i]] <- ggplot(data = data.frame(beta = methyl_df[,i],
                                                    norm_beta = methyl_df_normalized[,i]),
                                  aes(x = beta,
                                      y = norm_beta)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
      methyl_read_plots[[i]] <- ggplot(data = data.frame(methyl_reads = methyl_df_reads[,i],
                                                         norm_methyl_reads = transformabx(a, b, methyl_df_reads[,i])),
                                       aes(x = methyl_reads,
                                           y = norm_methyl_reads)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = 'gray') +
        geom_smooth(method = 'lm', color = 'red')
    }
  }
  suppressMessages(
    out_df <- bind_cols(df_noNA$Position, as.data.frame(coverage_df_normalized), as.data.frame(methyl_df_normalized))
  )
  colnames(out_df)[1] <- 'Position'
  if (plots == TRUE) {
    return(list(data = out_df,
                models = models,
                model_plots = model_plots,
                cov_plots = cov_plots,
                r_squared_cov = r_squared_cov,
                r_squared_methyl = r_squared_methyl,
                methyl_plots = methyl_plots,
                methyl_read_plots = methyl_read_plots,
                RMSE = rmse_list))
  } else {
    return(list(data = out_df,
                models = models,
                r_squared_cov = r_squared_cov,
                r_squared_methyl = r_squared_methyl,
                RMSE = rmse_list))
  }
  
}


# Simulate methylation random walk
methylSimRandom <- function(starting_df = NULL,
                            methyl_distribution,
                            generations,
                            epimutation_rate) {
  if (is.null(starting_df)) {
    starting_df <- sample(methyl_distribution,
                          length(methyl_distribution),
                          replace = TRUE)
  }
  gens <- lapply(1:generations,
                 function(x) sample(c(FALSE,TRUE),
                                    size = length(starting_df),
                                    replace = TRUE,
                                    prob = c(1 - epimutation_rate,
                                             epimutation_rate)))
  for (i in 1:generations) {
    starting_df[gens[[i]]] <- sample(methyl_distribution,
                                     length(starting_df[gens[[i]]]),
                                     replace = TRUE)
  }
  return(starting_df)
}

# Calculate mutation distance from nearest GATC site
#distanceMutationGATC <- function(methyl_df, mutation_vector, 
#                                 position_col, methyl_col, genome_size = 4641652)  {
#  position_vector <- methyl_df[[position_col]]
#  methyl_vector <- methyl_df[[methyl_col]]
#  # min(abs(methyl position - mutation position)) for each mutation position
#  # or something like that
#  
#}

#annotateTSS <- function(methyl_df, meta_df, location, size) {
#  meta_df <- meta_df %>% 
#    filter(Type == 'Transcription-Units')
#  for (position in methyl_df[[location]]) {
#    if (nrow(meta_df[(meta_df$Strand == '+' &
#                      meta_df$Left - size <= position &
#                      meta_df$Left + size >= position) |
#                     (meta_df$Strand == '-' &
#                      meta_df$Right - size <= position &
#                      meta_df$Right + size >= position), ]) == 0) {
#      methyl_df[methyl_df[[location]] == position, 'NoTSS'] <- 'X'
#      next
#    }
#    SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
#                                     meta_df$Left - size <= position &
#                                     meta_df$Left + size >= position, ]
#    AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
#                                         meta_df$Right - size <= position &
#                                         meta_df$Right + size >= position, ]
#    for (i in 1:nrow(SenseTU_at_position)) {
#      if (nrow(SenseTU_at_position) == 0) {
#        next
#      }
#      methyl_df[methyl_df[[location]] == position, paste0('RelPos_+', i)] <-
#        position - SenseTU_at_position[i, 'Left']
#    }
#    for (i in 1:nrow(AntisenseTU_at_position)) {
#      if (nrow(AntisenseTU_at_position) == 0) {
#        next
#      }
#      methyl_df[methyl_df[[location]] == position, paste0('RelPos_-', i)] <-
#        AntisenseTU_at_position[i, 'Right'] - position
#    }
#  }
#  return(methyl_df)
#}

