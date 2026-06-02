# =========================================================
# Shiny App: Simulating Linear Regressions (Power Analysis)
# Requires: shiny only (all plots use base R graphics)
# =========================================================


# https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html
# Run the shinylive::export line to populate the docs folder 
# so that shinylive works from github
#shinylive::export(appdir = "../LinearModelPowerSimulation/", destdir = "docs")
#httpuv::runStaticServer("docs/", port = 8008)

max_continuous <- 2

# Colour palettes matching Dark2 brewer for groups
group_cols <- c(control = "#1B9E77", treat = "#D95F02")

# Colours for significance
sig_cols <- c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")

# ---------------------------
# Helper: run one lm and extract tidy slope row
# Robust to singular fits and missing rows
# ---------------------------
tidy_slope <- function(mod) {
  cf <- summary(mod)$coefficients
  # Drop the intercept row; keep the first non-intercept term
  non_int <- rownames(cf)[rownames(cf) != "(Intercept)"]
  if (length(non_int) == 0) {
    # Singular model - return NA row so it doesn't break rbind
    return(data.frame(term = NA_character_, estimate = NA_real_,
                      std.error = NA_real_, statistic = NA_real_,
                      p.value   = NA_real_,
                      stringsAsFactors = FALSE))
  }
  row <- cf[non_int[1], , drop = FALSE]
  data.frame(
    term      = rownames(row),
    estimate  = row[, "Estimate"],
    std.error = row[, "Std. Error"],
    statistic = row[, "t value"],
    p.value   = row[, "Pr(>|t|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ---------------------------
# UI
# ---------------------------
power_lm_ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      #run_sim {
        background-color: #2196F3;
        color: white;
        border: none;
        padding: 6px 14px;
        font-size: 14px;
        margin-top: 4px;
        .row { margin-top: 0px; }
        .col-sm-4, .col-sm-6 { padding-top: 0px; }
      }
      #run_sim:hover   { background-color: #1769aa; }
      #run_sim:active  { transform: scale(0.97); }
    "))
  ),
  
  titlePanel("Sampling Variability and Statistical Power"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("type", "Type of predictor variable",
                   c("Categorical", "Continuous"), "Categorical"),
      sliderInput("n", "Number of observations:", 2, 100, 10, 1),
      sliderInput("delta", "True difference in means:", 0, 10, 5),
      sliderInput("sd", "Standard deviation of each group:", 0, 10, 5),
      # shinyWidgets::sliderTextInput replaced with selectInput
      selectInput("p", "p level:",
                  choices = c("0.001", "0.01", "0.05", "0.1"),
                  selected = "0.05"),
      br(),
      actionButton("run_sim", "Run Simulation"),
      helpText("_____________________"),
      helpText("Glenn Tattersall, PhD"),
      helpText("For use in BIOL 3P96 - Biostatistics")
    ),
    
    mainPanel(
      fluidRow(
        column(6, h4("Population distributions",  style = "text-align:center; margin-top:0px; margin-bottom:2px;"), plotOutput("pop_plot")),
        column(6, h4(textOutput("samp_size"),     style = "text-align:center; margin-top:0px; margin-bottom:2px;"), plotOutput("sample_plot"))
      ),
      fluidRow(
        column(4, h4(textOutput("uncertainty_text"), style = "text-align:center; margin-top:0px; margin-bottom:2px;"), plotOutput("uncertainty_plot")),
        column(4, h4("Effect size of 100 trials",    style = "text-align:center; margin-top:0px; margin-bottom:2px;"), plotOutput("effect_plot")),
        column(4, h4("P-values of 100 trials",       style = "text-align:center; margin-top:0px; margin-bottom:2px;"), plotOutput("p_plot"))
      ),
      hr(),
      h4("Interpretation"),
      verbatimTextOutput("interpretation")
    )
  )
)

# ---------------------------
# SERVER
# ---------------------------
power_lm_server <- function(input, output, session) {
  
  # Dynamic slider label
  observe({
    new_label <- if (input$type == "Categorical") {
      "True difference in means:"
    } else {
      "True slope:"
    }
    updateSliderInput(session, "delta", label = new_label)
  })
  
  output$samp_size <- renderText({
    if (mods()$type == "Categorical") {
      paste0("Sample N = ", input$n, " for each treatment")
    } else {
      paste0("Sample N = ", input$n, " observations")
    }
  })
  
  output$uncertainty_text <- renderText("95% CI of 100 trials")
  
  # ---------------------------
  # Population distribution plot
  # ---------------------------
  output$pop_plot <- renderPlot({
    
    p_alpha <- as.numeric(input$p)
    
    if (input$type == "Categorical") {
      
      upper <- max(0, input$delta) + input$sd * 3
      lower <- min(0, input$delta) - input$sd * 3
      xs    <- seq(lower, upper, length.out = 201)
      
      y_ctrl  <- dnorm(xs, mean = 0,           sd = max(input$sd, 1e-6))
      y_treat <- dnorm(xs, mean = input$delta,  sd = max(input$sd, 1e-6))
      ylim    <- c(0, max(y_ctrl, y_treat) * 1.1)
      
      plot(xs, y_ctrl,
           type = "n", ylim = ylim,
           xlab = "Value", ylab = "Density",
           main = "", las = 1)
      
      # Shaded areas
      polygon(c(xs, rev(xs)), c(y_ctrl,  rep(0, length(xs))),
              col = adjustcolor(group_cols["control"], 0.3), border = NA)
      polygon(c(xs, rev(xs)), c(y_treat, rep(0, length(xs))),
              col = adjustcolor(group_cols["treat"],   0.3), border = NA)
      
      # Mean lines
      abline(v = 0,           col = group_cols["control"], lwd = 2)
      abline(v = input$delta, col = group_cols["treat"],   lwd = 2)
      
      legend("topright",
             legend = c("Control", "Treatment"),
             fill   = adjustcolor(group_cols, 0.3),
             border = group_cols,
             bty    = "n", cex = 0.9)
      
    } else {
      
      xs  <- seq(0, max_continuous, length.out = 50)
      yc  <- xs * input$delta
      sd_ <- max(input$sd, 1e-6)
      
      plot(xs, yc,
           type = "n",
           ylim = range(yc + 2 * sd_, yc - 2 * sd_),
           xlab = "Predictor", ylab = "Response",
           main = "", las = 1)
      
      # ±2 SD band
      polygon(c(xs, rev(xs)),
              c(yc + 2 * sd_, rev(yc - 2 * sd_)),
              col = adjustcolor("grey40", 0.2), border = NA)
      # ±1 SD band
      polygon(c(xs, rev(xs)),
              c(yc + sd_, rev(yc - sd_)),
              col = adjustcolor("grey40", 0.25), border = NA)
      
      lines(xs, yc, lwd = 2)
    }
  })
  
  # ---------------------------
  # Simulations (run on button click)
  # ---------------------------
  mods <- eventReactive(input$run_sim, ignoreNULL = FALSE, {
    
    # isolate() all inputs so only the button click triggers this reactive,
    # preventing Shiny caching when sliders change between clicks
    p_alpha <- isolate(as.numeric(input$p))
    n       <- isolate(input$n)
    delta   <- isolate(input$delta)
    sd_     <- isolate(input$sd)
    type_   <- isolate(input$type)
    
    # Generate 100 datasets and fit lm to each
    results_list <- lapply(seq_len(100), function(i) {
      if (type_ == "Categorical") {
        # Use an explicit factor with fixed reference level so the
        # coefficient row is always named "treatmenttreat", never ambiguous
        trt <- factor(rep(c("control", "treat"), each = n),
                      levels = c("control", "treat"))
        dat <- data.frame(
          treatment = trt,
          value     = rnorm(n * 2,
                            mean = rep(c(0, delta), each = n),
                            sd   = sd_)
        )
      } else {
        x <- runif(n, 0, max_continuous)
        dat <- data.frame(
          treatment = x,
          value     = x * delta + rnorm(n, 0, sd_)
        )
      }
      tidy_slope(lm(value ~ treatment, data = dat))
    })
    
    result <- do.call(rbind, results_list)
    # Drop any rows from singular/failed fits before sorting
    result <- result[!is.na(result$estimate), ]
    result <- result[order(result$estimate), ]
    result$n           <- seq_len(nrow(result))
    result$significant <- result$p.value < p_alpha
    
    # First simulation dataset for sample_plot
    if (type_ == "Categorical") {
      trt <- factor(rep(c("control", "treat"), each = n),
                    levels = c("control", "treat"))
      sample_dat <- data.frame(
        treatment = trt,
        value     = rnorm(n * 2,
                          mean = rep(c(0, delta), each = n),
                          sd   = sd_)
      )
    } else {
      x <- runif(n, 0, max_continuous)
      sample_dat <- data.frame(
        treatment = x,
        value     = x * delta + rnorm(n, 0, sd_)
      )
    }
    
    list(result = result, sample = sample_dat, type = type_)
  })
  
  # ---------------------------
  # P-value histogram
  # ---------------------------
  output$p_plot <- renderPlot({
    req(mods())
    df      <- mods()$result
    p_alpha <- as.numeric(input$p)
    cols    <- ifelse(df$significant, sig_cols["TRUE"], sig_cols["FALSE"])
    
    h <- hist(df$p.value, breaks = seq(0, 1, length.out = 101),
              plot = FALSE)
    
    # Colour each bar by whether its midpoint < alpha
    bar_cols <- ifelse(h$mids < p_alpha, sig_cols["TRUE"], sig_cols["FALSE"])
    
    plot(h, col = bar_cols, border = "white",
         xlim = c(0, 1),
         xlab = "P value", ylab = "Count",
         main = "", las = 1)
    abline(v = p_alpha, col = "red", lty = 2, lwd = 2)
    
    legend("topright",
           legend = c("Significant (reject H0)", "Not significant"),
           fill   = c(sig_cols["TRUE"], sig_cols["FALSE"]),
           bty    = "n", cex = 0.85)
  })
  
  # ---------------------------
  # Effect size histogram
  # ---------------------------
  output$effect_plot <- renderPlot({
    req(mods())
    df      <- mods()$result
    p_alpha <- as.numeric(input$p)
    
    xlab_txt <- if (mods()$type == "Categorical") {
      "Estimated difference in means"
    } else {
      "Estimated slope"
    }
    
    h        <- hist(df$estimate, breaks = 25, plot = FALSE)
    # Assign colour by whether trials in that bin are mostly significant
    # (simplest faithful approach: colour each bar by sign of midpoint vs delta)
    # Instead we colour by significance of each observation using density approach:
    # Use a plain two-colour split: significant obs drive a red overlay
    plot(h, col = sig_cols["FALSE"], border = "white",
         xlab = xlab_txt, ylab = "Count",
         main = "", las = 1)
    
    # Overlay significant subset
    h_sig <- hist(df$estimate[df$significant], breaks = h$breaks, plot = FALSE)
    plot(h_sig, col = sig_cols["TRUE"], border = "white", add = TRUE)
    
    abline(v = input$delta, lwd = 2)
    
    legend("topright",
           legend = c("Significant (reject H0)", "Not significant"),
           fill   = c(sig_cols["TRUE"], sig_cols["FALSE"]),
           bty    = "n", cex = 0.85)
  })
  
  # ---------------------------
  # CI / uncertainty plot
  # ---------------------------
  output$uncertainty_plot <- renderPlot({
    req(mods())
    df      <- mods()$result
    p_alpha <- as.numeric(input$p)
    mult    <- qnorm(1 - p_alpha / 2)
    
    df$xmin <- df$estimate - mult * df$std.error
    df$xmax <- df$estimate + mult * df$std.error
    
    pt_cols <- ifelse(df$significant, sig_cols["TRUE"], sig_cols["FALSE"])
    
    xlab_txt <- if (mods()$type == "Categorical") {
      "Estimated difference in means"
    } else {
      "Estimated slope"
    }
    
    xlim <- range(c(df$xmin, df$xmax))
    
    plot(df$estimate, df$n,
         col  = pt_cols, pch = 16, cex = 0.7,
         xlim = xlim,
         xlab = xlab_txt, ylab = "Simulation number",
         main = "", las = 1)
    
    # Error bars
    segments(df$xmin, df$n, df$xmax, df$n,
             col = pt_cols, lwd = 0.8)
    
    abline(v = input$delta, lwd = 2)
    abline(v = 0, lty = 2, lwd = 1.5)
  })
  
  # ---------------------------
  # Sample data plot
  # ---------------------------
  output$sample_plot <- renderPlot({
    req(mods())
    dat <- mods()$sample
    
    if (mods()$type == "Categorical") {
      
      ctrl  <- dat$value[dat$treatment == "control"]
      treat <- dat$value[dat$treatment == "treat"]
      
      # Side-by-side boxplots with jitter
      grp_list <- list(control = ctrl, treat = treat)
      cols_fill <- adjustcolor(group_cols, 0.3)
      
      bp <- boxplot(value ~ treatment, data = dat,
                    col    = cols_fill,
                    border = group_cols,
                    outline = FALSE,
                    xlab = "Treatment", ylab = "Response",
                    main = "", las = 1,
                    names = c("Control", "Treatment"))
      
      # Jitter points
      points(jitter(rep(1, length(ctrl)),  amount = 0.1), ctrl,
             col = group_cols["control"], pch = 16, cex = 0.9)
      points(jitter(rep(2, length(treat)), amount = 0.1), treat,
             col = group_cols["treat"],   pch = 16, cex = 0.9)
      
    } else {
      
      mod_s <- lm(value ~ treatment, data = dat)
      xs    <- seq(min(dat$treatment), max(dat$treatment), length.out = 100)
      pred  <- predict(mod_s,
                       newdata  = data.frame(treatment = xs),
                       interval = "confidence")
      
      plot(dat$treatment, dat$value,
           pch  = 16, col = "grey30",
           xlab = "Predictor", ylab = "Response",
           main = "", las = 1)
      
      polygon(c(xs, rev(xs)),
              c(pred[, "lwr"], rev(pred[, "upr"])),
              col = adjustcolor("steelblue", 0.25), border = NA)
      
      lines(xs, pred[, "fit"], col = "steelblue", lwd = 2)
    }
  })
  
  # ---------------------------
  # Interpretation text
  # ---------------------------
  output$interpretation <- renderText({
    req(mods())
    df        <- mods()$result
    p_alpha   <- as.numeric(input$p)
    mean_est  <- round(mean(df$estimate), 2)
    prop_sign <- round(mean(df$significant) * 100, 1)
    bias      <- round(mean_est - input$delta, 2)
    
    if (input$type == "Categorical") {
      body <- paste0(
        "Across 100 simulated trials:\n",
        "- Average estimated difference in means = ", mean_est, "\n",
        "- Bias relative to true difference = ", bias, "\n",
        "- ", prop_sign, "% of trials would lead a researcher correctly to reject H0",
        " at p < ", p_alpha, ".\n"
      )
    } else {
      body <- paste0(
        "Across 100 simulated trials:\n",
        "- Average estimated slope = ", mean_est, "\n",
        "- Bias relative to true slope = ", bias, "\n",
        "- ", prop_sign, "% of trials would lead a researcher to reject H0",
        " at p < ", p_alpha, ".\n"
      )
    }
    
    paste0(
      body,
      "\nNote: In the histogram and CI plots, colour indicates statistical decision:\n",
      "'Significant' = reject null hypothesis (p < chosen alpha),\n",
      "'Not significant' = fail to reject null hypothesis."
    )
  })
}

shinyApp(ui = power_lm_ui, server = power_lm_server)