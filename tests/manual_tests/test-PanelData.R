test_that("testing PanelData basics", {
  
  d3 <- PanelData(dem,
                  'wbcode2',
                  'year',
                  'dem',
                  'y')
  expect_true(inherits(d3, "PanelData"))
  expect_output(print(d3))
  expect_true(attr(d3, "time.id") == "year")
  expect_true(attr(d3, "unit.id") == "wbcode2")
  expect_true(attr(d3, "treatment") == "dem")
  expect_true(attr(d3, "outcome") == "y")
  expect_true(attr(d3, "is.sorted") == TRUE)
  expect_true(attr(d3, "is.balanced") ==TRUE)
  #incorrect specification
  expect_error(PanelData(dem, "wbcode2", "year", "asdf", "y"))
  
  
})

test_that("test balancing of PanelData", {
  
  set.seed(123)  # For reproducibility
  
  # Parameters
  num_units <- 10
  max_time_periods <- 10
  
  # Function to create imbalanced panel data
  create_imbalanced_panel <- function(num_units, max_time_periods) {
    data_list <- lapply(1:num_units, function(unit_id) {
      # Randomly determine the number of time periods for each unit
      num_time_periods <- sample(1:max_time_periods, 1)
      time_ids <- 1:num_time_periods
      
      # Create data frame for each unit
      data.frame(
        unit_id = unit_id,
        time_id = time_ids,
        treatment = sample(0:1, num_time_periods, replace = TRUE),
        outcome = rnorm(num_time_periods)
      )
    })
    
    # Combine into one data frame
    do.call(rbind, data_list)
  }
  
  # Create the panel data
  panel_data <- create_imbalanced_panel(num_units, max_time_periods)
  
  d2 <- PanelData(panel_data,
                  'unit_id',
                  'time_id',
                  'treatment',
                  'outcome')
  expect_true(all(table(d2$unit_id) == 10))
  expect_true(nrow(d2) == 100)
})


test_that("test sorting", {
  # Create the panel data
  set.seed(123)  # For reproducibility
  unit_id <- rep(1:5, each = 4)  # 5 units, each observed 4 times
  time_id <- rep(1:4, times = 5)  # 4 time periods for each unit
  value <- rnorm(20)  # Random values for the observed variable
  
  # Create the treatment column (0 or 1, assigned randomly to each unit)
  treatment <- rep(sample(c(0, 1), 5, replace = TRUE), each = 4)
  
  # Create the outcome column with random data
  outcome <- rnorm(20)
  
  # Combine into a data frame
  panel_data <- data.frame(unit_id, time_id, value, treatment, outcome)
  
  # Shuffle the rows to make the data out of order
  panel_data <- panel_data[sample(nrow(panel_data)), ]
  pd <- PanelData(panel.data = panel_data, unit.id = "unit_id", time.id = "time_id", treatment = "treatment", outcome = "outcome")
  sorted_panel_data <- panel_data[order(panel_data$unit_id, panel_data$time_id), ]
  
  expect_true(all(pd[,1:2] == sorted_panel_data[, 1:2]))
})

test_that("testing PanelData Error Checking", {
  
  dem$trash <- sample(letters, nrow(dem), replace = TRUE)
  expect_error(PanelData(dem,
                  'trash',
                  'year',
                  'dem',
                  'y'))
  expect_error(PanelData(dem,
                         'wbcode2',
                         'trash',
                         'dem',
                         'y'))
  
  expect_error(PanelData(dem,
                         'wbcode2',
                         'year',
                         'trash',
                         'y'))
  
  expect_error(PanelData(dem,
                         'wbcode2',
                         'year',
                         'dem',
                         'trash'))
  
  dem$year.numeric <- as.numeric(dem$year)
  
  expect_warning(d.convert <- PanelData(dem,
            'wbcode2',
            'year.numeric',
            'dem',
            'y'))
  expect_false(is.null(attr(d.convert, 
                            "time.data.map")))
})


test_that("testing PanelData methods", {
  dem$rdata <- rnorm(nrow(dem))
  d <- PanelData(dem, "wbcode2", "year", "dem", "y")
  expect_output(print(d))
  dt <- summary(d)
  expect_true(all(dim(dt) == c(1,2)))
  expect_true(all(colnames(dt) == c("num.units", "num.periods")))
  plot(d)
  plot(d, plotting.variable = "rdata")
})