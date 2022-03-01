#" Code for radar plotting adapted from Ricardo Bion
#' devtools::source_url("https://raw.githubusercontent.com/ricardo-bion/ggradar/master/R/CalculateGroupPath.R")
#" devtools::source_url("https://raw.githubusercontent.com/ricardo-bion/ggradar/master/R/CalculateAxisPath.R")
#' devtools::source_url("https://raw.githubusercontent.com/ricardo-bion/ggradar/master/R/funcCircleCoords.R")
#'
#' Calculate Group Path
#'
#' Converts variable values into a set of radial x-y coordinates
#'
#' @param df a dataframe with Col 1 is group ('unique' cluster / group ID of entity) and Col 2-n are  v1.value to vn.value - values (e.g. group/cluser mean or median) of variables v1 to v.n
#'
#' @return a dataframe of the calculated axis paths
#' 
#' @source 
#' Code adapted from a solution posted by Tony M to \url{http://stackoverflow.com/questions/9614433/creating-radar-chart-a-k-a-star-plot-spider-plot-using-ggplot2-in-r}.
CalculateGroupPath <- function(df) {
  path <- df[, 1]
  
  ## find increment
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / (ncol(df) - 1))
  ## create graph data frame
  graphData <- data.frame(seg = "", x = 0, y = 0)
  graphData <- graphData[-1, ]
  
  for (i in levels(path)) {
    pathData <- subset(df, df[, 1] == i)
    for (j in c(2:ncol(df))) {
      # pathData[,j]= pathData[,j]
      
      
      graphData <- rbind(graphData, data.frame(
        group = i,
        x = pathData[, j] * sin(angles[j - 1]),
        y = pathData[, j] * cos(angles[j - 1])
      ))
    }
    ## complete the path by repeating first pair of coords in the path
    graphData <- rbind(graphData, data.frame(
      group = i,
      x = pathData[, 2] * sin(angles[1]),
      y = pathData[, 2] * cos(angles[1])
    ))
  }
  # Make sure that name of first column matches that of input data (in case !="group")
  colnames(graphData)[1] <- colnames(df)[1]
  graphData$group <- factor(graphData$group, levels=levels(df[, 1]) ) # keep group order
  graphData # data frame returned by function
}

#' Calculate Axis Path
#'
#' Calculates x-y coordinates for a set of radial axes (one per variable being plotted in radar plot)
#'
#' @param var.names list of variables to be plotted on radar plot
#' @param min MININUM value required for the plotted axes (same value will be applied to all axes)
#' @param max MAXIMUM value required for the plotted axes (same value will be applied to all axes)
#'
#' @return a dataframe of the calculated axis paths
CalculateAxisPath <- function(var.names, min, max) {
  # var.names <- c("v1","v2","v3","v4","v5")
  n.vars <- length(var.names) # number of vars (axes) required
  # Cacluate required number of angles (in radians)
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / n.vars)
  # calculate vectors of min and max x+y coords
  min.x <- min * sin(angles)
  min.y <- min * cos(angles)
  max.x <- max * sin(angles)
  max.y <- max * cos(angles)
  # Combine into a set of uniquely numbered paths (one per variable)
  axisData <- NULL
  for (i in 1:n.vars) {
    a <- c(i, min.x[i], min.y[i])
    b <- c(i, max.x[i], max.y[i])
    axisData <- rbind(axisData, a, b)
  }
  # Add column names + set row names = row no. to allow conversion into a data frame
  colnames(axisData) <- c("axis.no", "x", "y")
  rownames(axisData) <- seq(1:nrow(axisData))
  # Return calculated axis paths
  as.data.frame(axisData)
}

#' Generate circle coordinates
#' 
#' Generate coordinates to draw a circle. 
#'
#' @param center coordinate for centroid
#' @param r radius
#' @param npoints number of coordinates to generate
#'
#' @return a dataframe
#' @source 
#' Adapted from Joran's response to \url{http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2}.
funcCircleCoords <- function(center = c(0, 0), r = 1, npoints = 100) {
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

plot_radar <- function(mda_1, mda_2, factor_n, add_title = FALSE) {
  
  factor_name <- paste0("Factor", factor_n)
  
  grid.min = -1
  grid.mid = 0
  grid.max = 1
  centre.y = grid.min - ((1 / 9) * (grid.max - grid.min))
  plot.extent.x.sf = 1
  plot.extent.y.sf = 1.2
  x.centre.range = 0.02 * (grid.max - centre.y)
  axis.label.offset = 1.15
  gridline.label.offset = -0.1 * (grid.max - centre.y)
  
  
  df.1 <- attributes(mda_1)$loadings %>%
    rownames_to_column("Category") %>%
    dplyr::select(Category, !!as.name(factor_name)) %>%
    rename(MICUSP = !!as.name(factor_name))
  
  
  df.2 <- attributes(mda_2)$loadings %>%
    rownames_to_column("Category") %>%
    dplyr::select(Category, !!as.name(factor_name)) %>%
    rename(BAWE = !!as.name(factor_name))
  
  plot.data <- full_join(df.1, df.2, by = "Category") %>%
    filter(MICUSP > .35 | BAWE > .35 | MICUSP < -.35 | BAWE < -.35) %>%
    column_to_rownames("Category") %>%
    t() %>%
    data.frame() %>%
    rownames_to_column("Corpus")
  
  axis.labels <- str_replace_all(colnames(plot.data[,-1]), "([a-z])([A-Z])", "\\1\n\\2")
  
  
  if(!is.factor(plot.data[, 1])) {
    plot.data[, 1] <- as.factor(as.character(plot.data[, 1]))
  }
  
  names(plot.data)[1] <- "group"
  
  var.names <- colnames(plot.data)[-1] # Short version of variable names
  # axis.labels [if supplied] is designed to hold 'long version' of variable names
  # with line-breaks indicated using \n
  
  # calculate total plot extent as radius of outer circle x a user-specifiable scaling factor
  plot.extent.x <- (grid.max + abs(centre.y)) * plot.extent.x.sf
  plot.extent.y <- (grid.max + abs(centre.y)) * plot.extent.y.sf
  
  # Check supplied data makes sense
  if (length(axis.labels) != ncol(plot.data) - 1) {
    stop("'axis.labels' contains the wrong number of axis labels", call. = FALSE)
  }
  if (min(plot.data[, -1], na.rm = T) < centre.y) {
    stop("plot.data' contains value(s) < centre.y", call. = FALSE)
  }
  if (max(plot.data[, -1], na.rm = T) > grid.max) {
    stop("'plot.data' contains value(s) > grid.max", call. = FALSE)
  }
  
  ### Convert supplied data into plottable format
  # (a) add abs(centre.y) to supplied plot data
  # [creates plot centroid of 0,0 for internal use, regardless of min. value of y
  # in user-supplied data]
  plot.data.offset <- plot.data
  plot.data.offset[, 2:ncol(plot.data)] <- plot.data[, 2:ncol(plot.data)] + abs(centre.y)
  # print(plot.data.offset)
  # (b) convert into radial coords
  group <- NULL
  group$path <- CalculateGroupPath(plot.data.offset)
  
  # print(group$path)
  # (c) Calculate coordinates required to plot radial variable axes
  axis <- NULL
  axis$path <- CalculateAxisPath(var.names, grid.min + abs(centre.y), grid.max + abs(centre.y))
  # print(axis$path)
  # (d) Create file containing axis labels + associated plotting coordinates
  # Labels
  axis$label <- data.frame(text = axis.labels, x = NA, y = NA)
  # print(axis$label)
  # axis label coordinates
  
  n.vars <- length(var.names)
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / n.vars)
  axis$label$x <- sapply(1:n.vars, function(i, x) {
    ((grid.max + abs(centre.y)) * axis.label.offset) * sin(angles[i])
  })
  axis$label$y <- sapply(1:n.vars, function(i, x) {
    ((grid.max + abs(centre.y)) * axis.label.offset) * cos(angles[i])
  })
  
  gridline <- NULL
  gridline$min$path <- funcCircleCoords(c(0, 0), grid.min + abs(centre.y), npoints = 360)
  gridline$mid$path <- funcCircleCoords(c(0, 0), grid.mid + abs(centre.y), npoints = 360)
  gridline$max$path <- funcCircleCoords(c(0, 0), grid.max + abs(centre.y), npoints = 360)
  # print(head(gridline$max$path))
  # gridline labels
  gridline$min$label <- data.frame(x = gridline.label.offset, y = grid.min + abs(centre.y),
                                   text = as.character(grid.min))
  
  gridline$max$label <- data.frame(x = gridline.label.offset, y = grid.max + abs(centre.y),
                                   text = as.character(grid.max))
  
  gridline$mid$label <- data.frame(x = gridline.label.offset, y = grid.mid + abs(centre.y),
                                   text = as.character(grid.mid))
  
  radar <- ggplot(axis$label) + xlab(NULL) + ylab(NULL) + coord_equal() +
    geom_text(data = subset(axis$label, axis$label$x < (-x.centre.range)),
              aes(x = x, y = y, label = text), size = 2.5, hjust = 1) +
    scale_x_continuous(limits = c(-1.5 * plot.extent.x, 1.5 * plot.extent.x)) +
    scale_y_continuous(limits = c(-plot.extent.y, plot.extent.y)) +
    geom_path(data = gridline$min$path, aes(x = x, y = y), 
              lty = "solid", colour = "gray80", size = .25) +
    geom_path(data = gridline$mid$path, aes(x = x, y = y),
              lty = "dashed", colour = "#21908CFF", size = .25) +
    geom_path(data = gridline$max$path, aes(x = x, y = y),
              lty = "solid", colour = "gray80", size = .25) +
    geom_text(data = subset(axis$label, abs(axis$label$x) <= x.centre.range),
              aes(x = x, y = y, label = text), size = 2.5, hjust = 0.5) +
    geom_text(data = subset(axis$label, axis$label$x > x.centre.range),
              aes(x = x, y = y, label = text), size = 2.5, hjust = 0) +
    geom_path(data = axis$path, aes(x = x, y = y, group = axis.no),
              colour = "gray80", size = .25) +
    geom_polygon(data = group$path, aes(x = x, y = y, group = group, fill = group),
                 alpha = .35) +
    geom_path(data = group$path, aes(x = x, y = y, group = group, color = group),
              size = .25) +
    geom_point(data = group$path, aes(x = x, y = y, group = group, fill = group), 
               size = .75, shape = 21, alpha = .5) +
    viridis::scale_color_viridis(discrete = T, direction = 1) +
    viridis::scale_fill_viridis(discrete = T, direction = 1) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size =4.5),
          legend.position="bottom")
  
  if(add_title == T) {radar <- radar + ggtitle(paste0("Dimension ", factor_n)) +
    theme(plot.title = element_text(size=10))
  }
  
  return(radar)
  
}

#' A simple function for producing the stick plots that are common in visualizing the location of category means along a given dimension.
#' @param mda_data An mda data.frame produced by the mda_loadings() function.
#' @param n_factor The factor to be plotted.
#' @return A stick plot showing category means long a positve/negative cline.
#' @export
stickplot_left <- function(mda_data, n_factor=1) {
  
  if (class(mda_data)[1] != "mda") stop ("Your mda_data must be an mda object.")
  scores <- attributes(mda_data)$group_means
  
  group_summ <- mda_data %>%
    group_by(group) %>%
    summarise(mean = mean(Factor1) %>% round(2),
              sd = sd(Factor1) %>% round(2)) %>%
    unite("values", mean:sd, sep = ", ") %>%
    mutate(values = str_replace(values, "(^.*?$)", "\\(\\1\\)"))
  
  scores <- scores %>%
    right_join(group_summ, by = "group") %>%
    select(group, values, everything()) %>%
    unite("group", group:values, sep = " ")
  
  factor_n <- paste0("Factor", n_factor)
  
  scores <- dplyr::mutate(scores, pos_neg = ifelse(!!as.name(factor_n) > 0, "High", "Low"))
  
  
  p1 <- ggplot2::ggplot(scores, ggplot2::aes(y = !!as.name(factor_n), x = 1, label = group, fill = pos_neg)) +
    ggplot2::geom_point(shape = 21) +
    scale_fill_manual(values = c("white", "gray20")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    ) +
    ggrepel::geom_text_repel(
      nudge_x      = 0.25,
      direction    = "y",
      hjust        = 0,
      segment.size = 0.1,
      size = 2
    ) +
    ggplot2::xlim(1, 2) +
    ggplot2::theme(legend.position = "none")
  
  return(p1)
}

stickplot_dual <- function(mda_1, mda_2, n_factor=1, title_a = "Corpus 1", title_b = "Corpus 2") {
  
  if (class(mda_1)[1] != "mda") stop ("Your mda_data must be an mda object.")
  if (class(mda_2)[1] != "mda") stop ("Your mda_data must be an mda object.")
  
  factor_n <- paste0("Factor", n_factor)
  
  scores_1 <- attributes(mda_1)$group_means
  
  summ_1 <- mda_1 %>%
    group_by(group) %>%
    summarise(mean = mean(!!as.name(factor_n)) %>% round(2) %>% format(nsmall = 2),
              sd = sd(!!as.name(factor_n)) %>% round(2) %>% format(nsmall = 2)) %>%
    unite("values", mean:sd, sep = ", ") %>%
    mutate(values = str_replace(values, "(^.*?$)", "\\(\\1\\)"))
  
  scores_1 <- scores_1 %>%
    right_join(summ_1, by = "group") %>%
    select(group, values, everything()) %>%
    unite("group", group:values, sep = " ")
  
  scores_1 <- dplyr::mutate(scores_1, pos_neg = ifelse(!!as.name(factor_n) > 0, "High", "Low"))
  
  scores_2 <- attributes(mda_2)$group_means
  
  summ_2 <- mda_2 %>%
    group_by(group) %>%
    summarise(mean = mean(!!as.name(factor_n)) %>% round(2) %>% format(nsmall = 2),
              sd = sd(!!as.name(factor_n)) %>% round(2) %>% format(nsmall = 2)) %>%
    unite("values", mean:sd, sep = ", ") %>%
    mutate(values = str_replace(values, "(^.*?$)", "\\(\\1\\)"))
  
  scores_2 <- scores_2 %>%
    right_join(summ_2, by = "group") %>%
    select(group, values, everything()) %>%
    unite("group", group:values, sep = " ")
  
  scores_2 <- dplyr::mutate(scores_2, pos_neg = ifelse(!!as.name(factor_n) > 0, "High", "Low"))
  
  max_y <- max(scores_1[factor_n], scores_2[factor_n]) + 0.25
  min_y <- min(scores_1[factor_n], scores_2[factor_n]) - 0.25
  
  p1 <- ggplot2::ggplot(scores_1, ggplot2::aes(y = !!as.name(factor_n), x = 1, label = group, fill = pos_neg)) +
    ggplot2::geom_point(shape = 21) +
    scale_fill_manual(values = c("white", "gray20")) +
    ggtitle(title_a) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    ) +
    ggrepel::geom_text_repel(
      nudge_x      = 0.25,
      direction    = "y",
      hjust        = 0,
      segment.size = 0.1,
      size = 2.5
    ) +
    ggplot2::xlim(1, 2) +
    ggplot2::ylim(min_y, max_y) +
    ggplot2::theme(plot.margin = margin(l = 0, r=-1)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.35, vjust = -.5, size = 10, face = "bold")) +
    ggplot2::theme(legend.position = "none")
  
  p2 <- ggplot2::ggplot(scores_2, ggplot2::aes(y = !!as.name(factor_n), x = 1, label = group, fill = pos_neg)) +
    ggplot2::geom_point(shape = 21) +
    scale_fill_manual(values = c("white", "gray40")) +
    scale_y_continuous(position = "right", limits = c(min_y, max_y)) +
    ggtitle(title_b) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    ) +
    ggrepel::geom_text_repel(
      nudge_x      = -.5,
      direction    = "y",
      hjust        = 1,
      segment.size = 0.1,
      size = 2.5
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::theme(plot.margin = margin(l = -1, r=0)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.65, vjust = -.5, size = 10, face = "bold")) +
    ggplot2::theme(legend.position = "none")
  
  p3 <- ggpubr::ggarrange(p1, p2, nrow = 1)
  
  return(p3)
}

