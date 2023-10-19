#' @include Algorithm.SDM.R Ensemble.SDM.R Stacked.SDM.R
#' @importFrom raster stack crop extent aggregate reclassify
#' @importFrom shiny h1 h2 h3 p reactiveValues observeEvent icon brushOpts shinyApp htmlOutput
#' @importFrom shiny fluidPage fluidRow tabPanel actionButton
#' @importFrom shiny selectInput sliderInput radioButtons checkboxInput checkboxGroupInput
#' @importFrom shiny textOutput uiOutput plotOutput dataTableOutput tableOutput
#' @importFrom shiny renderText renderUI renderPlot renderDataTable renderTable
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody box tabBox
#' @importFrom shinydashboard sidebarMenu sidebarMenuOutput renderMenu menuItem menuSubItem tabItems tabItem
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient theme_minimal theme element_blank element_text
#' @importFrom reshape2 melt
#' @importFrom scales muted
#' @importFrom leaflet leafletOutput renderLeaflet colorNumeric leaflet addTiles addRasterImage addLegend
#' @importFrom sf st_as_sf
#' @importFrom raster stack crop extent aggregate reclassify
NULL

#'Plot SDMs, ensemble SDMs, and SSDMs
#'
#'Allows to plot S4 \linkS4class{Algorithm.SDM}, \linkS4class{Ensemble.SDM} and
#'\linkS4class{Stacked.SDM} class objects.
#'
#'@param x Object to be plotted (S4 Algorithm.SDM, Ensemble.SDM or Stacked.SDM
#'  object).
#'@param y,... Plot-based parameter not used.
#'
#'@return Open a window with a shiny app rendering all the results (habitat
#'  suitability map, binary map, evaluation table, variable importance and/or
#'  between-algorithm variance map, and/or algorithm evaluation,  and/or
#'  algorithm correlation matrix and/or local species richness map) in a
#'  user-friendly interface.
#'
#'@name plot.model
#'
NULL

#' @rdname plot.model
#' @export
setMethod('plot', 'Stacked.SDM', function(x, y, ...) {
  choices = list()
  for (i in seq_len(length(x@esdms))) {choices[[i]] = strsplit(x@esdms[[i]]@name, '.Ensemble.SDM', fixed = TRUE)[[1]][1]}
  if (inherits(x@esdms[[1]], 'Algorithm.SDM')) {full = FALSE} else {full = TRUE}

  ui <- dashboardPage(
    dashboardHeader(title = x@name, titleWidth = 450),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Stacked species", tabName = "stack", icon = icon("dashboard")),
        menuItem("Ensemble model", tabName = "esdm", icon = icon("pagelines")),
        selectInput('esdmchoice', 'Species:', choices, selected = NULL, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      )
    ),
    dashboardBody(
      tabItems(
        # Main page beginning #
        tabItem(tabName = "stack",
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel(htmlOutput('diversity.title'),
                                  leaflet::leafletOutput("Diversity"),
                                  textOutput('diversity.info'),
                                  title = 'Local species richness'),
                         tabPanel(leaflet::leafletOutput('endemism'),
                                  title = 'Endemism map',
                                  textOutput('endemism.info')),
                         tabPanel(leaflet::leafletOutput('uncertainty'),
                                  title = 'Uncertainty'),
                         tabPanel(tableOutput('summary'),
                                  title = 'Summary')
                  ),
                  tabBox(title = 'Variable importance',
                         tabPanel(plotOutput('varimp.barplot'),
                                  textOutput('varimp.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('varimp.table'), title = 'Table'),
                         tabPanel(tableOutput('varimplegend'), title = 'Legend')
                  )
                ),
                fluidRow(
                  tabBox(title = 'Model evaluation',
                         tabPanel(plotOutput('evaluation.barplot'),
                                  textOutput('evaluation.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('evaluation.table'), title = 'Table')
                  ),
                  if (length(x@algorithm.correlation) > 0) {
                    tabBox(title = 'Algorithm correlation',
                           tabPanel(plotOutput('algo.corr.heatmap'), title = 'Heatmap'),
                           tabPanel(tableOutput('algo.corr.table'), title = 'Table')
                    )
                  }
                )
        ),
        # Main page end #
        # ESDM page beginning #
        tabItem(tabName = 'esdm',
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel(
                           htmlOutput('probability.title'),
                           leaflet::leafletOutput('probability'),
                                   title = 'Habitat suitability'),
                         tabPanel(
                           leaflet::leafletOutput('niche'),
                                  textOutput('esdm.binary.info'),
                                  title = 'Binary map'),
                         tabPanel(leaflet::leafletOutput('esdm.uncertainty'), title = 'Uncertainty'),
                         tabPanel(tableOutput('esdm.summary'), title = 'Summary')
                  ),
                  tabBox(title = 'Variable importance',
                         tabPanel(plotOutput('esdm.varimp.barplot'),
                                  textOutput('esdm.varimp.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('esdm.varimp.table'), title = 'Table'),
                         tabPanel(tableOutput('esdmvarimplegend'), title = 'Legend')
                  )
                ),
                fluidRow(
                  tabBox(title = 'Model evaluation',
                         tabPanel(plotOutput('esdm.evaluation.barplot'),
                                  textOutput('esdm.evaluation.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('esdm.evaluation.table'), title = 'Table')
                  ),

                  tabBox(title = 'Algorithm correlation',
                         tabPanel(plotOutput('esdm.algo.corr.heatmap'), title = 'Heatmap'),
                         tabPanel(tableOutput('esdm.algo.corr.table'), title = 'Table')
                  )
                )
        )
        # esdm page ending #
      )
    )
  )

  server <- function(input, output) {

    # Main page beginning #
    # Maps
    output$diversity.title <- renderText({paste0('<b>Mean species richness error: ', round(x@evaluation['mean','species.richness.error'], 3), "<br>")})
    output$Diversity <- leaflet::renderLeaflet({
      diversity <- x@diversity.map
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(diversity), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(diversity, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(diversity), title = "Diversity")
    })
    output$endemism <- leaflet::renderLeaflet({
      endemism = x@endemism.map
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(endemism), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(endemism, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(endemism), title = "Endemism")
    })
    output$uncertainty <- leaflet::renderLeaflet({
      uncert = x@uncertainty
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(uncert), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(uncert, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(uncert), title = "Uncertainty")
      })
    # Evaluation
    output$evaluation.barplot <- renderPlot({
      evaluation <- as.matrix(x@evaluation['mean',])
      colnames(evaluation)[1:2] <- c('richness', 'success')
      barplot(evaluation, col = rainbow(length(evaluation)), beside=TRUE)
    })
    output$evaluation.table <- renderTable({
      evaluation <- x@evaluation
      names(evaluation)[1:2] <- c('richness', 'success')
      evaluation
    })

    # Algorithms correlation
    if (length(x@algorithm.correlation) > 0) {
      output$algo.corr.table <- renderTable({x@algorithm.correlation})
      output$algo.corr.heatmap <- renderPlot({
        m <- as.matrix(x@algorithm.correlation)
        # heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
        #           cellnote = round(m,2), notecol = "black", notecex = 2,
        #           trace = "none", key = FALSE, margins = c(7, 11), na.rm = TRUE,
        #           col = rev(heat.colors(1000)))
        .heatmap(m)
      })
    }
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@variable.importance))
      names(varimp) = 'Axes.evaluation'
      bar = barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)),
                    ylim = c(0,(max(varimp$Axes.evaluation)+max((varimp[2])))),
                    las = 2, ylab = 'Variable relative contribution (%)')
      arrows(bar,varimp$Axes.evaluation+varimp[,2], bar, varimp$Axes.evaluation-varimp[,2], angle=90, code=3, length=0.1)
    })
    output$varimp.table <- renderTable({x@variable.importance})
    output$varimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(x@variable.importance)), 'Variable' = names(x@variable.importance))})
    # Parameters
    output$summary <- renderTable({
      summary = data.frame(matrix(nrow = 7, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Final number of species', 'Original algorithms', 'Number of repetitions',
                             'Pseudo-absences selection', 'Cross validation method', 'Cross validation parameters')
      algo.info = character()
      for (i in seq_len(length(strsplit(x@parameters$algorithms, '.', fixed = TRUE)[[1]][-1]))) {
        algo.info = paste(algo.info, strsplit(x@parameters$algorithms, '.', fixed = TRUE)[[1]][-1][i])
      }
      if (x@parameters$PA) {PA = 'default'}
      if(x@parameters$cv == 'LOO') {cv.param = 'None'}
      if(x@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                         strsplit(x@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                         'rep =',
                                                         strsplit(x@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
      if(x@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                        strsplit(x@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                        'rep =',
                                                        strsplit(x@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
      summary$Summary = c(x@parameters$data, length(x@esdms), algo.info, x@parameters$rep, PA, x@parameters$cv, cv.param)
      if(!is.null(x@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = x@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$diversity.info <- renderText({
      x@parameters$method = switch(x@parameters$method,
                                   'pSSDM' = 'summing habitat suitability map probabilities.',
                                   'T' = paste('summing habitat suitability binary map after thresholding with',x@parameters$metric),
                                   'B' = paste('drawing repeatdly bernoulli repetitions with',x@parameters$rep.B))
      text = paste('Local species richness map realized by', x@parameters$method)
      text
    })
    output$endemism.info <- renderText({
      x@parameters$endemism = switch(x@parameters$endemism,
                                   'Any' = 'Endemism map not built (unactivated)',
                                   'WEI' = 'Endemism map built with the Weighted Endemism Index (WEI)',
                                   'B' = 'Endemism map built with the Corrected Weighted Endemism Index (CWEI)')
      x@parameters$endemism
    })
    output$varimp.info <- renderText({
      varimp.info = 'Variable relative contribution evaluated with '
      for (i in seq_len(length(x@parameters$axes.metric))) {
        if (i == 1) {
          varimp.info = paste(varimp.info, x@parameters$axes.metric[i])
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', x@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', x@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in seq_len(length(strsplit(x@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1]))) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        }
      }
      if (x@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
    # Main page ending #

    # ESDM beginning #
    # Maps
    output$probability.title <- renderText({
      paste0('<b>AUC :',round(x@esdms[[which(choices == input$esdmchoice)]]@evaluation$AUC,3),
             '  Kappa',round(x@esdms[[which(choices == input$esdmchoice)]]@evaluation$Kappa,3), "<br>")
      })
    output$probability <- leaflet::renderLeaflet({
      proba <- x@esdms[[which(choices == input$esdmchoice)]]@projection
      obs <- data.frame(
        X = x@esdms[[which(choices == input$esdmchoice)]]@data$X[which(x@esdms[[which(choices == input$esdmchoice)]]@data$Presence == 1)],
        Y = x@esdms[[which(choices == input$esdmchoice)]]@data$Y[which(x@esdms[[which(choices == input$esdmchoice)]]@data$Presence == 1)]
        )
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(proba), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(proba, colors = pal, opacity = 0.8) %>%
        leaflet::addCircles(data = obs, lng = ~ X, lat = ~ Y, radius = 1) %>%
        leaflet::addLegend(pal = pal, values = values(proba), title = "Habitat\nsuitability")
    })
    output$niche <- leaflet::renderLeaflet({
      niche.map = x@esdms[[which(choices == input$esdmchoice)]]@binary
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(niche.map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(niche.map, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(niche.map), title = "Habitat\nsuitability")
      })
    output$esdm.uncertainty <- leaflet::renderLeaflet({
      uncert.map = x@esdms[[which(choices == input$esdmchoice)]]@uncertainty
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(uncert.map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(uncert.map, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(uncert.map), title = "Uncertainty")
      })
    # Evaluation
    output$esdm.evaluation.barplot <- renderPlot({
      for (i in seq_len(length(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)))) {row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i]), '.SDM', fixed = TRUE)[[1]][1]}
      for (i in seq_len(length(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)))) {row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i] = tail(strsplit(as.character(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i]), '.', fixed = TRUE)[[1]], n = 1)}
      evaluation = x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
      table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
      barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
    })
    output$esdm.evaluation.table <- renderTable({
      for (i in seq_len(length(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)))) {row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i]), '.SDM')[[1]][1]}
      for (i in seq_len(length(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)))) {row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i] = tail(strsplit(as.character(row.names(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation)[i]), '.', fixed = TRUE)[[1]], n = 1)}
      x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation[c(2,4:8)]
    })
    # Algorithms correlation
    output$esdm.algo.corr.table <- renderTable({
      if (length(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation) > 0) {
        x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation[upper.tri(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation, diag = TRUE)] = NA
        x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation
      }
    })
    output$esdm.algo.corr.heatmap <- renderPlot({
      if (length(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation) > 0) {
        x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation[upper.tri(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation, diag = TRUE)] = NA
        m <- as.matrix(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.correlation)
        # heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
        #           cellnote = round(m,3), notecol = "black", notecex = 2,
        #           trace = "none", key = FALSE, margins = c(7, 11), na.rm = TRUE,
        #           col = rev(heat.colors(1000)))
        .heatmap(m)
      }
    })
    # Variable importance
    output$esdm.varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@esdms[[which(choices == input$esdmchoice)]]@variable.importance))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$esdm.varimp.table <- renderTable({x@esdms[[which(choices == input$esdmchoice)]]@variable.importance})
    output$varimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(x@esdms[[which(choices == input$esdmchoice)]]@variable.importance)), 'Variable' = names(x@esdms[[which(choices == input$esdmchoice)]]@variable.importance))})
    # Parameters
    output$esdm.summary <- renderTable({
      summary = data.frame(matrix(nrow = 8, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Occurrences number', 'Final number of species', 'Original algorithms', 'Number of repetitions',
                             'Pseudo-absences selection', 'Cross validation method', 'Cross validation parameters')
      algo.info = character()
      for (i in seq_len(length(strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$algorithms, '.', fixed = TRUE)[[1]][-1]))) {
        algo.info = paste(algo.info, strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$algorithms, '.', fixed = TRUE)[[1]][-1][i])
      }
      if (x@esdms[[which(choices == input$esdmchoice)]]@parameters$PA) {PA = 'default'}
      if (x@esdms[[which(choices == input$esdmchoice)]]@parameters$data == "presence-only data set") {
        nb.occ =  length(as.factor(x@esdms[[which(choices == input$esdmchoice)]]@data$Presence[which(x@esdms[[which(choices == input$esdmchoice)]]@data$Presence==1)])) / sum(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation$kept.model)
      } else {
        nb.occ =  length(as.factor(x@esdms[[which(choices == input$esdmchoice)]]@data$Presence)) / sum(x@esdms[[which(choices == input$esdmchoice)]]@algorithm.evaluation$kept.model)
      }
      if(x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv == 'LOO') {cv.param = 'None'}
      if(x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                         strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                         'rep =',
                                                         strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
      if(x@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                        strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                        'rep =',
                                                        strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
      summary$Summary = c(x@esdms[[which(choices == input$esdmchoice)]]@parameters$data, nb.occ, length(x@esdms), algo.info, x@esdms[[which(choices == input$esdmchoice)]]@parameters$rep, PA, x@esdms[[which(choices == input$esdmchoice)]]@parameters$cv, cv.param)
      if(!is.null(x@esdms[[which(choices == input$esdmchoice)]]@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = x@esdms[[which(choices == input$esdmchoice)]]@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$esdm.binary.info <- renderText({
      x@esdms[[which(choices == input$esdmchoice)]]@parameters$metric = switch(x@esdms[[which(choices == input$esdmchoice)]]@parameters$metric,
                                                                             'Kappa' = 'maximizing Kappa',
                                                                             'CCR' = 'maximizing correct proportion (CCR)',
                                                                             'TSS' = 'maximizing sensitivity and specificity sum (TSS)',
                                                                             'SES' = 'equalizing sensitivity and specificity',
                                                                             'LW' = 'taking the minimum occurence prediction',
                                                                             'ROC' = 'calculating the minimum ROC plot distance')
      text = paste('Binary map realized by', x@esdms[[which(choices == input$esdmchoice)]]@parameters$metric,
                   'with a final threshold of',  round(x@esdms[[which(choices == input$esdmchoice)]]@evaluation$threshold, digits = 3))
      text
    })
    output$esdm.varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in seq_len(length(x@esdms[[which(choices == input$esdmchoice)]]@parameters$axes.metric))) {
        if (i == 1) {
          varimp.info = paste(varimp.info, x@esdms[[which(choices == input$esdmchoice)]]@parameters$axes.metric[i])
        } else if (i == length(x@esdms[[which(choices == input$esdmchoice)]]@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', x@esdms[[which(choices == input$esdmchoice)]]@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', x@esdms[[which(choices == input$esdmchoice)]]@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$esdm.evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in seq_len(length(strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1]))) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],'(>',strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        } else if (i == length(x@esdms[[which(choices == input$esdmchoice)]]@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],'(>',strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],'(>',strsplit(x@esdms[[which(choices == input$esdmchoice)]]@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        }
      }
      if (x@esdms[[which(choices == input$esdmchoice)]]@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
    # ESDM end #
  }

  shinyApp(ui, server)
})

#' @rdname plot.model
#' @export
setMethod('plot', 'SDM', function(x, y, ...) {
  if (inherits(x, "Algorithm.SDM")) {
    full <- FALSE
  } else {
    full <- TRUE
  }
  ui <- dashboardPage(dashboardHeader(title = x@name, titleWidth = 450),
                      dashboardSidebar(disable = TRUE), dashboardBody(fluidRow(tabBox(title = "Maps",
                                                                                      tabPanel(htmlOutput('probability.title'),
                                                                                               leaflet::leafletOutput("probability"),
                                                                                               title = "Habitat suitability"), if (full) {
                                                                                                 tabPanel(leaflet::leafletOutput("niche"),
                                                                                                          title = "Binary map")
                                                                                               }, if (full) {
                                                                                                 tabPanel(leaflet::leafletOutput("uncertainty"), title = "uncertainty")
                                                                                               }, tabPanel(tableOutput("summary"), title = "Summary")), tabBox(title = "Variables importance",
                                                                                                                                                               tabPanel(plotOutput("varimp.barplot"), textOutput("varimp.info"),
                                                                                                                                                                        title = "Barplot"), tabPanel(tableOutput("varimp.table"), title = "Table"),
                                                                                                                                                               tabPanel(tableOutput("varimplegend"), title = "Legend"))), if (full) {
                                                                                                                                                                 fluidRow(tabBox(title = "Model evaluation", tabPanel(plotOutput("evaluation.barplot"),
                                                                                                                                                                                                                      textOutput("evaluation.info"), title = "Barplot"), tabPanel(tableOutput("evaluation.table"),
                                                                                                                                                                                                                                                                                  title = "Table")), if (length(x@algorithm.correlation) > 0) {
                                                                                                                                                                                                                                                                                    tabBox(title = "Algorithms correlation", tabPanel(plotOutput("algo.corr.heatmap"),
                                                                                                                                                                                                                                                                                                                                      title = "Heatmap"), tabPanel(tableOutput("algo.corr.table"),
                                                                                                                                                                                                                                                                                                                                                                   title = "Table"))
                                                                                                                                                                                                                                                                                  })
                                                                                                                                                               }))


  server <- function(input, output) {
    set.seed(122)

    if (full) {
      for (i in seq_len(length(row.names(x@algorithm.evaluation)))) {
        row.names(x@algorithm.evaluation)[i] <- strsplit(as.character(row.names(x@algorithm.evaluation)[i]),
                                                         ".SDM", fixed = TRUE)[[1]][1]
      }
      for (i in seq_len(length(row.names(x@algorithm.evaluation)))) {
        row.names(x@algorithm.evaluation)[i] <- tail(strsplit(as.character(row.names(x@algorithm.evaluation)[i]),
                                                              ".", fixed = TRUE), n = 1)
      }
      if (length(x@algorithm.correlation) > 0) {
        for (i in seq_len(length(row.names(x@algorithm.correlation)))) {
          row.names(x@algorithm.correlation)[i] <- strsplit(as.character(row.names(x@algorithm.correlation)[i]),
                                                            ".SDM", fixed = TRUE)[[1]][1]
        }
        for (i in seq_len(length(row.names(x@algorithm.correlation)))) {
          row.names(x@algorithm.correlation)[i] <- tail(strsplit(as.character(row.names(x@algorithm.correlation)[i]),
                                                                 ".", fixed = TRUE), n = 1)
        }
        x@algorithm.correlation[upper.tri(x@algorithm.correlation,
                                          diag = TRUE)] <- NA
        names(x@algorithm.correlation) <- row.names(x@algorithm.correlation)
      }
    }
    # Maps
    output$probability.title <- renderText({
      paste0('<b>AUC :', round(x@evaluation$AUC, 3),
             '  Kappa', round(x@evaluation$Kappa, 3), "<br>")
    })
    output$probability <- leaflet::renderLeaflet({
        proba.map <- x@projection
        obs <- data.frame(
          X = x@data$X[which(x@data$Presence == 1)],
          Y = x@data$Y[which(x@data$Presence == 1)]
        )
        pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                     values(proba.map), na.color = "transparent")
        leaflet::leaflet() %>%
          leaflet::addTiles() %>%
          leaflet::addRasterImage(proba.map, colors = pal, opacity = 0.8) %>%
          leaflet::addCircles(data = obs, lng = ~ X, lat = ~ Y, radius = 1) %>%
          leaflet::addLegend(pal = pal, values = values(proba.map), title = "Habitat\nsuitability")
    })
    output$niche <- leaflet::renderLeaflet({
        niche.map <- reclassify(x@projection, c(-Inf, x@evaluation$threshold,
                                                0, x@evaluation$threshold, Inf, 1))
        pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                     values(niche.map), na.color = "transparent")
        leaflet::leaflet() %>%
          leaflet::addTiles() %>%
          leaflet::addRasterImage(niche.map, colors = pal, opacity = 0.8) %>%
          leaflet::addLegend(pal = pal, values = values(niche.map), title = "Habitat\nsuitability")
    })
    if (full) {
      output$uncertainty <- leaflet::renderLeaflet({
          uncert.map <- x@uncertainty
          pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                       values(uncert.map), na.color = "transparent")
          leaflet::leaflet() %>%
            leaflet::addTiles() %>%
            leaflet::addRasterImage(uncert.map, colors = pal, opacity = 0.8) %>%
            leaflet::addLegend(pal = pal, values = values(uncert.map), title = "Uncertainty")
      })
      output$evaluation.barplot <- renderPlot({
        evaluation <- x@algorithm.evaluation
        if (!is.null(x@parameters$rep)) {
          evaluation$kept.model <- evaluation$kept.model/as.numeric(x@parameters$rep)
        } else {
          evaluation$kept.model <- evaluation$kept.model/max(evaluation$kept.model)
        }
        metrics <- "% kept.model"
        metrics.nb <- c(which(names(evaluation) == "kept.model"))
        for (i in seq_len(length(strsplit(x@parameters$ensemble.metric,
                                          ".", fixed = TRUE)[[1]][-1]))) {
          metrics <- c(metrics, strsplit(x@parameters$ensemble.metric,
                                         ".", fixed = TRUE)[[1]][-1][i])
          metrics.nb <- c(metrics.nb, which(names(evaluation) ==
                                              strsplit(x@parameters$ensemble.metric, ".", fixed = TRUE)[[1]][-1][i]))
        }
        table <- t(evaluation[metrics.nb])
        barplot(table, col = rainbow(length(metrics)), names.arg = row.names(evaluation),
                beside = TRUE)
        legend("bottomright", metrics, fill = rainbow(length(metrics)))
      })
      output$evaluation.table <- renderTable({
        x@algorithm.evaluation[c(2, 4:8)]
      })
      if (length(x@algorithm.correlation) > 0) {
        # Algorithms correlation
        output$algo.corr.table <- renderTable({
          x@algorithm.correlation
        })
        output$algo.corr.heatmap <- renderPlot({
          m <- as.matrix(x@algorithm.correlation)
          # heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          #           cellnote = round(m, 3), notecol = "black", notecex = 2,
          #           trace = "none", key = FALSE, margins = c(7, 11), na.rm = TRUE,
          #           col = rev(heat.colors(1000)))
          .heatmap(m)
        })
      }
    }
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp <- as.data.frame(t(x@variable.importance))
      names(varimp) <- "Axes.evaluation"
      barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)),
              las = 2)
    })
    output$varimp.table <- renderTable({
      x@variable.importance
    })
    output$varimplegend <- renderTable({
      data.frame(Abbreviation = abbreviate(names(x@variable.importance)),
                 Variable = names(x@variable.importance))
    })
    # Parameters
    output$summary <- renderTable({
      summary <- data.frame(matrix(nrow = 4, ncol = 1))
      names(summary) <- "Summary"
      row.names(summary) <- c("Occurrences type", "Pseudo-absences selection",
                              "Cross validation method", "Cross validation parameters")
      if (x@parameters$PA) {
        PA <- "default"
      }
      if (x@parameters$cv == "LOO") {
        cv.param <- "None"
      }
      if (x@parameters$cv == "holdout") {
        cv.param <- paste("fraction =", strsplit(x@parameters$cv.param,
                                                 "|", fixed = TRUE)[[1]][2], "rep =", strsplit(x@parameters$cv.param,
                                                                                               "|", fixed = TRUE)[[1]][3])
      }
      if (x@parameters$cv == "k-fold") {
        cv.param <- paste("k =", strsplit(x@parameters$cv.param, "|",
                                          fixed = TRUE)[[1]][2], "rep =", strsplit(x@parameters$cv.param,
                                                                                   "|", fixed = TRUE)[[1]][3])
      }
      summary$Summary <- c(x@parameters$data, PA, x@parameters$cv, cv.param)
      if (!is.null(x@parameters$algorithms)) {
        algo.info <- character()
        for (i in seq_len(length(strsplit(x@parameters$algorithms,
                                          ".", fixed = TRUE)[[1]][-1]))) {
          algo.info <- paste(algo.info, strsplit(x@parameters$algorithms,
                                                 ".", fixed = TRUE)[[1]][-1][i])
        }
        summary <- rbind(summary, data.frame(Summary = c(algo.info,
                                                         x@parameters$rep), row.names = c("Original algorithms",
                                                                                          "Number of repetitions")))
      }
      summary
    })
    output$varimp.info <- renderText({
      varimp.info <- "Axes evaluated with the variation of "
      for (i in seq_len(length(x@parameters$axes.metric))) {
        if (i == 1) {
          varimp.info <- paste(varimp.info, x@parameters$axes.metric[i])
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          varimp.info <- paste(varimp.info, "and", x@parameters$axes.metric[i],
                               ".")
        } else {
          varimp.info <- paste(varimp.info, ",", x@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$evaluation.info <- renderText({
      evaluation.info <- "Models evaluated with"
      for (i in seq_len(length(strsplit(x@parameters$ensemble.metric,
                                        ".", fixed = TRUE)[[1]][-1]))) {
        if (i == 1) {
          evaluation.info <- paste(evaluation.info, strsplit(x@parameters$ensemble.metric,
                                                             ".", fixed = TRUE)[[1]][-1][i], "(>", strsplit(x@parameters$ensemble.thresh,
                                                                                                            "|", fixed = TRUE)[[1]][-1][i], ")")
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          evaluation.info <- paste(evaluation.info, "and", strsplit(x@parameters$ensemble.metric,
                                                                    ".", fixed = TRUE)[[1]][-1][i], "(>", strsplit(x@parameters$ensemble.thresh,
                                                                                                                   "|", fixed = TRUE)[[1]][-1][i], ")", ".")
        } else {
          evaluation.info <- paste(evaluation.info, ",", strsplit(x@parameters$ensemble.metric,
                                                                  ".", fixed = TRUE)[[1]][-1][i], "(>", strsplit(x@parameters$ensemble.thresh,
                                                                                                                 "|", fixed = TRUE)[[1]][-1][i], ")")
        }
      }
      if (x@parameters$weight) {
        evaluation.info <- paste(evaluation.info, ", and then weighted with the previous metrics means")
      }
      evaluation.info
    })
  }

  shinyApp(ui, server)
})

.heatmap <- function(m) {
  Var1 <- NULL ; Var2 <- NULL ; value <- NULL
  ggplot(melt(m, id.vars = NULL),
         aes(Var1, Var2, fill = value, label = round(value, 2))) +
    geom_tile() +
    geom_text(col = "white") +
    scale_fill_gradient(guide = "none", na.value = "white",
                        low = muted("blue"), high = muted("red")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90))
}
