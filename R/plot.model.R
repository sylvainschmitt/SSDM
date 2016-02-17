#' @include Algorithm.SDM.R Ensemble.SDM.R Stacked.SDM.R
#' @importFrom raster stack crop extent aggregate reclassify
#' @importFrom shiny h1 h2 h3 p reactiveValues observeEvent icon brushOpts shinyApp
#' @importFrom shiny fluidPage fluidRow tabPanel actionButton
#' @importFrom shiny selectInput sliderInput radioButtons checkboxInput checkboxGroupInput
#' @importFrom shiny textOutput uiOutput plotOutput dataTableOutput tableOutput
#' @importFrom shiny renderText renderUI renderPlot renderDataTable renderTable
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody box tabBox
#' @importFrom shinydashboard sidebarMenu sidebarMenuOutput renderMenu menuItem menuSubItem tabItems tabItem
#' @importFrom gplots heatmap.2
#' @importFrom sp spplot SpatialPoints
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
  for (i in 1:length(x@enms)) {choices[[i]] = strsplit(x@enms[[i]]@name, '.', fixed = T)[[1]][1]}
  if (inherits(x@enms[[1]], 'Algorithm.SDM')) {full = F} else {full = T}

  ui <- dashboardPage(
    dashboardHeader(title = x@name, titleWidth = 450),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Stacked species", tabName = "stack", icon = icon("dashboard")),
        menuItem("Ensemble model", tabName = "enm", icon = icon("pagelines")),
        selectInput('enmchoice', 'Species:', choices, selected = NULL, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      )
    ),
    dashboardBody(
      tabItems(
        # Main page beginning #
        tabItem(tabName = "stack",
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel( actionButton('unzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                                   plotOutput('Diversity', dblclick = "plot1_dblclick", brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)),
                                   textOutput('diversity.info'),
                                   title = 'Local specie richness'),
                         tabPanel(plotOutput('endemism'), title = 'Endemism map', textOutput('endemism.info')),
                         tabPanel(plotOutput('uncertainty'), title = 'Uncertainty'),
                         tabPanel(tableOutput('summary'), title = 'Summary')
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
        # ENM page beginning #
        tabItem(tabName = 'enm',
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel( actionButton('enmunzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                                   plotOutput('probability', dblclick = "proba_dblclick", brush = brushOpts(id = "proba_brush", resetOnNew = TRUE)),
                                   title = 'Habitat suitability'),
                         tabPanel(plotOutput('niche'),
                                  textOutput('enm.binary.info'),
                                  title = 'Binary map'),
                         tabPanel(plotOutput('enm.uncertainty'), title = 'Uncertainty'),
                         tabPanel(tableOutput('enm.summary'), title = 'Summary')
                  ),
                  tabBox(title = 'Variable importance',
                         tabPanel(plotOutput('enm.varimp.barplot'),
                                  textOutput('enm.varimp.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('enm.varimp.table'), title = 'Table'),
                         tabPanel(tableOutput('enmvarimplegend'), title = 'Legend')
                  )
                ),
                fluidRow(
                  tabBox(title = 'Model evaluation',
                         tabPanel(plotOutput('enm.evaluation.barplot'),
                                  textOutput('enm.evaluation.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('enm.evaluation.table'), title = 'Table')
                  ),

                  tabBox(title = 'Algorithm correlation',
                         tabPanel(plotOutput('enm.algo.corr.heatmap'), title = 'Heatmap'),
                         tabPanel(tableOutput('enm.algo.corr.table'), title = 'Table')
                  )
                )
        )
        # ENM page ending #
      )
    )
  )

  server <- function(input, output) {

    # Main page beginning #
    # Maps
    # Single zoomable plot
    ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        if(!is.null(ranges$x)){ref = crop(x@diversity.map, c(ranges$x, ranges$y))} else {ref = x@diversity.map}
        ranges$x <- c(brush$xmin, brush$xmax) * (extent(ref)[2] - extent(ref)[1]) + extent(ref)[1]
        ranges$y <- c(brush$ymin, brush$ymax) * (extent(ref)[4] - extent(ref)[3]) + extent(ref)[3]
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    observeEvent(input$unzoom, {
      ranges$x <- NULL
      ranges$y <- NULL
    })
    eval = 'Mean'
    ensemble.metric = strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1]
    for (i in 1:length(ensemble.metric)) {
      eval = paste(eval, paste(ensemble.metric[i],':',round(x@evaluation[1,which(names(x@evaluation) == ensemble.metric[i])], digits = 3)))
      if (i < length(ensemble.metric)) {eval = paste(eval, ',')}
    }
    output$Diversity <- renderPlot({
      if (!is.null(ranges$x)) {diversity = crop(x@diversity.map, c(ranges$x, ranges$y))} else {diversity = x@diversity.map}
      spplot(diversity,
           main = eval,
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)))
    })
    output$endemism <- renderPlot({
      if (!is.null(ranges$x)) {endemism = crop(x@endemism.map, c(ranges$x, ranges$y))} else {endemism = x@endemism.map}
      spplot(endemism,
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)',
             col.regions = rev(terrain.colors(10000)))
    })
    output$uncertainty <- renderPlot({
      if (!is.null(ranges$x)) {uncert = crop(x@uncertainty, c(ranges$x, ranges$y))} else {uncert = x@uncertainty}
      spplot(uncert,
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)',
             col.regions = rev(terrain.colors(10000)))
      })
    # Evaluation
    output$evaluation.barplot <- renderPlot({
      evaluation = x@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / as.numeric(x@parameters$rep)
      metrics = '% kept.model'
      metrics.nb = c(which(names(evaluation) == 'kept.model'))
      for (i in 1:length(strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        metrics = c(metrics, strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i])
        metrics.nb = c(metrics.nb, which(names(evaluation) == strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i]))
      }
      table <- t(evaluation[metrics.nb])
      barplot(table, col = rainbow(length(metrics)), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', metrics, fill = rainbow(length(metrics)))
    })
    output$evaluation.table <- renderTable({x@algorithm.evaluation[c(2,4:8)]})
    if (length(x@algorithm.correlation) > 0) {
      # Algorithms correlation
      output$algo.corr.table <- renderTable({x@algorithm.correlation})
      output$algo.corr.heatmap <- renderPlot({
        m <- as.matrix(x@algorithm.correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,2), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T,
                  col = rev(heat.colors(1000)))
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
      for (i in 1:length(strsplit(x@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(x@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (x@parameters$PA) {PA = 'default'}
      if(x@parameters$cv == 'LOO') {cv.param = 'None'}
      if(x@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                         strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                         'rep =',
                                                         strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][3])}
      if(x@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                        strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                        'rep =',
                                                        strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][3])}
      summary$Summary = c(x@parameters$data, length(x@enms), algo.info, x@parameters$rep, PA, x@parameters$cv, cv.param)
      if(!is.null(x@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = x@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$diversity.info <- renderText({
      x@parameters$method = switch(x@parameters$method,
                                   'P' = 'summing habitat suitability map probabilities.',
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
      for (i in 1:length(x@parameters$axes.metric)) {
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
      for (i in 1:length(strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (x@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
    # Main page ending #

    # ENM beginning #
    # Maps
    # Single zoomable plot
    observeEvent(input$proba_dblclick, {
      brush <- input$proba_brush
      if (!is.null(brush)) {
        if(!is.null(ranges$x)){ref = crop(x@diversity.map, c(ranges$x, ranges$y))} else {ref = x@diversity.map}
        ranges$x <- c(brush$xmin, brush$xmax) * (extent(ref)[2] - extent(ref)[1]) + extent(ref)[1]
        ranges$y <- c(brush$ymin, brush$ymax) * (extent(ref)[4] - extent(ref)[3]) + extent(ref)[3]
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    observeEvent(input$enmunzoom, {
      ranges$x <- NULL
      ranges$y <- NULL
    })
    output$probability <- renderPlot({
      if (!is.null(ranges$x)) {proba = crop(x@enms[[which(choices == input$enmchoice)]]@projection, c(ranges$x, ranges$y))} else {proba = x@enms[[which(choices == input$enmchoice)]]@projection}
      spplot(proba,
           main = paste('AUC :',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$AUC,3),'  Kappa',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)),
           sp.layout=list(SpatialPoints(data.frame(X = x@enms[[which(choices == input$enmchoice)]]@data$X[which(x@enms[[which(choices == input$enmchoice)]]@data$Presence == 1)],
                                                   Y = x@enms[[which(choices == input$enmchoice)]]@data$Y[which(x@enms[[which(choices == input$enmchoice)]]@data$Presence == 1)])),
                          pch = 16, cex = 0.7, col = 'black'))
    })
    output$niche <- renderPlot({
      niche.map = reclassify(x@enms[[which(choices == input$enmchoice)]]@projection, c(-Inf,x@enms[[which(choices == input$enmchoice)]]@evaluation$threshold,0, x@enms[[which(choices == input$enmchoice)]]@evaluation$threshold,Inf,1))
      if (!is.null(ranges$x)) {niche.map = crop(niche.map, c(ranges$x, ranges$y))}
      spplot(niche.map,
           main = paste('AUC :',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$AUC,3),'  Kappa',round(x@enms[[which(choices == input$enmchoice)]]@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)))
      })
    output$enm.uncertainty <- renderPlot({
      if (!is.null(ranges$x)) {uncert.map = crop(x@enms[[which(choices == input$enmchoice)]]@uncertainty, c(ranges$x, ranges$y))} else {uncert.map = x@enms[[which(choices == input$enmchoice)]]@uncertainty}
      spplot(uncert.map,
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)',
             col.regions = rev(terrain.colors(10000)))
      })
    # Evaluation
    output$enm.evaluation.barplot <- renderPlot({
      for (i in 1:length(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation))) {row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i]), '.SDM', fixed = T)[[1]][1]}
      for (i in 1:length(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation))) {row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i] = tail(strsplit(as.character(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i]), '.', fixed = T)[[1]], n = 1)}
      evaluation = x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
      table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
      barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
    })
    output$enm.evaluation.table <- renderTable({
      for (i in 1:length(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation))) {row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i]), '.SDM')[[1]][1]}
      for (i in 1:length(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation))) {row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i] = tail(strsplit(as.character(row.names(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation)[i]), '.', fixed = T)[[1]], n = 1)}
      x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation[c(2,4:8)]
    })
    # Algorithms correlation
    output$enm.algo.corr.table <- renderTable({
      if (length(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation) > 0) {
        x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation[upper.tri(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation, diag = T)] = NA
        x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation
      }
    })
    output$enm.algo.corr.heatmap <- renderPlot({
      if (length(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation) > 0) {
        x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation[upper.tri(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation, diag = T)] = NA
        m <- as.matrix(x@enms[[which(choices == input$enmchoice)]]@algorithm.correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,3), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T,
                  col = rev(heat.colors(1000)))
      }
    })
    # Variable importance
    output$enm.varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@enms[[which(choices == input$enmchoice)]]@variable.importance))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$enm.varimp.table <- renderTable({x@enms[[which(choices == input$enmchoice)]]@variable.importance})
    output$varimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(x@enms[[which(choices == input$enmchoice)]]@variable.importance)), 'Variable' = names(x@enms[[which(choices == input$enmchoice)]]@variable.importance))})
    # Parameters
    output$enm.summary <- renderTable({
      summary = data.frame(matrix(nrow = 8, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Occurrences number', 'Final number of species', 'Original algorithms', 'Number of repetitions',
                             'Pseudo-absences selection', 'Cross validation method', 'Cross validation parameters')
      algo.info = character()
      for (i in 1:length(strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (x@enms[[which(choices == input$enmchoice)]]@parameters$PA) {PA = 'default'}
      if (x@enms[[which(choices == input$enmchoice)]]@parameters$data == "presence-only data set") {
        nb.occ =  length(as.factor(x@enms[[which(choices == input$enmchoice)]]@data$Presence[which(x@enms[[which(choices == input$enmchoice)]]@data$Presence==1)])) / sum(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation$kept.model)
      } else {
        nb.occ =  length(as.factor(x@enms[[which(choices == input$enmchoice)]]@data$Presence)) / sum(x@enms[[which(choices == input$enmchoice)]]@algorithm.evaluation$kept.model)
      }
      if(x@enms[[which(choices == input$enmchoice)]]@parameters$cv == 'LOO') {cv.param = 'None'}
      if(x@enms[[which(choices == input$enmchoice)]]@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                         strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                         'rep =',
                                                         strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$cv.param, '|', fixed = T)[[1]][3])}
      if(x@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                        strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                        'rep =',
                                                        strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$cv.param, '|', fixed = T)[[1]][3])}
      summary$Summary = c(x@enms[[which(choices == input$enmchoice)]]@parameters$data, nb.occ, length(x@enms), algo.info, x@enms[[which(choices == input$enmchoice)]]@parameters$rep, PA, x@enms[[which(choices == input$enmchoice)]]@parameters$cv, cv.param)
      if(!is.null(x@enms[[which(choices == input$enmchoice)]]@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = x@enms[[which(choices == input$enmchoice)]]@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$enm.binary.info <- renderText({
      x@enms[[which(choices == input$enmchoice)]]@parameters$metric = switch(x@enms[[which(choices == input$enmchoice)]]@parameters$metric,
                                                                             'Kappa' = 'maximizing Kappa',
                                                                             'CCR' = 'maximizing correct proportion (CCR)',
                                                                             'TSS' = 'maximizing sensitivity and specificity sum (TSS)',
                                                                             'SES' = 'equalizing sensitivity and specificity',
                                                                             'LW' = 'taking the minimum occurence prediction',
                                                                             'ROC' = 'calculating the minimum ROC plot distance')
      text = paste('Binary map realized by', x@enms[[which(choices == input$enmchoice)]]@parameters$metric,
                   'with a final threshold of',  round(x@enms[[which(choices == input$enmchoice)]]@evaluation$threshold, digits = 3))
      text
    })
    output$enm.varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in 1:length(x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric)) {
        if (i == 1) {
          varimp.info = paste(varimp.info, x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric[i])
        } else if (i == length(x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$enm.evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in 1:length(strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(x@enms[[which(choices == input$enmchoice)]]@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@enms[[which(choices == input$enmchoice)]]@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (x@enms[[which(choices == input$enmchoice)]]@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
    # ENM end #
  }

  shinyApp(ui, server)
})

#' @rdname plot.model
#' @export
setMethod('plot', 'SDM', function(x, y, ...) {
  if (inherits(x, 'Algorithm.SDM')) {full = F} else {full = T}
  ui <- dashboardPage(
    dashboardHeader(title = x@name, titleWidth = 450),
    dashboardSidebar(disable = T),
    dashboardBody(
      fluidRow(
        tabBox(title = 'Maps',
               tabPanel( actionButton('unzoom', 'unzoom', icon = icon('search-minus'), width = NULL, ...),
                         plotOutput('probability', dblclick = "plot1_dblclick", brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)),
                         title = 'Habitat suitability'),
               if(full) {tabPanel(plotOutput('niche'), title = 'Binary map')},
               if(full) {tabPanel(plotOutput('uncertainty'), title = 'uncertainty')},
               tabPanel(tableOutput('summary'), title = 'Summary')
        ),
        tabBox(title = 'Variables importance',
               tabPanel(plotOutput('varimp.barplot'),
                        textOutput('varimp.info'),
                        title = 'Barplot'),
               tabPanel(tableOutput('varimp.table'), title = 'Table'),
               tabPanel(tableOutput('varimplegend'), title = 'Legend')
        )
      ),
      if(full) {
        fluidRow(
          tabBox(title = 'Model evaluation',
                 tabPanel(plotOutput('evaluation.barplot'),
                          textOutput('evaluation.info'),
                          title = 'Barplot'),
                 tabPanel(tableOutput('evaluation.table'), title = 'Table')
          ),
          if (length(x@algorithm.correlation) > 0) {
            tabBox(title = 'Algorithms correlation',
                   tabPanel(plotOutput('algo.corr.heatmap'), title = 'Heatmap'),
                   tabPanel(tableOutput('algo.corr.table'), title = 'Table')
            )}
        )}
    )
  )


  server <- function(input, output) {
    set.seed(122)

    if(full) {
      for (i in 1:length(row.names(x@algorithm.evaluation))) {row.names(x@algorithm.evaluation)[i] = strsplit(as.character(row.names(x@algorithm.evaluation)[i]), '.SDM', fixed = T)[[1]][1]}
      for (i in 1:length(row.names(x@algorithm.evaluation))) {row.names(x@algorithm.evaluation)[i] = tail(strsplit(as.character(row.names(x@algorithm.evaluation)[i]), '.', fixed = T), n = 1)}
      if (length(x@algorithm.correlation) > 0) {
        for (i in 1:length(row.names(x@algorithm.correlation))) {row.names(x@algorithm.correlation)[i] = strsplit(as.character(row.names(x@algorithm.correlation)[i]), '.SDM', fixed = T)[[1]][1]}
        for (i in 1:length(row.names(x@algorithm.correlation))) {row.names(x@algorithm.correlation)[i] = tail(strsplit(as.character(row.names(x@algorithm.correlation)[i]), '.', fixed = T), n = 1)}
        x@algorithm.correlation[upper.tri(x@algorithm.correlation, diag = T)] = NA
        names(x@algorithm.correlation) = row.names(x@algorithm.correlation)
      }
    }
    # Maps

    # Single zoomable plot
    ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        if(!is.null(ranges$x)){ref = crop(x@projection, c(ranges$x, ranges$y))} else {ref = x@projection}
        ranges$x <- c(brush$xmin, brush$xmax) * (extent(ref)[2] - extent(ref)[1]) + extent(ref)[1]
        ranges$y <- c(brush$ymin, brush$ymax) * (extent(ref)[4] - extent(ref)[3]) + extent(ref)[3]
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    observeEvent(input$unzoom, {
      ranges$x <- NULL
      ranges$y <- NULL
    })
    output$probability <- renderPlot({
      if (!is.null(ranges$x)) {proba.map = crop(x@projection, c(ranges$x, ranges$y))} else {proba.map = x@projection}
      spplot(proba.map,
           main = paste('AUC :',round(x@evaluation$AUC,3),'  Kappa',round(x@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)),
           sp.layout=list(SpatialPoints(data.frame(X = x@data$X[which(x@data$Presence == 1)], Y = x@data$Y[which(x@data$Presence == 1)])),
                          pch = 16, cex = 0.7, col = 'black'))
    })
    output$niche <- renderPlot({
      if (!is.null(ranges$x)) {
        niche.map = crop(reclassify(x@projection, c(-Inf,x@evaluation$threshold,0, x@evaluation$threshold,Inf,1)), c(ranges$x, ranges$y))
      } else {niche.map = reclassify(x@projection, c(-Inf,x@evaluation$threshold,0, x@evaluation$threshold,Inf,1))}
      spplot(niche.map,
           main = paste('AUC :',round(x@evaluation$AUC,3),'  Kappa',round(x@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)))
    })
    if(full) {
      output$uncertainty <- renderPlot({
        if (!is.null(ranges$x)) {uncert.map = crop(x@uncertainty, c(ranges$x, ranges$y))} else {uncert.map = x@uncertainty}
        spplot(uncert.map,
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)',
             col.regions = rev(terrain.colors(10000)))
        })
      output$evaluation.barplot <- renderPlot({
        evaluation = x@algorithm.evaluation
        if(!is.null(x@parameters$rep)) {
          evaluation$kept.model = evaluation$kept.model / as.numeric(x@parameters$rep)
        } else {
          evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
        }
        metrics = '% kept.model'
        metrics.nb = c(which(names(evaluation) == 'kept.model'))
        for (i in 1:length(strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
          metrics = c(metrics, strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i])
          metrics.nb = c(metrics.nb, which(names(evaluation) == strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i]))
        }
        table <- t(evaluation[metrics.nb])
        barplot(table, col = rainbow(length(metrics)), names.arg = row.names(evaluation), beside=TRUE)
        legend('bottomright', metrics, fill = rainbow(length(metrics)))
      })
      output$evaluation.table <- renderTable({x@algorithm.evaluation[c(2,4:8)]})
      if (length(x@algorithm.correlation) > 0) {
        # Algorithms correlation
        output$algo.corr.table <- renderTable({x@algorithm.correlation})
        output$algo.corr.heatmap <- renderPlot({
          m <- as.matrix(x@algorithm.correlation)
          heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                    cellnote = round(m,3), notecol = "black", notecex = 2,
                    trace = "none", key = FALSE, margins = c(7, 11), na.rm = T,
                    col = rev(heat.colors(1000)))
        })
      }
    }
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(x@variable.importance))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)), las = 2)
    })
    output$varimp.table <- renderTable({x@variable.importance})
    output$varimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(x@variable.importance)), 'Variable' = names(x@variable.importance))})
    # Parameters
    output$summary <- renderTable({
      summary = data.frame(matrix(nrow = 4, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Pseudo-absences selection', 'Cross validation method', 'Cross validation parameters')
      if (x@parameters$PA) {PA = 'default'}
      if(x@parameters$cv == 'LOO') {cv.param = 'None'}
      if(x@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                         strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                         'rep =',
                                                         strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][3])}
      if(x@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                        strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                        'rep =',
                                                        strsplit(x@parameters$cv.param, '|', fixed = T)[[1]][3])}
      summary$Summary = c(x@parameters$data, PA, x@parameters$cv, cv.param)
      if(!is.null(x@parameters$algorithms)) {
        algo.info = character()
        for (i in 1:length(strsplit(x@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
          algo.info = paste(algo.info, strsplit(x@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
        }
        summary = rbind(summary,
                        data.frame(Summary = c(algo.info, x@parameters$rep)
                                   , row.names = c('Original algorithms','Number of repetitions')))
      }
      summary
    })
    output$varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in 1:length(x@parameters$axes.metric)) {
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
      for (i in 1:length(strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(x@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(x@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(x@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (x@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
  }

  shinyApp(ui, server)
})

