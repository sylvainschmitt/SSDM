#' @include load.occ.R load.var.R
#' @include Algorithm.Niche.Model.R Ensemble.Niche.Model.R Stack.Niche.Model.R
#' @include ensemble.R evaluate.R stacking.R
#' @include Modelling.R Ensemble.Modelling.R Stack.Modelling.R
#' @include save.model.R load.model.R plot.model.R
#' @import shiny shinydashboard raster
#' @importFrom gplots heatmap.2
#@importFrom raster stack reclassify extent
NULL

#' Number one package Global User Interface
#'
#' User interface to use the number one package
#'
#' @param browser logical. Option to plot or not the user in interface in you
#'   internet browser
#' @param maxmem numeric. Option to choose the maximum memory allocated to the
#'   user interface
#'
#' @return Open a window with a shiny app allowing to use the number one package
#'   with an easy user interface
#'
#' @details Due to the relatively important siez environmental data you should
#'   gave enough memory to the interface
#'
#' @examples
#' \dontrun{
#' NumberOneGUI()
#' }
#'
#' @export
NumberOneGUI = function (browser = T, maxmem = 10e9) {

  #### User Interface ####
  ui <- dashboardPage(
    dashboardHeader(title = 'NumberOne'),
    dashboardSidebar(
      sidebarMenu(id = 'actions',
                  menuItem('Welcome page', tabName = 'welcomepage'),
                  menuItem('Load',
                           menuSubItem('New data', tabName = 'newdata'),
                           menuSubItem('Previous model', tabName = 'previousmodel')),
                  uiOutput('modelling.menu'),
                  uiOutput('plotting.menu')
      )
    ),
    dashboardBody(
      tabItems(
        ### Welcome page ###
        tabItem('welcomepage',
                h2('Welcome in NumberOne')),

        ### Loading page ###

        ## Loading new data ##
        tabItem('newdata',
                fluidPage(
                  fluidRow(
                    box(title = 'Environment variables', height = 600,
                        selectInput('format', 'Format',
                                    list('.grd','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img'),
                                    multiple = T, selectize = T),
                        uiOutput('file'),
                        uiOutput('factors'),
                        checkboxGroupInput('load.var.options', 'loading options', list('Normalization'), selected = 'Normalization', inline = T),
                        actionButton('load', 'Load')
                    ),
                    box(title = 'Preview', height = 600,
                        uiOutput('layerchoice'),
                        plotOutput('env'))
                  ),
                  fluidRow(
                    box(title = 'Occurences table',
                        fileInput('Occ', 'Occurences', multiple = F, accept = '.csv'),
                        uiOutput('Xcol'),
                        uiOutput('Ycol'),
                        uiOutput('Spcol'),
                        checkboxInput('GeoRes', 'Geographic resampling', value = T),
                        sliderInput('reso', 'Resmpling gird coefficient', 1,10,1),
                        actionButton('load2', 'Load')),
                    box(title = 'Preview',
                        dataTableOutput('occ')
                    )
                  )
                )
        ),

        ## Loading Previous model ##
        tabItem('previousmodel',
                fluidPage(
                  fluidRow(
                    box(title = 'Load model',
                        uiOutput('directory'),
                        selectInput('model.type','Choose model type', c(' ', 'Ensemble model', 'Stack species model')),
                        actionButton('home', 'home', icon = icon('home')),
                        actionButton('back', 'back', icon = icon('long-arrow-left')),
                        actionButton('open', 'open', icon = icon('folder-open-o')),
                        actionButton('load.model', 'load', icon = icon('file'))
                    ),
                    box('Preview',
                        plotOutput('model.preview')
                    )
                  )
                )
        ),

        ### Modelling Page ####

        tabItem('algom',
                h2('Algorithm Modelling')),
        tabItem('ensemblem',
                h2('Ensemble Modelling')),
        tabItem('stackm',
                h2('Stack species Modelling')),

        ### Results Page ###
        tabItem('stack',
                fluidRow(
                  box(uiOutput('Stackname'), width = 12)
                ),
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel( actionButton('unzoom', 'unzoom', icon = icon('search-minus'), width = NULL),
                                   plotOutput('Diversity', dblclick = "plot1_dblclick", brush = brushOpts(id = "plot1_brush", resetOnNew = TRUE)),
                                   textOutput('diversity.info'),
                                   title = 'Diversity'),
                         tabPanel(plotOutput('Uncertainity'), title = 'Uncertainty'),
                         tabPanel(tableOutput('summary'), title = 'Summary')
                  ),
                  tabBox(title = 'Variable importance',
                         tabPanel(plotOutput('varimp.barplot'),
                                  textOutput('varimp.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('varimp.table'), title = 'Table')
                  )
                ),
                fluidRow(
                  tabBox(title = 'Model evaluation',
                         tabPanel(plotOutput('evaluation.barplot'),
                                  textOutput('evaluation.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('evaluation.table'), title = 'Table')
                  ),
                  tabBox(title = 'Algorithm correlation',
                         tabPanel(plotOutput('algo.corr.heatmap'), title = 'Heatmap'),
                         tabPanel(tableOutput('algo.corr.table'), title = 'Table')
                  )
                )
        ),
        tabItem('stackenm',
                fluidRow(
                  tabBox(title = 'Maps',
                         tabPanel(actionButton('enmunzoom', 'unzoom', icon = icon('search-minus'), width = NULL),
                                   plotOutput('probability', dblclick = "proba_dblclick", brush = brushOpts(id = "proba_brush", resetOnNew = TRUE)),
                                   title = 'Habitat suitability'),
                         tabPanel(plotOutput('niche'),
                                  textOutput('enm.binary.info'),
                                  title = 'Binary map'),
                         tabPanel(plotOutput('enm.uncertainity'), title = 'Uncertainty'),
                         tabPanel(tableOutput('enm.summary'), title = 'Summary')
                  ),
                  tabBox(title = 'Variable importance',
                         tabPanel(plotOutput('enm.varimp.barplot'),
                                  textOutput('enm.varimp.info'),
                                  title = 'Barplot'),
                         tabPanel(tableOutput('enm.varimp.table'), title = 'Table')
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
      )
    )
  )

  #### Server ####
  server <- function(input, output) {
    ### Server data ###
    data <- reactiveValues(Env = stack(), Occ = data.frame(), dir = getwd(),
                           AlgoModel = NULL, ENM = NULL, Stack = NULL)

    ### Load Menu ###

    ## Load new data page ##

    # Environmental variables loading
    load.var <- reactiveValues(factors = c(), formats = c(), norm = T,  vars = list())
    observeEvent(input$Env, {
      load.var$vars = list()
      for (i in 1:length(input$Env$name))
        load.var$vars[[i]] = input$Env$name[i]
    })
    output$file = renderUI({fileInput('Env', 'Environment', multiple = T, accept = input$format)})
    output$factors <- renderUI({
      selectInput('factors', 'Factors', load.var$vars, multiple = T, selectize = T)
    })
    observeEvent(input$load, {
      validate(
        need(length(load.var$vars) > 0, 'Choose environment variable files first !')
      )
      load.var$formats = c()
      for (i in 1 :length(load.var$vars)) {
        format = strsplit(load.var$vars[[i]], '.', fixed = T)[[1]][2]
        if (!(format %in% load.var$formats)) {load.var$formats = c(load.var$formats, format)}
      }
      load.var$factors = c()
      if (length(input$factors > 0)) {
        for (i in 1:length(input$factors)) {
          load.var$factors = c(load.var$factors,
                               paste0('X',gsub('/','.',
                                               as.character(input$Env$datapath[which(input$Env$name == input$factors[i])]),
                                               fixed = T)))
        }
      }
      if('Normalization' %in% input$load.var.options) {
        load.var$norm = T
      } else {
        load.var$norm = F
      }
      data$Env = withProgress(message = 'Variables loading',
                              load.var(directory = {},
                                       files = as.character(input$Env$datapath),
                                       format = load.var$formats,
                                       Norm = load.var$norm,
                                       tmp = F,
                                       factors = c(load.var$factors),
                                       verbose = F,
                                       GUI = T))
      for (i in 1 :length(load.var$vars)) {
        names(data$Env)[i] = strsplit(load.var$vars[[i]], '.', fixed = T)[[1]][1]
      }
      output$layerchoice <- renderUI({
        selectInput('layer', 'Variable', as.list(names(data$Env)), multiple = F, selectize = T)
      })
      output$env <- renderPlot({
        i = as.numeric(which(as.list(names(data$Env)) == input$layer))
        raster::image(data$Env[[i]])
        plot(data$Env[[i]],
             main = names(data$Env)[i],
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)')
      })
    })

    # Occurences loading
    load.occ <- reactiveValues(columns = c())
    observeEvent(input$Occ, {
      load.occ$columns = names(read.csv2(input$Occ$datapath))
    })
    output$Xcol <- renderUI({selectInput('Xcol', 'X column', load.occ$columns, multiple = F)})
    output$Ycol <- renderUI({selectInput('Ycol', 'Y column', load.occ$columns, multiple = F)})
    output$Spcol <- renderUI({selectInput('Spcol', 'Specie column', c('None', load.occ$columns), multiple = F)})
    observeEvent(input$load2, {
      validate(
        need(length(data$Env@layers) > 0, 'You need to load environmental variables before !'),
        need(length(input$Occ$datapath) > 0, 'Choose occurences file first !')
      )
      if (input$Spcol == 'None') {Spcol = NULL} else {Spcol = input$Spcol}
      data$Occ = withProgress(message = 'Occurences loading',
                              load.occ(directory = {},
                                       file = as.character(input$Occ$datapath),
                                       Xcol = input$Xcol,
                                       Ycol = input$Ycol,
                                       Spcol = Spcol,
                                       GeoRes = input$GeoRes,
                                       reso = max(res(data$Env@layers[[1]])) * as.numeric(input$reso),
                                       verbose = F,
                                       GUI = T))
    })
    output$occ <- renderDataTable({if(length(data$Occ) > 0) {data$Occ}})

    ## Load previous model page ##
    observeEvent(input$open, {data$dir = paste0(data$dir, '/', input$dir)})
    observeEvent(input$back, {data$dir = paste0(strsplit(data$dir, '/')[[1]][-length(strsplit(data$dir, '/')[[1]])], collapse = '/')})
    observeEvent(input$home, {data$dir = getwd()})
    observeEvent(input$load.model, {
      validate(
        need(input$model.type != ' ', 'You need to choose the model type first !')
      )
      if (input$model.type == 'Ensemble model') {data$ENM = load.enm(input$dir, directory = data$dir)}
      if (input$model.type == 'Stack species model') {
        data$Stack = withProgress(message = 'Model loading',
                                  load.stack(name = input$dir, directory = data$dir, GUI = T))
        }
      output$model.preview <- renderPlot({
        if (input$model.type == 'Ensemble model') {
          plot(data$ENM@projection,
               main = data$ENM@name,
               xlab = 'Longitude (\u02DA)',
               ylab = 'Latitude (\u02DA)',
               legend.args=list(text='Presence\nprobability', font = 3, line = 1))
        }
        if (input$model.type == 'Stack species model') {
          plot(data$Stack@diversity.map,
               main = data$Stack@name,
               xlab = 'Longitude (\u02DA)',
               ylab = 'Latitude (\u02DA)',
               legend.args=list(text='Local \nspecies \nrichness', font = 3, line = 1))
        }
      })
    })
    output$directory <- renderUI(selectInput('dir','Choose model directory', list.dirs(data$dir, full.names = F, recursive = F)))

    ### Modelling Menu ###
    output$modelling.menu <- renderUI({
      if(length(data$Occ) > 0) {
        sidebarMenu(
          menuItem('Modelling',
                   menuSubItem('Algorithm Modelling', tabName = 'algom'),
                   menuSubItem('Ensemble Modelling', tabName = 'ensemblem'),
                   menuSubItem('Stack Modelling', tabName = 'stackm'))
        )
      }
    })

    ## Algorithm Modelling ##

    ## Ensemble Modelling ##

    ## Stack Modelling ##

    ### Results Menu ###
    output$plotting.menu <- renderUI({
      if(!is.null(data$Stack)) {
        enms = list()
        for (i in 1:length(data$Stack@enms)) {enms[[i]] = strsplit(data$Stack@enms[[i]]@name, '.', fixed = T)[[1]][1]}
        sidebarMenu(
          menuItem('Results',
                   menuItem("Stacked species", tabName = "stack", icon = icon("dashboard")),
                   menuItem("Ensemble model", tabName = "stackenm", icon = icon("pagelines")),
                   selectInput('enmchoice', 'Ensemble model specie :', enms, selected = NULL, multiple = FALSE,
                               selectize = TRUE, width = NULL, size = NULL))
        )
      }
    })

    ## Stack results ##
    output$Stackname <- renderUI(h2(data$Stack@name, align = 'center'))
    ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    observeEvent(input$unzoom, {
      ranges$x <- NULL
      ranges$y <- NULL
    })
    output$Diversity <- renderPlot({
      eval = character()
      ensemble.metric = strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1]
      for (i in 1:length(ensemble.metric)) {
        eval = paste(eval, paste(ensemble.metric[i],':',round(data$Stack@evaluation[1,which(names(data$Stack@evaluation) == ensemble.metric[i])], digits = 3)))
        if (i < length(ensemble.metric)) {eval = paste(eval, ',')}
      }
      if (!is.null(ranges$x)) {diversity = crop(data$Stack@diversity.map, c(ranges$x, ranges$y))} else {diversity = data$Stack@diversity.map}
      plot(diversity,
           main = eval,
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           legend.args=list(text='Local \nspecies \nrichness', font = 3, line = 1))
    })
    output$Uncertainity <- renderPlot({
      eval = character()
      ensemble.metric = strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1]
      for (i in 1:length(ensemble.metric)) {
        eval = paste(eval, paste(ensemble.metric[i],':',round(data$Stack@evaluation[1,which(names(data$Stack@evaluation) == ensemble.metric[i])], digits = 3)))
        if (i < length(ensemble.metric)) {eval = paste(eval, ',')}
      }
      if (!is.null(ranges$x)) {uncert = crop(data$Stack@uncertainity, c(ranges$x, ranges$y))} else {uncert = data$Stack@uncertainity}
      plot(uncert, main = eval, legend.args=list(text='Models \nvariance', font = 3, line = 1))})
    # Evaluation
    output$evaluation.barplot <- renderPlot({
      evaluation = data$Stack@algorithm.evaluation
      evaluation$kept.model = evaluation$kept.model / as.numeric(data$Stack@parameters$rep)
      metrics = '% kept.model'
      metrics.nb = c(which(names(evaluation) == 'kept.model'))
      for (i in 1:length(strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        metrics = c(metrics, strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i])
        metrics.nb = c(metrics.nb, which(names(evaluation) == strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i]))
      }
      table <- t(evaluation[metrics.nb])
      barplot(table, col = rainbow(length(metrics)), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', metrics, fill = rainbow(length(metrics)))
    })
    output$evaluation.table <- renderTable({data$Stack@algorithm.evaluation[c(2,4:8)]})
    # Algorithms correlation
    output$algo.corr.table <- renderTable({data$Stack@algorithm.correlation})
    output$algo.corr.heatmap <- renderPlot({
      m <- as.matrix(data$Stack@algorithm.correlation)
      heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                cellnote = round(m,2), notecol = "black", notecex = 2,
                trace = "none", key = FALSE, margins = c(7, 11), na.rm = T)
    })
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(data$Stack@variables.importance))
      names(varimp) = 'Axes.evaluation'
      bar = barplot(varimp$Axes.evaluation, names.arg = strwrap(row.names(varimp)),
                    ylim = c(0,(max(varimp$Axes.evaluation)+max((varimp[2])))),
                    las = 2, ylab = 'Variable relative contribution (%)')
      arrows(bar,varimp$Axes.evaluation+varimp[,2], bar, varimp$Axes.evaluation-varimp[,2], angle=90, code=3, length=0.1)
    })
    output$varimp.table <- renderTable({data$Stack@variables.importance})
    # Parameters
    output$summary <- renderTable({
      summary = data.frame(matrix(nrow = 5, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurences type', 'Final number of species', 'Original algorithms', 'Number of repetitions', 'Pseudo-absences selection')
      algo.info = character()
      for (i in 1:length(strsplit(data$Stack@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(data$Stack@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (data$Stack@parameters$PA) {PA = 'default'}
      summary$Summary = c(data$Stack@parameters$data, length(data$Stack@enms), algo.info, data$Stack@parameters$rep, PA)
      if(!is.null(data$Stack@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = data$Stack@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$diversity.info <- renderText({
      text = switch(data$Stack@parameters$method,
                    'P' = 'summing habitat suitability map probabilities.',
                    'T' = paste('summing habitat suitability binary map after thresholding with',data$Stack@parameters$metric),
                    'B' = paste('drawing repeatdly bernoulli repetitions with',data$Stack@parameters$rep.B))
      text = paste('Local species richness map realized by', text)
      text
    })
    output$varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in 1:length(data$Stack@parameters$axes.metric)) {
        if (i == 1) {
          varimp.info = paste(varimp.info, data$Stack@parameters$axes.metric[i])
        } else if (i == length(data$Stack@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', data$Stack@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', data$Stack@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in 1:length(strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(data$Stack@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (data$Stack@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })

    ## ENM Stack Result ##
    observeEvent(input$enmchoice, {
      enms = list()
      for (i in 1:length(data$Stack@enms)) {enms[[i]] = strsplit(data$Stack@enms[[i]]@name, '.', fixed = T)[[1]][1]}
      data$ENM = data$Stack@enms[[which(enms == input$enmchoice)]]
    })
    observeEvent(input$proba_dblclick, {
      brush <- input$proba_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

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
      if (!is.null(ranges$x)) {proba = crop(data$ENM@projection, c(ranges$x, ranges$y))} else {proba = data$ENM@projection}
      plot(proba,
           main = paste('AUC :',round(data$ENM@evaluation$AUC,3),'  Kappa',round(data$ENM@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           legend.args=list(text='Presence\nprobability', font = 3, line = 1))
      points(data$ENM@data$X[which(data$ENM@data$Presence == 1)],
             data$ENM@data$Y[which(data$ENM@data$Presence == 1)],
             pch = 16, cex = 0.7)
    })
    output$niche <- renderPlot({
      niche.map = reclassify(data$ENM@projection, c(-Inf,data$ENM@evaluation$threshold,0, data$ENM@evaluation$threshold,Inf,1))
      if (!is.null(ranges$x)) {niche.map = crop(niche.map, c(ranges$x, ranges$y))}
      plot(niche.map, main = paste('AUC :',round(data$ENM@evaluation$AUC,3),'  Kappa',round(data$ENM@evaluation$Kappa,3)))})
    output$enm.uncertainity <- renderPlot({
      if (!is.null(ranges$x)) {uncert.map = crop(data$ENM@uncertainity, c(ranges$x, ranges$y))} else {uncert.map = data$ENM@uncertainity}
      plot(uncert.map, main = paste('AUC :',round(data$ENM@evaluation$AUC,3),'  Kappa',round(data$ENM@evaluation$Kappa,3)), legend.args=list(text='Models \nvariance', font = 3, line = 1))})
    # Evaluation
    output$enm.evaluation.barplot <- renderPlot({
      evaluation = data$ENM@algorithm.evaluation
      for (i in 1:length(row.names(data$ENM@algorithm.evaluation))) {row.names(evaluation)[i] = strsplit(as.character(row.names(data$ENM@algorithm.evaluation)[i]), '.Niche')[[1]][1]}
      evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
      table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
      barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
    })
    output$enm.evaluation.table <- renderTable({
      algo.eval = data$ENM@algorithm.evaluation
      for (i in 1:length(row.names(data$ENM@algorithm.evaluation))) {row.names(algo.eval)[i] = strsplit(as.character(row.names(data$ENM@algorithm.evaluation)[i]), '.Niche')[[1]][1]}
      algo.eval[c(2,4:8)]
    })
    # Algorithms correlation
    output$enm.algo.corr.table <- renderTable({
      if (length(data$ENM@algorithm.correlation) > 0) {
        data$ENM@algorithm.correlation[upper.tri(data$ENM@algorithm.correlation, diag = T)] = NA
        data$ENM@algorithm.correlation
      }
    })
    output$enm.algo.corr.heatmap <- renderPlot({
      if (length(data$ENM@algorithm.correlation) > 0) {
        data$ENM@algorithm.correlation[upper.tri(data$ENM@algorithm.correlation, diag = T)] = NA
        m <- as.matrix(data$ENM@algorithm.correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,3), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T)
      }
    })
    # Variable importance
    output$enm.varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(data$ENM@variables.importance[-1]))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = row.names(varimp), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$enm.varimp.table <- renderTable({data$ENM@variables.importance[-1]})
    # Parameters
    output$enm.summary <- renderTable({
      summary = data.frame(matrix(nrow = 5, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurences type', 'Final number of species', 'Original algorithms', 'Number of repetitions', 'Pseudo-absences selection')
      algo.info = character()
      for (i in 1:length(strsplit(data$ENM@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(data$ENM@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (data$ENM@parameters$PA) {PA = 'default'}
      summary$Summary = c(data$ENM@parameters$data, length(data$Stack@enms), algo.info, data$ENM@parameters$rep, PA)
      if(!is.null(data$ENM@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = data$ENM@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$enm.binary.info <- renderText({
      data$ENM@parameters$metric = switch(data$ENM@parameters$metric,
                                                                             'Kappa' = 'maximizing Kappa',
                                                                             'CCR' = 'maximizing correct proportion (CCR)',
                                                                             'TSS' = 'maximizing sensitivity and specificity sum (TSS)',
                                                                             'SES' = 'equalizing sensitivity and specificity',
                                                                             'LW' = 'taking the minimum occurence prediction',
                                                                             'ROC' = 'calculating the minimum ROC plot distance')
      text = paste('Binary map realized by', data$ENM@parameters$metric,
                   'with a final threshold of',  round(data$ENM@evaluation$threshold, digits = 3))
      text
    })
    output$enm.varimp.info <- renderText({
      varimp.info = 'Axes evaluated with the variation of '
      for (i in 1:length(data$ENM@parameters$axes.metric)) {
        if (i == 1) {
          varimp.info = paste(varimp.info, data$ENM@parameters$axes.metric[i])
        } else if (i == length(data$ENM@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', data$ENM@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', data$ENM@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$enm.evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in 1:length(strsplit(data$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste(evaluation.info, strsplit(data$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(data$ENM@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(data$ENM@parameters$axes.metric) && i != 1) {
          evaluation.info = paste(evaluation.info, 'and', strsplit(data$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(data$ENM@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste(evaluation.info, ',', strsplit(data$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],'(>',strsplit(data$ENM@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (data$ENM@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })
  }

  #### Launching server ####
  #options(shiny.launch.browser = browser)
  options(shiny.maxRequestSize = maxmem)
  shinyApp(ui, server)
}
