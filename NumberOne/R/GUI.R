#' @include load.occ.R load.var.R
#' @include Algorithm.Niche.Model.R Ensemble.Niche.Model.R Stack.Niche.Model.R
#' @include ensemble.R evaluate.R stacking.R
#' @include Modelling.R Ensemble.Modelling.R Stack.Modelling.R
#' @include save.model.R load.model.R plot.model.R
#' @import shiny
#' @import shinydashboard
NULL

#' @export
NumberOneGUI = function () {

  #### User Interface ####
  ui <- dashboardPage(
    dashboardHeader(title = 'NumberOne'),
    dashboardSidebar(
      sidebarMenu(id = 'actions',
                  menuItem('Welcome page', tabName = 'welcomepage'),
                  menuItem('Load',
                           menuSubItem('New data', tabName = 'newdata'),
                           menuSubItem('Previous model', tabName = 'previousmodel')),
                  uiOutput('modelling.menu')
      )
    ),
    dashboardBody(
      tabItems(
        tabItem('welcomepage',
                h2('Welcome in NumberOne')),
        tabItem('newdata',
                fluidPage(
                  fluidRow(
                    box(title = 'Environment variables',
                        selectInput('format', 'Format',
                                    list('.grd','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img'),
                                    multiple = T, selectize = T),
                        uiOutput('file'),
                        uiOutput('factors'),
                        checkboxGroupInput('load.var.options', 'loading options', list('Normalization'), selected = 'Normalization', inline = T),
                        actionButton('load', 'Load')
                    ),
                    box(title = 'Preview',
                        uiOutput('layer'),
                        renderText('debug'))
                        #uiOutput('env'))
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
        tabItem('previousmodel',
                h2('Load previous realized models')),
        tabItem('algom',
                h2('Algorithm Modelling')),
        tabItem('ensemblem',
                h2('Ensemble Modelling')),
        tabItem('stackm',
                h2('Stack species Modelling'))
      )
    )
  )

  #### Server ####
  server <- function(input, output) {

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
    data <- reactiveValues(Env = stack(), Occ = data.frame())
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
    })
    output$layer <- renderUI({
      if(length(data$Env@layers) > 0) {
        selectInput('layer', 'Variable', as.list(names(data$Env)), multiple = F)
      }
    })
#     output$env <- renderUI({
#       if(length(data$Env@layers) > 0) {
#         plotOutput(renderPlot(plot(aggregate(data$Env[[which(names(data$Env)) == input$layer]], 10))))
#       }
#     })
    output$debug <- renderText(input$layer)

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

  }

  #### Launching server ####
  options(shiny.launch.browser = T)
  options(shiny.maxRequestSize = 10e9)
  shinyApp(ui, server)
}
