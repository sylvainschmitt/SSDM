#' @include load_occ.R load_var.R
#' @include Algorithm.SDM.R Ensemble.SDM.R Stacked.SDM.R
#' @include ensemble.R stacking.R
#' @include modelling.R ensemble_modelling.R stack_modelling.R
#' @include save.model.R load_model.R plot.model.R
#' @importFrom sp spplot SpatialMultiPoints
#' @importFrom raster stack crop extent aggregate reclassify
#' @importFrom shiny h1 h2 h3 p reactiveValues observeEvent icon brushOpts shinyApp runApp
#' @importFrom shiny fluidPage fluidRow tabPanel actionButton validate need withProgress
#' @importFrom shiny textInput selectInput sliderInput radioButtons checkboxInput checkboxGroupInput
#' @importFrom shiny textOutput uiOutput plotOutput dataTableOutput tableOutput
#' @importFrom shiny renderText renderUI renderPlot renderDataTable renderTable
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody box tabBox
#' @importFrom shinydashboard sidebarMenu sidebarMenuOutput renderMenu menuItem menuSubItem tabItems tabItem
#' @importFrom shinyFiles shinyDirChoose shinyFileChoose shinyDirButton shinyFilesButton
#' @importFrom gplots heatmap.2
NULL

#' SSDM package Global User Interface
#'
#' User interface of the SSDM package.
#'
#' @param browser logical. Option to plot the user global interface in your
#'   internet browser.
#' @param maxmem numeric. Option to choose the maximum memory allocated to the
#'   user interface.
#' @param col character. User interface color theme. One of "blue", "black",
#'   "purple", "green", "red", or "yellow".
#'
#' @return Open a window with a shiny app to use the SSDM package with an
#'   user-friendly interface.
#'
#' @details If your environmental variables have an important size, you should
#'   gave enough memory to the interface with the (\code{maxmem} parameter).
#'
#' @examples
#' \dontrun{
#' gui()
#' }
#'
#' @export
gui = function (browser = F, maxmem = 10e9, col = 'blue') {

  #### User Interface ####
  ui <- dashboardPage(skin = col,
    dashboardHeader(title = 'SSDM'),
    dashboardSidebar(sidebarMenuOutput('menu')),
    dashboardBody(
      tabItems(
        ### Welcome page ###
        tabItem('welcomepage',
                h1('SSDM: Stacked species distribution modelling', align = 'center'),
                h1(' '),
                h2('Welcome to the user-friendly interface', align = 'center'),
                h2(' '),
                p('SSDM is a package to map species richness and endemism based on stacked species distribution models (SSDM). Individual SDMs can be created using a single or multiple algorithms (ensemble SDMs). For each species, an SDM can yield a habitat suitability map, a binary map, a between-algorithm variance map, and can assess variable importance, algorithm accuracy, and between-algorithm correlation. Methods to stack individual SDMs include summing individual probabilities and thresholding then summing. Thresholding can be based on a specific evaluation metric or by drawing repeatedly from a Bernouilli distribution.'),
                h3(' '),
                p('This interface allows you to load new data or previously built SDMs, ensemble SDMs or SSDMs.'),
                h3(' '),
                p('Please refer to the package help (?SSDM) for further information on available functionalities.')
                ),

        ### Loading page ###

        ## Loading new data ##
        tabItem('newdata',
                fluidPage(
                  fluidRow(
                    box(title = 'Environmental variable', height = 600,
                        uiOutput('Envbug'),
                        shinyFilesButton('envfiles', 'Rasters selection', 'Please select rasters', FALSE, multiple = T),
                        uiOutput('factors'),
                        p('Which variable should be considered as a categorical variable'),
                        checkboxGroupInput('load.var.options', 'loading options', list('Normalization'), selected = 'Normalization', inline = T),
                        actionButton('load', 'Load')
                    ),
                    box(title = 'Preview', height = 600,
                        uiOutput('layerchoice'),
                        plotOutput('env'))
                  ),
                  fluidRow(
                    box(title = 'Occurrence table',
                        uiOutput('Occbug'),
                        shinyFilesButton('Occ', 'Occurrence selection', 'Please select occurrence file', FALSE),
                        radioButtons('sep', 'Separator',
                                     c(Comma = ',',
                                       Semicolon = ';',
                                       Tab = '\t',
                                       'White space' = ' '),
                                     ',', inline = T),
                        radioButtons('dec', 'Decimal',
                                     c(Point ='.',
                                       Comma = ','),
                                     '.', inline = T),
                        uiOutput('Xcol'),
                        uiOutput('Ycol'),
                        uiOutput('Spcol'),
                        uiOutput('Pcol'),
                        checkboxInput('GeoRes', 'Geographic resampling', value = T),
                        uiOutput('reso'),
                        p('Geographical thinning can be performed on occurrences to limit geographical biases in the SDMs.'),
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
                    box(title = 'Load SDM',
                        shinyDirButton('prevmodel', 'Select previous model folder', 'Please select previous model folder', FALSE),
                        selectInput('model.type','Choose model type', c(' ', 'Ensemble SDM', 'SSDM')),
                        actionButton('load.model', 'load', icon = icon('file'))
                    ),
                    box('Preview',
                        textOutput('prevmodelbug'),
                        plotOutput('model.preview')
                    )
                  )
                )
        ),

        ### Modelling Page ####
        tabItem('modelling',
                fluidRow(
                  tabBox(
                    tabPanel('Basic',
                             uiOutput('specie'),
                             uiOutput('algoUI'),
                             uiOutput('repUI'),
                             uiOutput('nameUI'),
                             uiOutput('uncertUI'),
                             uiOutput('uncertinfoUI'),
                             uiOutput('endemismUI'),
                             uiOutput('endemisminfoUI'),
                             uiOutput('metricUI'),
                             uiOutput('metricinfoUI'),
                             uiOutput('methodUI'),
                             uiOutput('methodinfoUI'),
                             uiOutput('repBslide')
                    ),
                    tabPanel('Intermediate',
                             checkboxInput('PA', 'Automatic Pseudo-Absences', value = T),
                             uiOutput('PAnbUI'),
                             uiOutput('PAstratUI'),
                             selectInput('cval', 'Model evaluation method', c('holdout','k-fold','LOO'), selected = 'holdout'),
                             uiOutput('cvalinfoUI'),
                             uiOutput('cvalparam1UI'),
                             uiOutput('cvalparam2UI'),
                             selectInput('axesmetric', 'Variable importance evaluation metric', c('Pearson','AUC','Kappa','sensitivity','specificity','prop.correct'), selected = 'Pearson'),
                             uiOutput('axesmetricinfoUI'),
                             uiOutput('weightUI'),
                             uiOutput('ensemblemetricUI'),
                             uiOutput('AUCUI'),
                             uiOutput('KappaUI'),
                             uiOutput('sensitivityUI'),
                             uiOutput('specificityUI'),
                             uiOutput('propcorrectUI'),
                             uiOutput('ensembleinfoUI'),
                             uiOutput('rangeUI'),
                             uiOutput('rangevalUI'),
                             uiOutput('rangeinfoUI')
                    ),
                    tabPanel('Advanced',
                             sliderInput('thresh','Threshold precision',11,10001,101, step = 10),
                             p('Value representing the number of equal interval threshold values between 0 and 1'),
                             uiOutput('tmpUI'),
                             uiOutput('tmpinfoUI'),
                             checkboxGroupInput('algoparam', 'Algorithm parameters', c('GLM','GAM','MARS','GBM','CTA','RF','ANN','SVM'), inline = T),
                             uiOutput('testUI'),
                             uiOutput('epsilonUI'),
                             uiOutput('maxitUI'),
                             uiOutput('threshshrinkUI'),
                             uiOutput('treesUI'),
                             uiOutput('finalleaveUI'),
                             uiOutput('cvUI'),
                             uiOutput('degreeUI')
                    )
                  ),
                  box(title = 'Preview', textOutput('modelfailed'), plotOutput('modelprev'))
                ),
                fluidRow(
                  box(
                    actionButton('model','Model')
                  )
                )
                ),

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
                                   title = 'Local species richness'),
                         tabPanel(plotOutput('endemism'), title = 'Endemism map', textOutput('endemism.info')),
                         tabPanel(plotOutput('Uncertainity'), title = 'Uncertainty'),
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
                uiOutput('algoevalcorr')
        ),

        ## Save model ##
        tabItem('save',
                fluidPage(
                  fluidRow(
                    box(title = 'Save results',
                        textOutput('savegetwd'),
                        shinyDirButton('save', 'Folder selection', 'Please select folder to save the model', FALSE),
                        actionButton('savemodel', 'save', icon = icon('floppy-o'))
                    )
                  )
                )
        )
      )
    )
  )

  #### Server ####
  server <- function(input, output, session) {
    ### Server data ###
    data <- reactiveValues(Env = stack(), Occ = data.frame(), dir = getwd(), ENM = NULL, enms = list(), Stack = NULL)
    result <- reactiveValues(ENM = NULL)

    ### Menu ###
    output$menu <- renderMenu({
      sidebarMenu(
        id = 'actions',
        menuItem('Welcome page', tabName = 'welcomepage', selected = T),
        menuItem('Load',
                 menuSubItem('New data', tabName = 'newdata'),
                 menuSubItem('Previous model', tabName = 'previousmodel')
        ),
        if(length(data$Occ) > 0) {
          menuItem('Modelling',
                   menuSubItem('Modelling tab', 'modelling'),
                   selectInput('modellingchoice', 'Modelling type ',
                               c('Algorithm modelling', 'Ensemble modelling', 'Stack modelling'),
                               selected = 'Stack modelling')
          )
        },
        if(!is.null(data$Stack) || !is.null(data$ENM)) {
          if(inherits(data$ENM,'Algorithm.SDM')) {tabname = 'Algorithm SDM'} else {tabname = 'Ensemble SDM'}
          menuItem('Results',
                   if(!is.null(data$Stack)){menuItem("SSDM", tabName = "stack", icon = icon("dashboard"))},
                   menuItem(tabname, tabName = "stackenm", icon = icon("pagelines")),
                   if(!is.null(data$Stack)){selectInput('enmchoice', 'Species:', data$enms, selectize = TRUE)},
                   if(!inherits(data$ENM,'Algorithm.SDM')) {menuItem('Save', tabName = "save", icon = icon("floppy-o"))}
          )
        }
      )
    })

    ### Load Menu ###

    ## Load new data page ##

    # Environmental variable loading
    load.var <- reactiveValues(factors = c(), formats = c(), norm = T,  vars = list())
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyFileChoose(input, 'envfiles', session=session,
                     roots=c( wd='.', home = '/home', root = '/'),
                     filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = T)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyFileChoose(input, 'envfiles', session=session, roots=c( wd='.', disks),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    } else {
      shinyFileChoose(input, 'envfiles', session=session, roots=c( wd='.', home = '/user', root = '/'),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    }
    observeEvent(input$envfiles,{
      load.var$vars = list()
      for(i in 1:length(input$envfiles$files)){
        load.var$vars[[i]] = as.character(input$envfiles$files[[i]][length(input$envfiles$files[[i]])])
      }
    })
    output$factors <- renderUI({
      selectInput('factors', 'Categorical', load.var$vars, multiple = T, selectize = T)
    })
    observeEvent(input$load, {
      validate(
        need(length(load.var$vars) > 0, 'Choose environment variable files first !')
      )
      if(Sys.info()[['sysname']] == 'Linux') {
        path = switch(input$envfiles$root,
                      'wd' = getwd(),
                      'home' = '/home',
                      'root' = '/')
      } else if (Sys.info()[['sysname']] == 'Windows') {
        path = switch(input$envfiles$root, 'wd' = getwd())
      } else {
        path = switch(input$envfiles$root,
                      'wd' = getwd(),
                      'home' = '/home',
                      'root' = '/')
      }
      for(i in 2:(length(input$envfiles$files[[1]]))-1){
        path = paste0(path, '/', input$envfiles$files[[1]][i])
      }
      load.var$formats = c()
      for (i in 1 :length(load.var$vars)) {
        format = paste0('.',strsplit(load.var$vars[[i]], '.', fixed = T)[[1]][2])
        if (!(format %in% load.var$formats)) {load.var$formats = c(load.var$formats, format)}
      }
      if('Normalization' %in% input$load.var.options) {
        load.var$norm = T
      } else {
        load.var$norm = F
      }
      a = try(withProgress(message = 'Variables loading',
                              load_var(path,
                                       files = unlist(load.var$vars),
                                       format = load.var$formats,
                                       Norm = load.var$norm,
                                       tmp = F,
                                       factors = load.var$factors,
                                       verbose = F,
                                       GUI = T)))
      if(inherits(a, 'try-error')){
        output$Envbug <- renderUI(p('Environmental variables loading failed, please check your inputs and try again'))
      } else {
        output$Envbug <- renderUI(p())
        data$Env = a
        for (i in 1 :length(load.var$vars)) {
          names(data$Env)[i] = strsplit(load.var$vars[[i]], '.', fixed = T)[[1]][1]
        }
        output$layerchoice <- renderUI({
          selectInput('layer', 'Variable', as.list(names(data$Env)), multiple = F, selectize = T)
        })
        output$env <- renderPlot({
          if(!is.null(input$layer)){
            i = as.numeric(which(as.list(names(data$Env)) == input$layer))
            if(data$Env[[i]]@data@isfactor) {
              map = !as.factor(data$Env[[i]])
            } else {
              map = data$Env[[i]]
            }
            spplot(map,
                   main = names(data$Env)[i],
                   xlab = 'Longitude (\u02DA)',
                   ylab = 'Latitude (\u02DA)',
                   col.regions = rev(terrain.colors(10000)))
          }
        })
      }
      })

    # Occurrences loading
    load.occ <- reactiveValues(columns = c())
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyFileChoose(input, 'Occ', session=session,
                      roots=c( wd='.', home = '/home', root = '/'),
                      filetypes=c('',"csv", "txt"))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = T)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyFileChoose(input, 'Occ', session=session, roots=c( wd='.', disks),
                      filetypes=c('',"csv", "txt"))
    } else {
      shinyFileChoose(input, 'Occ', session=session, roots=c( wd='.', home = '/user', root = '/'),
                      filetypes=c('',"csv", "txt"))
    }
    observeEvent(input$Occ, {
      file = paste0(switch(input$Occ$root,
               'wd' = getwd(),
               'home' = '/home',
               'root' = '/'), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
      load.occ$columns = names(read.csv2(file))
    })
    observeEvent(input$sep, {
      if(!is.null(input$Occ)) {
        file = paste0(switch(input$Occ$root,
                             'wd' = getwd(),
                             'home' = '/home',
                             'root' = '/'), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
        load.occ$columns = names(read.csv2(file, sep = input$sep, nrows = 0))
        }
    })
    observeEvent(input$Occ, {
      if(!is.null(input$Occ)) {
        file = paste0(switch(input$Occ$root,
                             'wd' = getwd(),
                             'home' = '/home',
                             'root' = '/'), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
        load.occ$columns = names(read.csv2(file, sep = input$sep, nrows = 0))
        }
    })
    output$Xcol <- renderUI({selectInput('Xcol', 'X column', load.occ$columns, multiple = F)})
    output$Ycol <- renderUI({selectInput('Ycol', 'Y column', load.occ$columns, multiple = F)})
    output$Pcol <- renderUI({selectInput('Pcol', 'Presence (1) / absence (0) column', c('None', load.occ$columns), multiple = F)})
    output$Spcol <- renderUI({selectInput('Spcol', 'Species column', c('None', load.occ$columns), multiple = F)})
    output$reso <- renderUI({if(input$GeoRes){sliderInput('reso', 'Resampling grid coefficient', 1,10,1)}})
    observeEvent(input$load2, {
      validate(
        need(length(data$Env@layers) > 0, 'You need to load environmental variable before !'),
        need(length(input$Occ) > 0, 'Choose occurrences file first !')
      )
      if(is.null(input$dec)) {dec = ","} else {dec = input$dec}
      if (input$Spcol == 'None') {Spcol = NULL} else {Spcol = input$Spcol}
      if(is.null(input$sep)) {sep = ""} else {sep = input$sep}
      file = paste0(switch(input$Occ$root,
                           'wd' = getwd(),
                           'home' = '/home',
                           'root' = '/'), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
      a = try(withProgress(message = 'Occurrences loading',
                              load_occ(path = {},
                                       data$Env,
                                       file,
                                       Xcol = input$Xcol,
                                       Ycol = input$Ycol,
                                       Spcol = Spcol,
                                       GeoRes = input$GeoRes,
                                       reso = max(res(data$Env@layers[[1]])) * as.numeric(input$reso),
                                       verbose = F,
                                       GUI = T,
                                       sep = sep,
                                       dec = dec)))
      if(inherits(a, 'try-error')){
        output$Occbug <- renderUI(p('Occurrences loading failed, please check your inputs and try again'))
      } else {
        output$Occbug <- renderUI(p(' '))
        data$Occ = a
      }
    })
    output$occ <- renderDataTable({if(length(data$Occ) > 0) {data$Occ}})

    ## Load previous model page ##
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyDirChoose(input, 'prevmodel', session=session, roots=c( wd='.', home = '/home', root = '/'), filetypes=c(''))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = T)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyDirChoose(input, 'prevmodel', session=session, roots=c( wd='.', disks), filetypes=c(''))
    } else {
      shinyDirChoose(input, 'prevmodel', session=session, roots=c( wd='.', home = '/user', root = '/'), filetypes=c(''))
    }
    observeEvent(input$load.model, {
      validate(
        need(input$model.type != ' ', 'You need to choose the model type first !')
      )
      if(Sys.info()[['sysname']] == 'Linux') {
        path = switch(input$prevmodel$root,
                      'wd' = getwd(),
                      'home' = '/home/',
                      'root' = '/')
      } else if (Sys.info()[['sysname']] == 'Windows') {
        path = switch(input$prevmodel$root, 'wd' = getwd())
      } else {
        path = switch(input$prevmodel$root,
                      'wd' = getwd(),
                      'home' = '/home',
                      'root' = '/')
      }
      for(i in 2:(length(input$prevmodel$path)-1)){
        path = paste0(path, '/', input$prevmodel$path[[i]][1])
      }
      name = input$prevmodel$path[[length(input$prevmodel$path)]][1]
      if (input$model.type == 'Ensemble SDM') {
        a = try(load_enm(name, path))
      }
      if (input$model.type == 'SSDM') {
        a = try(withProgress(message = 'Model loading', load_stack(name, path, GUI = T)))
      }
      if(inherits(a, 'try-error')){
        output$prevmodelbug <- renderText('Previous model loading failed, please check your inputs and try again')
      } else {
        output$prevmodelbug <- renderText(' ')
        if (input$model.type == 'Ensemble SDM') {
          data$ENM = a
          if(!is.null(data$ENM)){result$ENM = data$ENM}
          output$model.preview <- renderPlot({spplot(data$ENM@projection,
                                                     main = 'Habitat suitability map',
                                                     xlab = 'Longitude (\u02DA)',
                                                     ylab = 'Latitude (\u02DA)',
                                                     col.regions = rev(terrain.colors(10000)))})
        }
        if (input$model.type == 'SSDM') {
          data$Stack = a
          for (i in 1:length(names(data$Stack@enms))){data$enms[[i]] = strsplit(names(data$Stack@enms), '.', fixed = T)[[i]][1]}
          output$model.preview <- renderPlot({spplot(data$Stack@diversity.map,
                                                     main = data$Stack@name,
                                                     xlab = 'Longitude (\u02DA)',
                                                     ylab = 'Latitude (\u02DA)',
                                                     col.regions = rev(terrain.colors(10000)))})
        }
      }
    })

    ## Basic modelling parameters ##
    output$specie <-renderUI({
      if(input$modellingchoice != 'Stack modelling'){
        if (input$Spcol != 'None') {
          selectInput('specie', 'Specie', as.list(levels(data$Occ[,which(names(data$Occ) == input$Spcol)])))
        }
      }
    })
    output$algoUI <- renderUI({
      if(input$modellingchoice == 'Algorithm modelling'){
        title = 'Algorithm'
        multiple = F
      } else {
        title = 'Algorithms'
        multiple = T
      }
      selectInput('algo', title, c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM'), multiple = multiple, selectize = T)})
    output$repUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){sliderInput('rep','Repetitions',1,50,10, step = 1)}})
    output$nameUI <- renderUI({textInput('name','Name')})
    output$uncertUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){checkboxInput('uncert', 'Uncertainty mapping', value = T)}})
    output$uncertinfoUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){p('Between-algorithm variance map')}})
    output$endemismUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){selectInput('endemism', 'Endemism mapping', c('None', 'WEI', 'CWEI'), selected = 'WEI')}})
    output$endemisminfoUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){
      p(switch(input$endemism,
               'None' = 'No endemism map will be built',
               'WEI' = 'Endemism map will be built by counting all species in each cell and weighting each by the inverse of its number of occurrences',
               'CWEI' = 'Endemism map will be built by dividing the weighted endemism index by the total count of species in the cell'))}})
    output$metricUI <- renderUI({selectInput('metric', 'Evaluation metric', c('Kappa','CCR','TSS','SES','LW','ROC'), selected = 'SES')})
    output$metricinfoUI <- renderUI({
      p(switch(input$metric,
               'Kappa' = 'Maximizes the Kappa',
               'CCR' = 'Maximizes the sum of sensitivity and specificity',
               'TSS' = '(True Skill Statistic) maximizes the sum of sensitivity and specificity',
               'SES' = 'Uses the sensitivity-specificity equality',
               'LW' = 'Uses the lowest occurrence prediction probability',
               'ROC' = 'Minimizes the distance between the ROC plot (receiving operating curve) and the upper left corner (1,1).'))
    })
    output$methodUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){selectInput('method', 'Diversity mapping method', c('Probability','Random Bernoulli','Threshold'), selected = 'Probability')}})
    output$methodinfoUI <- renderUI({
      if(input$modellingchoice == 'Stack modelling'){
      p(switch(input$method,
               'Probability' = 'Sum probabilities of habitat suitability maps',
               'Random Bernoulli' = 'Drawing repeatedly from a Bernoulli distribution',
               'Threshold' = 'Sum the binary map obtained with the thresholding (depending on the metric, see metric parameter)'))
      }
    })
    output$repBslide <- renderUI({if(input$modellingchoice == 'Stack modelling'){if(!is.null(input$method)){if(input$method=='Random Bernoulli'){sliderInput('repB','Bernoulli repetitions',1,10000,1000, step = 1)}}}})

    ## Intermediate modelling parameters ##
    output$PAnbUI <- renderUI(if(!input$PA){sliderInput('PAnb','PA number',1,10000,1000, step = 100)})
    output$PAstratUI <- renderUI(if(!input$PA){selectInput('PAstrat', 'PA strategy', c('random','disk'), selected = 'random')})
    output$cvalinfoUI <- renderUI({
        p(switch(input$cval,
                 'holdout' = 'Data are partitioned into a training set and an evaluation set using a fraction and the operation can be repeated',
                 'k-fold' = 'Data are partitioned into k folds being k-1 times in the training set and once the evaluation set and the operation can be repeated',
                 'LOO' = '(Leave One Out) each point is successively taken as evaluation data'))
    })
    output$cvalparam1UI <- renderUI({
      if(input$cval == 'holdout'){
        title = 'Holdout fraction'
        min = 0
        max = 1
        val = 0.7
        step = 0.05
      }
      if(input$cval == 'k-fold'){
        title = 'K fold number'
        min = 1
        max = 50
        val = 10
        step = 1
      }
      if(input$cval != 'LOO'){sliderInput('cvalparam1', title, min, max, val, step)}
    })
    output$cvalparam2UI <- renderUI({if(input$cval != 'LOO'){sliderInput('cvalrep', 'Cross-validation repetitions', 1,20,1, step = 1)}})
    output$axesmetricinfoUI <- renderUI({
      p(paste('Metric used to evaluate the variable relative importance (difference between a full model and one with each variable successively omitted):',
              switch(input$axesmetric,
               'Pearson' = 'Computes a simple Pearson\'s correlation r between predictions of the full model and the one without a variable, and returns the score 1-r: the highest the value, the more influence the variable has on the model',
               'AUC' = 'Area under the receiving operating characteristic curve',
               'Kappa' = 'Kappa',
               'sensitivity' = 'Sensitivity',
               'specificity' = 'Specificity',
               'prop.correct' = 'Proportion of correctly predicted occurrences')))
    })
    output$ensemblemetricUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){selectInput('ensemblemetric', 'Ensemble selection metric', c('AUC','Kappa','sensitivity','specificity','prop.correct'), selected = 'AUC', multiple = T, selectize = T)}})
    output$ensembleinfoUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){p('Metric(s) and threshold(s) used to select the best SDMs that will be kept in the ensemble SDM')}})
    output$AUCUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('AUC' %in% input$ensemblemetric){sliderInput('AUC','AUC threshold',0.5,1,0.75, step = 0.05)}}})
    output$KappaUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('Kappa' %in% input$ensemblemetric){sliderInput('Kappa','Kappa threshold',0,1,0.5, step = 0.05)}}})
    output$sensitivityUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('sensitivity' %in% input$ensemblemetric){sliderInput('sensitivity','Sensitivity threshold',0.5,1,0.75, step = 0.05)}}})
    output$specificityUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('specificity' %in% input$ensemblemetric){sliderInput('specificity','Specificity threshold',0.5,1,0.75, step = 0.05)}}})
    output$propcorrectUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('prop.correct' %in% input$ensemblemetric){sliderInput('propcorrect','Correct proportion threshold',0.5,1,0.75, step = 0.05)}}})
    output$weightUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){checkboxInput('weight', 'Ensemble weighted', value = T)}})
    output$rangeUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){checkboxInput('range', 'Range restriction', value = F)}})
    output$rangevalUI <- renderUI({if(input$modellingchoice == 'Stack modelling' && input$range){sliderInput('rangeval', 'Range restriction value', 1,100,5, step = 1)}})
    output$rangeinfoUI <- renderUI({if(input$modellingchoice == 'Stack modelling' && input$range){p('Set a value of range restriction (in pixels) around presence occurrences on habitat suitability maps (all further points will have a null probability)')}})

    ## Advanced modelling parameters ##
    output$tmpUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){checkboxInput('tmp', 'Temporary files', value = F)}})
    output$tmpinfoUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){p('Rasters are saved in temporary files to release memory')}})
    output$testUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam){selectInput('test','Test (GLM/GAM)',c('AIC', 'BIC'))}})
    output$epsilonUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'SVM' %in% input$algoparam){sliderInput('epsilon','Epsilon (GLM/GAM/SVM) = 1e-X',1,20,8, step = 1)}})
    output$maxitUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'ANN' %in% input$algoparam){sliderInput('maxit','Maximum number of iteration (GLM/GAM/ANN)',100,1000,500, step = 100)}})
    output$degreeUI <- renderUI({if('MARS' %in% input$algoparam){sliderInput('degree','Degree of interaction (MARS)',1,10,2, step = 1)}})
    output$threshshrinkUI <- renderUI({if('GBM' %in% input$algoparam){sliderInput('threshshrink','Shrinkage threshold (GBM)',1e-05,1,1e-03, step = 1e-04)}})
    output$treesUI <- renderUI({if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam){sliderInput('trees','Maximum number of trees (GBM/RF)',500,10000,2500, step = 500)}})
    output$finalleaveUI <- renderUI({if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam || 'CTA' %in% input$algoparam){sliderInput('finalleave','Final leaves size (GBM/CTA/RF)',1,20,1, step = 1)}})
    output$cvUI <- renderUI({if('GBM' %in% input$algoparam || 'CTA' %in% input$algoparam || 'SVM' %in% input$algoparam){sliderInput('cv','Cross-validation number (GBM/CTA/SVM)',1,20,3, step = 1)}})

    ## Modelling ##
    observeEvent(input$model,{
      validate(
        need(length(input$algo) > 0, 'Choose algorithm(s) to be run !')
      )
      data$ENM = NULL
      data$Stack = NULL
      data$enms = list()
      result$ENM = NULL
      output$modelprev = renderPlot({plot.new()})
      output$modelfailed = renderText({''})
      if (input$Pcol == 'None') {Pcol = NULL} else {Pcol = input$Pcol}
      if (input$Spcol == 'None') {Spcol = NULL} else {Spcol = input$Spcol}
      if (input$name == "") {name = NULL} else {name = input$name}
      algo = c()
      for (i in 1:length(input$algo)){algo = c(algo, input$algo[i])}
      if(!is.null(input$PA)) {PA = T}
      if(input$PA) {PA = NULL} else {PA = list('nb' = input$PAnb, 'strat' = input$PAstrat)}
      if(is.null(input$ensemblemetric)){
        ensemble.metric = c('AUC')
        ensemble.thresh = c(0.75)
      } else {
        ensemble.metric = c()
        ensemblethresh = c(input$AUC, input$Kappa, input$sensitivity, input$specificity, input$propcorrect)
        ensemble.thresh = c()
        for (i in 1:length(input$ensemblemetric)){
          ensemble.metric = c(ensemble.metric, input$ensemblemetric[i])
          ensemble.thresh = c(ensemble.thresh, as.numeric(ensemblethresh[i]))
        }
      }
      if(is.null(input$range) || !input$range) {range = NULL} else {range = input$rangeval}
      if(is.null(input$endemism) || input$endemism == 'None') {endemism = NULL} else {endemism = input$endemism}
      algoparam = list()
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam){algoparam$test = input$test} else {algoparam$test = 'AIC'}
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'SVM' %in% input$algoparam){algoparam$epsilon = as.numeric(paste0('1e-', as.character(input$epsilon)))} else {algoparam$epsilon = 1e-08}
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'ANN' %in% input$algoparam){algoparam$maxit = as.numeric(input$maxit)} else {algoparam$maxit = 500}
      if('MARS' %in% input$algoparam){algoparam$degree = as.numeric(input$degree)} else {algoparam$degree = 3}
      if('GBM' %in% input$algoparam){algoparam$thresh.shrink = as.numeric(input$threshshrink)} else {algoparam$threshshrink = 1e-03}
      if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam){algoparam$trees = as.numeric(input$trees)} else {algoparam$trees = 2500}
      if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam || 'CTA' %in% input$algoparam){algoparam$final.leave = as.numeric(input$finalleave)}  else {algoparam$finalleave = 1}
      if('GBM' %in% input$algoparam || 'CTA' %in% input$algoparam || 'SVM' %in% input$algoparam){algoparam$cv = as.numeric(input$cv)}  else {algoparam$cv = 3}
      method = switch(input$method,
                      'Probability' = 'P',
                      'Random Bernoulli' = 'B',
                      'Threshold' = 'T')
      if(is.null(input$repB)) {rep.B = 1000} else {rep.B = as.numeric(input$repB)}
      if(is.null(input$cval)) {cval = 'holdout'} else {cval = input$cval}
      if(length(c(as.numeric(input$cvalparam1), as.numeric(input$cvalrep))) < 2) {cv.param = c(0.7, 1)} else {cv.param = c(as.numeric(input$cvalparam1), as.numeric(input$cvalrep))}
      if(!inherits(input$tmp,'logical')) {tmp = F} else {tmp = input$tmp}
      if(!inherits(input$weight,'logical')) {weight = T} else {weight = input$weight}
      if(input$modellingchoice == 'Algorithm modelling'){
        if (input$Spcol != 'None') {
          Occ = data$Occ[which(data$Occ[,which(names(data$Occ) == input$Spcol)] == input$specie),]
        } else {
          Occ = data$Occ
        }
        data$ENM = withProgress(message = 'SDM',
                                modelling(algo,
                                          Occ, data$Env,
                                          Xcol = input$Xcol,
                                          Ycol = input$Ycol,
                                          Pcol = Pcol,
                                          name = name,
                                          save = F,
                                          path = 'nowhere',
                                          PA = PA,
                                          cv = cval,
                                          cv.param = cv.param,
                                          thresh = as.numeric(input$thresh),
                                          axes.metric = input$axesmetric,
                                          select.metric = c('AUC'),
                                          select.thresh = c(0),
                                          select = F,
                                          metric = input$metric,
                                          verbose = F,
                                          GUI = T,
                                          test = algoparam$test,
                                          maxit = algoparam$maxit,
                                          epsilon = algoparam$epsilon,
                                          degree = algoparam$degree,
                                          threshshrink = algoparam$threshshrink,
                                          trees = algoparam$trees,
                                          finalleave = algoparam$finalleave,
                                          algocv = algoparam$cv))
        output$modelprev = renderPlot(spplot(data$ENM@projection,
                                             main = 'Habitat suitability map',
                                             xlab = 'Longitude (\u02DA)',
                                             ylab = 'Latitude (\u02DA)',
                                             col.regions = rev(terrain.colors(10000))))
        result$ENM = data$ENM
      }
      if(!inherits(input$uncert,'logical')) {uncert = T} else {uncert = input$uncert}
      if(input$modellingchoice == 'Ensemble modelling'){
        if (input$Spcol != 'None') {
          Occ = data$Occ[which(data$Occ[,which(names(data$Occ) == input$Spcol)] == input$specie),]
        } else {
          Occ = data$Occ
        }
        data$ENM = withProgress(message = 'Ensemble SDM',
                                  ensemble_modelling(algo,
                                                  Occ, data$Env,
                                                  Xcol = input$Xcol,
                                                  Ycol = input$Ycol,
                                                  Pcol = Pcol,
                                                  rep = as.numeric(input$rep),
                                                  name = name,
                                                  save = F,
                                                  path = 'nowhere',
                                                  PA = PA,
                                                  cv = cval,
                                                  cv.param = cv.param,
                                                  thresh = as.numeric(input$thresh),
                                                  axes.metric = input$axesmetric,
                                                  uncertainty =  uncert,
                                                  tmp =  tmp,
                                                  ensemble.metric = ensemble.metric,
                                                  ensemble.thresh = ensemble.thresh,
                                                  weight = weight,
                                                  metric = input$metric,
                                                  verbose = F,
                                                  GUI = T,
                                                  test = algoparam$test,
                                                  maxit = algoparam$maxit,
                                                  epsilon = algoparam$epsilon,
                                                  degree = algoparam$degree,
                                                  threshshrink = algoparam$threshshrink,
                                                  trees = algoparam$trees,
                                                  finalleave = algoparam$finalleave,
                                                  algocv = algoparam$cv))
        if(!is.null(data$ENM)){
        output$modelprev = renderPlot(spplot(data$ENM@projection,
                                             main = 'Habitat suitability map',
                                             xlab = 'Longitude (\u02DA)',
                                             ylab = 'Latitude (\u02DA)',
                                             col.regions = rev(terrain.colors(10000))))
        result$ENM = data$ENM
        } else {
          output$modelfailed = renderText('No ensemble SDM were kept, maybe you should try lower ensemble threshold(s) ?')
        }
      }
      if(input$modellingchoice == 'Stack modelling'){
        data$Stack = withProgress(message = 'SSDM',
                                  stack_modelling(algo,
                                                  data$Occ, data$Env,
                                                  Xcol = input$Xcol,
                                                  Ycol = input$Ycol,
                                                  Pcol = Pcol,
                                                  Spcol = Spcol,
                                                  rep = as.numeric(input$rep),
                                                  name = name,
                                                  save = F,
                                                  path = 'nowhere',
                                                  PA = PA,
                                                  cv = cval,
                                                  cv.param = cv.param,
                                                  thresh = as.numeric(input$thresh),
                                                  axes.metric = input$axesmetric,
                                                  uncertainty =  uncert,
                                                  tmp =  tmp,
                                                  ensemble.metric = ensemble.metric,
                                                  ensemble.thresh = ensemble.thresh,
                                                  weight = weight,
                                                  method = method,
                                                  metric = input$metric,
                                                  rep.B =  rep.B,
                                                  range = range,
                                                  endemism = endemism,
                                                  verbose = T,
                                                  GUI = T,
                                                  test = algoparam$test,
                                                  maxit = algoparam$maxit,
                                                  epsilon = algoparam$epsilon,
                                                  degree = algoparam$degree,
                                                  threshshrink = algoparam$threshshrink,
                                                  trees = algoparam$trees,
                                                  finalleave = algoparam$finalleave,
                                                  algocv = algoparam$cv))
        if(!is.null(data$Stack)){
        output$modelprev = renderPlot(spplot(data$Stack@diversity.map,
                                             main = 'Local species richness',
                                             xlab = 'Longitude (\u02DA)',
                                             ylab = 'Latitude (\u02DA)',
                                             col.regions = rev(terrain.colors(10000))))
        for (i in 1:length(names(data$Stack@enms))){data$enms[[i]] = strsplit(names(data$Stack@enms), '.', fixed = T)[[i]][1]}
        } else {
          output$modelfailed = renderText('You have less than two remaining ensemble SDMs, maybe you should try lower ensemble threshold(s) ?')
        }
      }
      })

    ## Stack results ##
    output$Stackname <- renderUI(h2(data$Stack@name, align = 'center'))
    ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        if(!is.null(ranges$x)){ref = crop(data$Stack@diversity.map, c(ranges$x, ranges$y))} else {ref = data$Stack@diversity.map}
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
    output$Diversity <- renderPlot({
      eval = 'Mean'
      ensemble.metric = strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1]
      for (i in 1:length(ensemble.metric)) {
        eval = paste(eval, paste(ensemble.metric[i],':',round(data$Stack@evaluation[1,which(names(data$Stack@evaluation) == ensemble.metric[i])], digits = 3)))
        if (i < length(ensemble.metric)) {eval = paste(eval, ',')}
      }
      if (!is.null(ranges$x)) {diversity = crop(data$Stack@diversity.map, c(ranges$x, ranges$y))} else {diversity = data$Stack@diversity.map}
      spplot(diversity,
           main = eval,
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)))
    })
    output$endemism <- renderPlot({
      if (!is.null(ranges$x)) {endemism = crop(data$Stack@endemsim.map, c(ranges$x, ranges$y))} else {endemism = data$Stack@endemism.map}
      spplot(endemism,
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)',
             col.regions = rev(terrain.colors(10000)))
    })
    output$Uncertainity <- renderPlot({
      if (!is.null(ranges$x)) {uncert = crop(data$Stack@uncertainty, c(ranges$x, ranges$y))} else {uncert = data$Stack@uncertainty}
      spplot(uncert,
             xlab = 'Longitude (\u02DA)',
             ylab = 'Latitude (\u02DA)',
             col.regions = rev(terrain.colors(10000)))
      })
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
                trace = "none", key = FALSE, margins = c(7, 11), na.rm = T,
                col = rev(heat.colors(1000)))
    })
    # Variable importance
    output$varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(data$Stack@variable.importance))
      names(varimp) = 'Axes.evaluation'
      bar = barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)),
                    ylim = c(0,(max(varimp$Axes.evaluation)+max((varimp[2])))),
                    las = 2, ylab = 'Variable relative contribution (%)')
      arrows(bar,varimp$Axes.evaluation+varimp[,2], bar, varimp$Axes.evaluation-varimp[,2], angle=90, code=3, length=0.1)
    })
    output$varimp.table <- renderTable({data$Stack@variable.importance})
    output$varimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(data$Stack@variable.importance)), 'Variable' = names(data$Stack@variable.importance))})
    # Parameters
    output$summary <- renderTable({
      summary = data.frame(matrix(nrow = 7, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Final number of species', 'Original algorithms', 'Number of repetitions',
                             'Pseudo-absences selection', 'Cross-validation method', 'Cross-validation parameters')
      algo.info = character()
      for (i in 1:length(strsplit(data$Stack@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(data$Stack@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (data$Stack@parameters$PA) {PA = 'default'}
      if(data$Stack@parameters$cv == 'LOO') {cv.param = 'None'}
      if(data$Stack@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                         strsplit(data$Stack@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                         'rep =',
                                                         strsplit(data$Stack@parameters$cv.param, '|', fixed = T)[[1]][3])}
      if(data$Stack@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                        strsplit(data$Stack@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                        'rep =',
                                                        strsplit(data$Stack@parameters$cv.param, '|', fixed = T)[[1]][3])}
      summary$Summary = c(data$Stack@parameters$data, length(data$Stack@enms), algo.info, data$Stack@parameters$rep, PA, data$Stack@parameters$cv, cv.param)
      if(!is.null(data$Stack@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = data$Stack@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
    })
    output$diversity.info <- renderText({
      text = switch(data$Stack@parameters$method,
                    'P' = 'summing habitat suitability probabilities',
                    'T' = paste('summing habitat suitability binary maps after thresholding with',data$Stack@parameters$metric),
                    'B' = paste('drawing repeatdly from a Bernoulli distribution with',data$Stack@parameters$rep.B, 'repetitions'))
      text = paste('Local species richness map realized by', text)
      text
    })
    output$endemism.info <- renderText({
      text = switch(data$Stack@parameters$endemism,
                                     'Any' = 'Endemism map not built (unactivated)',
                                     'WEI' = 'Endemism map built with the Weighted Endemism Index (WEI)',
                                     'B' = 'Endemism map built with the Corrected Weighted Endemism Index (CWEI)')
      text
    })
    output$varimp.info <- renderText({
      if( data$Stack@parameters$axes.metric != 'Pearson') {
        varimp.info = paste('Variable relative contribution evaluated with', data$Stack@parameters$axes.metric)
      } else {
        varimp.info = paste('Variable relative contribution evaluated with', data$Stack@parameters$axes.metric, '\'s correlation coefficient')
      }
      varimp.info
    })
    output$evaluation.info <- renderText({
      evaluation.info = 'Evaluation of models with'
      for (i in 1:length(strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste0(evaluation.info, ' ', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],' (>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(data$Stack@parameters$axes.metric) && i != 1) {
          evaluation.info = paste0(evaluation.info, ' and ', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],' (>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste0(evaluation.info, ', ', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],' (>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      evaluation.info
    })

    ## ENM Stack Result ##
    observeEvent(input$enmchoice, {
      if(!is.null(data$Stack)){
        if(length(data$enms) > 0){result$ENM = data$Stack@enms[[which(data$enms == input$enmchoice)]]}
      }
    })
    observeEvent(input$proba_dblclick, {
      brush <- input$proba_brush
      if (!is.null(brush)) {
        if(!is.null(ranges$x)){ref = crop(result$ENM@projection, c(ranges$x, ranges$y))} else {ref = result$ENM@projection}
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
      if (!is.null(ranges$x)) {proba = crop(result$ENM@projection, c(ranges$x, ranges$y))} else {proba = result$ENM@projection}
      spplot(proba,
           main = paste('AUC :',round(result$ENM@evaluation$AUC,3),'  Kappa',round(result$ENM@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           col.regions = rev(terrain.colors(10000)),
           sp.layout=list(SpatialPoints(data.frame(X = result$ENM@data$X[which(result$ENM@data$Presence == 1)],
                                                   Y = result$ENM@data$Y[which(result$ENM@data$Presence == 1)])),
                          pch = 16, cex = 0.7, col = 'black'))
    })
    output$niche <- renderPlot({
      niche.map = reclassify(result$ENM@projection, c(-Inf,result$ENM@evaluation$threshold,0, result$ENM@evaluation$threshold,Inf,1))
      if (!is.null(ranges$x)) {niche.map = crop(niche.map, c(ranges$x, ranges$y))}
      spplot(niche.map,
             main = paste('AUC :',round(result$ENM@evaluation$AUC,3),'  Kappa',round(result$ENM@evaluation$Kappa,3)),
             col.regions = rev(terrain.colors(10000)))
      })
    output$enm.uncertainty <- renderPlot({
      if(!inherits(result$ENM,'Algorithm.SDM')) {
        if (!is.null(ranges$x)) {uncert.map = crop(result$ENM@uncertainty, c(ranges$x, ranges$y))} else {uncert.map = result$ENM@uncertainty}
        spplot(uncert.map, col.regions = rev(terrain.colors(10000)))
      }
    })
    # Evaluation
    output$enm.evaluation.barplot <- renderPlot({
      evaluation = result$ENM@algorithm.evaluation
      for (i in 1:length(row.names(result$ENM@algorithm.evaluation))) {row.names(evaluation)[i] = strsplit(as.character(row.names(result$ENM@algorithm.evaluation)[i]), '.SDM')[[1]][1]}
      for (i in 1:length(row.names(evaluation))) {row.names(evaluation)[i] = tail(strsplit(as.character(row.names(evaluation)[i]), '.', fixed = T)[[1]], n = 1)}
      evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
      table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
      barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
      legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
    })
    output$enm.evaluation.table <- renderTable({
      algo.eval = result$ENM@algorithm.evaluation
      for (i in 1:length(row.names(result$ENM@algorithm.evaluation))) {row.names(algo.eval)[i] = strsplit(as.character(row.names(result$ENM@algorithm.evaluation)[i]), '.SDM')[[1]][1]}
      for (i in 1:length(row.names(result$ENM@algorithm.evaluation))) {row.names(algo.eval)[i] = tail(strsplit(as.character(row.names(algo.eval)[i]), '.', fixed = T)[[1]], n = 1)}
      algo.eval[c(2,4:8)]
    })
    # Algorithms correlation
    output$enm.algo.corr.table <- renderTable({
      correlation = result$ENM@algorithm.correlation
      for (i in 1:length(row.names(result$ENM@algorithm.correlation))) {row.names(correlation)[i] = strsplit(as.character(row.names(result$ENM@algorithm.correlation)[i]), '.SDM')[[1]][1]}
      names(correlation) = row.names(correlation)
      if (length(correlation) > 0) {
        correlation[upper.tri(correlation, diag = T)] = NA
        correlation
      }
    })
    output$enm.algo.corr.heatmap <- renderPlot({
      correlation = result$ENM@algorithm.correlation
      for (i in 1:length(row.names(result$ENM@algorithm.correlation))) {row.names(correlation)[i] = strsplit(as.character(row.names(result$ENM@algorithm.correlation)[i]), '.SDM')[[1]][1]}
      names(correlation) = row.names(correlation)
      if (length(correlation) > 0) {
        correlation[upper.tri(correlation, diag = T)] = NA
        m <- as.matrix(correlation)
        heatmap.2(x = m, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
                  cellnote = round(m,3), notecol = "black", notecex = 2,
                  trace = "none", key = FALSE, margins = c(7, 11), na.rm = T,
                  col = rev(heat.colors(1000)))
      }
    })
    # Algo Eval Corr UI
    output$algoevalcorr <- renderUI({
      if(!inherits(result$ENM,'Algorithm.SDM')) {
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
      }
    })
    # Variable importance
    output$enm.varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(result$ENM@variable.importance))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$enm.varimp.table <- renderTable({result$ENM@variable.importance})
    output$enmvarimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(result$ENM@variable.importance)), 'Variable' = names(result$ENM@variable.importance))})
    # Parameters
    output$enm.summary <- renderTable({
      if(!inherits(result$ENM,'Algorithm.SDM')) {
      summary = data.frame(matrix(nrow = 7, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Type of occurrences', 'Number of occurrences', 'Originally selected algorithms', 'Number of repetitions',
                             'Pseudo-absence selection method', 'Cross-validation method', 'Cross-validation parameters')
      algo.info = character()
        for (i in 1:length(strsplit(result$ENM@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
          algo.info = paste(algo.info, strsplit(result$ENM@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
        }
      if (result$ENM@parameters$PA) {PA = 'default'}
      if (result$ENM@parameters$data == "presence-only data set") {
        nb.occ =  length(as.factor(result$ENM@data$Presence[which(result$ENM@data$Presence==1)])) / sum(result$ENM@algorithm.evaluation$kept.model)
      } else {
        nb.occ =  length(as.factor(result$ENM@data$Presence)) / sum(result$ENM@algorithm.evaluation$kept.model)
      }
      if(result$ENM@parameters$cv == 'LOO') {cv.param = 'None'}
      if(result$ENM@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                                  strsplit(result$ENM@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                                  'rep =',
                                                                  strsplit(result$ENM@parameters$cv.param, '|', fixed = T)[[1]][3])}
      if(result$ENM@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                                 strsplit(result$ENM@parameters$cv.param, '|', fixed = T)[[1]][2],
                                                                 'rep =',
                                                                 strsplit(result$ENM@parameters$cv.param, '|', fixed = T)[[1]][3])}
      summary$Summary = c(result$ENM@parameters$data, nb.occ, algo.info, result$ENM@parameters$rep, PA, result$ENM@parameters$cv, cv.param)
      if(!is.null(result$ENM@parameters$sp.nb.origin)) {
        summary = rbind(summary,
                        data.frame(Summary = result$ENM@parameters$sp.nb.origin, row.names = 'Original number of species'))
      }
      summary
      }
      else {
        data.frame('Not computed')
      }
    })
    output$enm.binary.info <- renderText({
      metric = switch(result$ENM@parameters$metric,
                                          'Kappa' = 'maximizing the Kappa',
                                          'CCR' = 'maximizing the proportion of correctly predicted occurrences (CCR)',
                                          'TSS' = 'maximizing the sensitivity and specificity sum (TSS)',
                                          'SES' = 'using the sensitivity and specificity equality',
                                          'LW' = 'using the lowest occurrence prediction probability',
                                          'ROC' = 'minimizing the distance between the ROC plot and the upper left corner')
      text = paste('Binary map realized by', metric,
                   'with a final threshold of',  round(result$ENM@evaluation$threshold, digits = 3))
      text
    })
    output$enm.varimp.info <- renderText({
      varimp.info = 'Variable relative contribution evaluated with '
      for (i in 1:length(result$ENM@parameters$axes.metric)) {
        if (i == 1) {
          varimp.info = paste(varimp.info, result$ENM@parameters$axes.metric[i])
        } else if (i == length(result$ENM@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', result$ENM@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', result$ENM@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$enm.evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in 1:length(strsplit(result$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1])) {
        if (i == 1) {
          evaluation.info = paste0(evaluation.info, ' ', strsplit(result$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],' (>',strsplit(result$ENM@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        } else if (i == length(result$ENM@parameters$axes.metric) && i != 1) {
          evaluation.info = paste0(evaluation.info, ' and ', strsplit(result$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],' (>',strsplit(result$ENM@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste0(evaluation.info, ' , ', strsplit(result$ENM@parameters$ensemble.metric, '.', fixed = T)[[1]][-1][i],' (>',strsplit(result$ENM@parameters$ensemble.thresh, '|', fixed = T)[[1]][-1][i],')')
        }
      }
      if (result$ENM@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })

    ### Save Menu ###

    ## Save model page ##
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyDirChoose(input, 'save', session=session, roots=c( wd='.', home = '/home', root = '/'), filetypes=c(''))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = T)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyDirChoose(input, 'save', session=session, roots=c( wd='.', disks), filetypes=c(''))
    } else {
      shinyDirChoose(input, 'save', session=session, roots=c( wd='.', home = '/user', root = '/'), filetypes=c(''))
    }
    observeEvent(input$savemodel, {
        path = switch(input$prevmodel$root,
                      'wd' = getwd(),
                      'home' = '/home',
                      'root' = '/')
      for(i in 2:length(input$prevmodel$path)){
        path = paste0(path, '/', input$prevmodel$path[[i]][1])
      }
      if(!is.null(data$ENM) && is.null(data$Stack)) {save.enm(data$ENM, name = data$ENM@name, path, verbose = F, GUI = T)}
        if(!is.null(data$Stack)) {save.stack(data$Stack, name = data$Stack@name, path, verbose = F, GUI = T)}
     })
  }


  #### Launching server ####
  #options(shiny.launch.browser = browser)
  options(shiny.maxRequestSize = maxmem)
  app <- shinyApp(ui, server)
  return(runApp(app))
}
