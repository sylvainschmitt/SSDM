#' @include load_occ.R load_var.R
#' @include Algorithm.SDM.R Ensemble.SDM.R Stacked.SDM.R
#' @include ensemble.R stacking.R
#' @include modelling.R ensemble_modelling.R stack_modelling.R
#' @include save.model.R load_model.R plot.model.R
#' @importFrom sp spplot SpatialMultiPoints
#' @importFrom raster crop extent aggregate
#' @import shiny shinydashboard raster
#' @importFrom gplots heatmap.2
#@importFrom raster stack reclassify extent
NULL

#' Number one package Global User Interface
#'
#' User interface to use the number one package
#'
#' @param browser logical. Option to plot or not the user global interface in
#'   your internet browser
#' @param maxmem numeric. Option to choose the maximum memory allocated to the
#'   user interface
#'
#' @return Open a window with a shiny app allowing to use the number one package
#'   with an easy user interface
#'
#' @details Due to the relatively important size of environmental data you
#'   should gave enough memory to the interface
#'
#' @examples
#' \dontrun{
#' gui()
#' }
#'
#' @export
gui = function (browser = F, maxmem = 10e9) {

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
                    box(title = 'Environment variable', height = 600,
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
                    box(title = 'Occurrences table',
                        fileInput('Occ', 'Occurrences', multiple = F, accept = '.csv'),
                        textInput('sep','Separator'),
                        textInput('dec','Decimal'),
                        uiOutput('Xcol'),
                        uiOutput('Ycol'),
                        uiOutput('Pcol'),
                        uiOutput('Spcol'),
                        checkboxInput('GeoRes', 'Geographic resampling', value = T),
                        uiOutput('reso'),
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
                fluidRow(
                  tabBox(
                    tabPanel('Easy'
                    ),
                    tabPanel('Intermediate'
                    ),
                    tabPanel('Professional'
                    )
                  )
                )
        ),
        tabItem('ensemblem',
                fluidRow(
                  tabBox(
                    tabPanel('Easy'
                    ),
                    tabPanel('Intermediate'
                    ),
                    tabPanel('Professional'
                    )
                  )
                )
                ),
        tabItem('stackm',
                fluidRow(
                  tabBox(
                    tabPanel('Easy',
                             selectInput('algo', 'Algorithms', c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM'), multiple = T, selectize = T),
                             sliderInput('rep','Repetitions',1,50,10, step = 1),
                             textInput('name','Name'),
                             checkboxInput('uncert', 'Uncertainty mapping', value = T),
                             selectInput('metric', 'Evaluation metric', c('Kappa','CCR','TSS','SES','LW','ROC'), selected = 'SES'),
                             selectInput('method', 'Diversity mapping method', c('Probability','Random Bernoulli','Threshold'), selected = 'Probability'),
                             uiOutput('repBslide')
                    ),
                    tabPanel('Intermediate',
                             checkboxInput('PA', 'Automatic Pseudo-Absences', value = T),
                             uiOutput('PAnbUI'),
                             uiOutput('PAstratUI'),
                             selectInput('cval', 'Cross validation', c('holdout','k-fold','LOO'), selected = 'holdout'),
                             uiOutput('cvalparam1UI'),
                             uiOutput('cvalparam2UI'),
                             selectInput('axesmetric', 'Axes evaluation metric', c('Pearson','AUC','Kappa','sensitivity','specificity','prop.correct'), selected = 'Pearson'),
                             selectInput('ensemblemetric', 'Ensemble selection metric', c('AUC','Kappa','sensitivity','specificity','prop.correct'), selected = 'AUC', multiple = T, selectize = T),
                             uiOutput('AUCUI'),
                             uiOutput('KappaUI'),
                             uiOutput('sensitivityUI'),
                             uiOutput('specificityUI'),
                             uiOutput('propcorrectUI'),
                             checkboxInput('weight', 'Ensemble weighted', value = T)
                    ),
                    tabPanel('Expert',
                             sliderInput('thresh','Threshold precision',11,10001,101, step = 10),
                             checkboxInput('tmp', 'Temporary files', value = F),
                             checkboxGroupInput('algoparam', 'Algortims parameters', c('GLM','GAM','MARS','GBM','CTA','RF','ANN','SVM'), inline = T),
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
                  box(title = 'Preview', plotOutput('modelprev'))
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
                         tabPanel(plotOutput('enm.uncertainty'), title = 'Uncertainty'),
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

    # Environmental variable loading
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
        format = paste0('.',strsplit(load.var$vars[[i]], '.', fixed = T)[[1]][2])
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
                              load_var(directory = {},
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
        if(!is.null(input$layer)){
          i = as.numeric(which(as.list(names(data$Env)) == input$layer))
          spplot(data$Env[[i]],
                        main = names(data$Env)[i],
                        xlab = 'Longitude (\u02DA)',
                        ylab = 'Latitude (\u02DA)')
        }
        })
      })

    # Occurrences loading
    load.occ <- reactiveValues(columns = c())
    observeEvent(input$Occ, {
      load.occ$columns = names(read.csv2(input$Occ$datapath))
    })
    observeEvent(input$sep, {
      if(is.null(input$sep)) {sep = ""} else {sep = input$sep}
      if(!is.null(input$Occ$datapath)) {load.occ$columns = names(read.csv2(input$Occ$datapath, sep = input$sep))}
    })
    output$Xcol <- renderUI({selectInput('Xcol', 'X column', load.occ$columns, multiple = F)})
    output$Ycol <- renderUI({selectInput('Ycol', 'Y column', load.occ$columns, multiple = F)})
    output$Pcol <- renderUI({selectInput('Pcol', 'Presence column', c('None', load.occ$columns), multiple = F)})
    output$Spcol <- renderUI({selectInput('Spcol', 'Specie column', c('None', load.occ$columns), multiple = F)})
    output$reso <- renderUI({if(input$GeoRes){sliderInput('reso', 'Resmpling gird coefficient', 1,10,1)}})
    observeEvent(input$load2, {
      validate(
        need(length(data$Env@layers) > 0, 'You need to load environmental variable before !'),
        need(length(input$Occ$datapath) > 0, 'Choose occurrences file first !')
      )
      if(is.null(input$dec)) {dec = ","} else {dec = input$dec}
      if (input$Spcol == 'None') {Spcol = NULL} else {Spcol = input$Spcol}
      if(is.null(input$sep)) {sep = ""} else {sep = input$sep}
      data$Occ = withProgress(message = 'Occurrences loading',
                              load_occ(directory = {},
                                       file = as.character(input$Occ$datapath),
                                       Xcol = input$Xcol,
                                       Ycol = input$Ycol,
                                       Spcol = Spcol,
                                       GeoRes = input$GeoRes,
                                       reso = max(res(data$Env@layers[[1]])) * as.numeric(input$reso),
                                       verbose = F,
                                       GUI = T,
                                       sep = sep,
                                       dec = dec))
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
      if (input$model.type == 'Ensemble model') {data$ENM = load_enm(input$dir, directory = data$dir)}
      if (input$model.type == 'Stack species model') {
        data$Stack = withProgress(message = 'Model loading',
                                  load_stack(name = input$dir, directory = data$dir, GUI = T))
        }
      output$model.preview <- renderPlot({
        if (input$model.type == 'Ensemble model') {
          spplot(data$ENM@projection,
               main = data$ENM@name,
               xlab = 'Longitude (\u02DA)',
               ylab = 'Latitude (\u02DA)',
               legend.args=list(text='Presence\nprobability', font = 3, line = 1))
        }
        if (input$model.type == 'Stack species model') {
          spplot(data$Stack@diversity.map,
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
    output$repBslide <- renderUI({if(input$method=='Random Bernoulli'){sliderInput('repB','Bernoulli repetitions',1,10000,1000, step = 1)}})
    output$PAnbUI <- renderUI(if(!input$PA){sliderInput('PAnb','PA number',1,10000,1000, step = 100)})
    output$PAstratUI <- renderUI(if(!input$PA){selectInput('PAstrat', 'PA strategy', c('random','disk'), selected = 'random')})
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
    output$cvalparam2UI <- renderUI({if(input$cval != 'LOO'){sliderInput('cvalrep', 'Cross validation repetitions', 1,20,1, step = 1)}})
    output$AUCUI <- renderUI({if('AUC' %in% input$ensemblemetric){sliderInput('AUC','AUC threshold',0.5,1,0.75, step = 0.05)}})
    output$KappaUI <- renderUI({if('Kappa' %in% input$ensemblemetric){sliderInput('Kappa','Kappa threshold',0,1,0.5, step = 0.05)}})
    output$sensitivityUI <- renderUI({if('sensitivity' %in% input$ensemblemetric){sliderInput('sensitivity','Sensitivity threshold',0.5,1,0.75, step = 0.05)}})
    output$specificityUI <- renderUI({if('specificity' %in% input$ensemblemetric){sliderInput('specificity','Specificity threshold',0.5,1,0.75, step = 0.05)}})
    output$propcorrectUI <- renderUI({if('prop.correct' %in% input$ensemblemetric){sliderInput('propcorrect','Correct proportion threshold',0.5,1,0.75, step = 0.05)}})
    output$testUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam){selectInput('test','Test (GLM/GAM)',c('AIC'))}})
    output$epsilonUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'SVM' %in% input$algoparam){sliderInput('epsilon','Epsilon (GLM/GAM/SVM) = 1e-X',1,20,8, step = 1)}})
    output$maxitUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'ANN' %in% input$algoparam){sliderInput('maxit','Maximum iteration (GLM/GAM/ANN)',100,1000,500, step = 100)}})
    output$degreeUI <- renderUI({if('MARS' %in% input$algoparam){sliderInput('degree','Degree of interaction (MARS)',1,10,2, step = 1)}})
    output$threshshrinkUI <- renderUI({if('GBM' %in% input$algoparam){sliderInput('threshshrink','Shrinkage threshold (GBM)',1e-05,1,1e-03, step = 1e-04)}})
    output$treesUI <- renderUI({if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam){sliderInput('trees','Maximum trees (GBM/RF)',500,10000,2500, step = 500)}})
    output$finalleaveUI <- renderUI({if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam || 'CTA' %in% input$algoparam){sliderInput('finalleave','Final leave size (GBM/CTA/RF)',1,20,1, step = 1)}})
    output$cvUI <- renderUI({if('GBM' %in% input$algoparam || 'CTA' %in% input$algoparam || 'SVM' %in% input$algoparam){sliderInput('cv','Cross validation number (GBM/CTA/SVM)',1,20,3, step = 1)}})
    observeEvent(input$model,{
      if (input$Pcol == 'None') {Pcol = NULL} else {Pcol = input$Pcol}
      if (input$Spcol == 'None') {Spcol = NULL} else {Spcol = input$Spcol}
      if (input$name == "") {name = NULL} else {name = input$name}
      algo = c()
      for (i in 1:length(input$algo)){algo = c(algo, input$algo[i])}
      if(!is.null(input$PA)) {PA = T}
      if(input$PA) {PA = NULL} else {PA = list('nb' = input$PAnb, 'strat' = input$PAstrat)}
      ensemble.metric = c()
      ensemblethresh = c(input$AUC, input$Kappa, input$sensitivity, input$specificity, input$propcorrect)
      ensemble.thresh = c()
      for (i in 1:length(input$ensemblemetric)){
        ensemble.metric = c(ensemble.metric, input$ensemblemetric[i])
        ensemble.thresh = c(ensemble.thresh, as.numeric(ensemblethresh[i]))
      }
      algoparam = list()
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam){algoparam$test = input$test}
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'SVM' %in% input$algoparam){algoparam$epsilon = as.numeric(input$epsilon)}
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'ANN' %in% input$algoparam){algoparam$maxit = as.numeric(input$maxit)}
      if('MARS' %in% input$algoparam){algoparam$degree = as.numeric(input$degree)}
      if('GBM' %in% input$algoparam){algoparam$thresh.shrink = as.numeric(input$threshshrink)}
      if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam){algoparam$trees = as.numeric(input$trees)}
      if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam || 'CTA' %in% input$algoparam){algoparam$final.leave = as.numeric(input$finalleave)}
      if('GBM' %in% input$algoparam || 'CTA' %in% input$algoparam || 'SVM' %in% input$algoparam){algoparam$cv = as.numeric(input$cv)}

      method = switch(input$method,
                      'Probability' = 'P',
                      'Random Bernoulli' = 'B',
                      'Threshold' = 'T')
      if(!is.null(as.numeric(input$repB))) {rep.B = 1000} else {rep.B = as.numeric(input$repB)}
      if(!is.null(c(as.numeric(input$cvalparam1), as.numeric(input$cvalrep)))) {cv.param = c(0.7, 1)} else {cv.param = c(as.numeric(input$cvalparam1), as.numeric(input$cvalrep))}
      data$Stack = stack_modelling(algo,
                                   data$Occ, data$Env,
                                   Xcol = input$Xcol,
                                   Ycol = input$Ycol,
                                   Pcol = Pcol,
                                   Spcol = Spcol,
                                   rep = as.numeric(input$rep),
                                   name = name,
                                   save = F,
                                   directory = 'nowhere',
                                   PA = PA,
                                   cv = input$cval,
                                   cv.param = cv.param,
                                   thresh = as.numeric(input$thresh),
                                   axes.metric = input$axesmetric,
                                   uncertainty =  input$uncert,
                                   tmp =  input$tmp,
                                   ensemble.metric = ensemble.metric,
                                   ensemble.thresh = ensemble.thresh,
                                   weight = input$weight,
                                   method = method,
                                   metric = input$metric,
                                   rep.B =  rep.B,
                                   verbose = T,
                                   GUI = F,
                                   algoparam)
      output$modelprev = renderPlot(image(data$Stack@diversity.map))
      })

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
      eval = character()
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
           legend.args=list(text='Local \nspecies \nrichness', font = 3, line = 1))
    })
    output$Uncertainity <- renderPlot({
      eval = character()
      ensemble.metric = strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = T)[[1]][-1]
      for (i in 1:length(ensemble.metric)) {
        eval = paste(eval, paste(ensemble.metric[i],':',round(data$Stack@evaluation[1,which(names(data$Stack@evaluation) == ensemble.metric[i])], digits = 3)))
        if (i < length(ensemble.metric)) {eval = paste(eval, ',')}
      }
      if (!is.null(ranges$x)) {uncert = crop(data$Stack@uncertainty, c(ranges$x, ranges$y))} else {uncert = data$Stack@uncertainty}
      spplot(uncert, main = eval, legend.args=list(text='Models \nvariance', font = 3, line = 1))})
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
      varimp = as.data.frame(t(data$Stack@variable.importance))
      names(varimp) = 'Axes.evaluation'
      bar = barplot(varimp$Axes.evaluation, names.arg = strwrap(row.names(varimp)),
                    ylim = c(0,(max(varimp$Axes.evaluation)+max((varimp[2])))),
                    las = 2, ylab = 'Variable relative contribution (%)')
      arrows(bar,varimp$Axes.evaluation+varimp[,2], bar, varimp$Axes.evaluation-varimp[,2], angle=90, code=3, length=0.1)
    })
    output$varimp.table <- renderTable({data$Stack@variable.importance})
    # Parameters
    output$summary <- renderTable({
      summary = data.frame(matrix(nrow = 5, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Final number of species', 'Original algorithms', 'Number of repetitions', 'Pseudo-absences selection')
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
        if(!is.null(ranges$x)){ref = crop(data$ENM@projection, c(ranges$x, ranges$y))} else {ref = data$ENM@projection}
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
      if (!is.null(ranges$x)) {proba = crop(data$ENM@projection, c(ranges$x, ranges$y))} else {proba = data$ENM@projection}
      spplot(proba,
           main = paste('AUC :',round(data$ENM@evaluation$AUC,3),'  Kappa',round(data$ENM@evaluation$Kappa,3)),
           xlab = 'Longitude (\u02DA)',
           ylab = 'Latitude (\u02DA)',
           legend.args=list(text='Presence\nprobability', font = 3, line = 1),
           col.regions = rev(terrain.colors(10000)))
      spplot(SpatialMultiPoints(c(data$ENM@data$X[which(data$ENM@data$Presence == 1)],
             data$ENM@data$Y[which(data$ENM@data$Presence == 1)])),
             pch = 16, cex = 0.7, add = T)
      # spplot(M@projection, sp.layout=list(SpatialPoints(points), pch = 16, cex = 0.7, col = 'black'), col.regions = rev(terrain.colors(10000)))
    })
    output$niche <- renderPlot({
      niche.map = reclassify(data$ENM@projection, c(-Inf,data$ENM@evaluation$threshold,0, data$ENM@evaluation$threshold,Inf,1))
      if (!is.null(ranges$x)) {niche.map = crop(niche.map, c(ranges$x, ranges$y))}
      spplot(niche.map, main = paste('AUC :',round(data$ENM@evaluation$AUC,3),'  Kappa',round(data$ENM@evaluation$Kappa,3)))})
    output$enm.uncertainty <- renderPlot({
      if (!is.null(ranges$x)) {uncert.map = crop(data$ENM@uncertainty, c(ranges$x, ranges$y))} else {uncert.map = data$ENM@uncertainty}
      spplot(uncert.map, main = paste('AUC :',round(data$ENM@evaluation$AUC,3),'  Kappa',round(data$ENM@evaluation$Kappa,3)), legend.args=list(text='Models \nvariance', font = 3, line = 1))})
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
      varimp = as.data.frame(t(data$ENM@variable.importance[-1]))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = row.names(varimp), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$enm.varimp.table <- renderTable({data$ENM@variable.importance[-1]})
    # Parameters
    output$enm.summary <- renderTable({
      summary = data.frame(matrix(nrow = 6, ncol = 1))
      names(summary) = 'Summary'
      row.names(summary) = c('Occurrences type', 'Occurrences number', 'Final number of species', 'Original algorithms', 'Number of repetitions', 'Pseudo-absences selection')
      algo.info = character()
      for (i in 1:length(strsplit(data$ENM@parameters$algorithms, '.', fixed = T)[[1]][-1])) {
        algo.info = paste(algo.info, strsplit(data$ENM@parameters$algorithms, '.', fixed = T)[[1]][-1][i])
      }
      if (data$ENM@parameters$PA) {PA = 'default'}
      if (data$ENM@parameters$data == "presence-only data set") {
        nb.occ =  length(as.factor(data$ENM@data$Presence[which(data$ENM@data$Presence==1)])) / sum(data$ENM@algorithm.evaluation$kept.model)
      } else {
        nb.occ =  length(as.factor(data$ENM@data$Presence)) / sum(data$ENM@algorithm.evaluation$kept.model)
      }
      summary$Summary = c(data$ENM@parameters$data, nb.occ, length(data$Stack@enms), algo.info, data$ENM@parameters$rep, PA)
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
