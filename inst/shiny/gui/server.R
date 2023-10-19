if (!requireNamespace("shinyFiles", quietly = TRUE)) {
  stop("shinyFiles package is needed for the Global User Interface. Please install it.")
}
library(shinydashboard)
library(shinyFiles)
library(raster)
library(ggplot2)


serverWD <- function(working.directory){
  function(input, output, session) {
    ### Server data ###
    data <- reactiveValues(Env = stack(), Occ = data.frame(), dir = getwd(), ESDM = NULL, esdms = list(), Stack = NULL)
    result <- reactiveValues(ESDM = NULL)

    ### Menu ###
    output$menu <- renderMenu({
      sidebarMenu(
        id = 'actions',
        menuItem('Welcome page', tabName = 'welcomepage',selected=TRUE),
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
        if(length(data$Env) > 0 & any(length(data$Stack) > 0, length(data$ESDM) > 0)){
          menuItem('Model forecasting',tabName='forecasting')
        },
        if(!is.null(data$Stack) || !is.null(data$ESDM)) {
          if(inherits(data$ESDM,'Algorithm.SDM')) {tabname = 'Algorithm SDM'} else {tabname = 'Ensemble SDM'}
          menuItem('Results',
                   if(!is.null(data$Stack)){menuItem("SSDM", tabName = "stack", icon = icon("dashboard"))},
                   menuItem(tabname, tabName = "stackesdm", icon = icon("pagelines")),
                   if(!is.null(data$Stack)){selectInput('esdmchoice', 'Species:', data$esdms, selectize = TRUE)},
                   if(!inherits(data$ESDM,'Algorithm.SDM')) {menuItem('Save model', tabName = "save", icon = icon("floppy-o"))},
                   if(!is.null(data$ESDM) | !is.null(data$Stack)) {
                     menuItem('Save maps', tabName = "savem", icon = icon("floppy-o"))
                     }
          )
        },
        menuItem('Quit', tabName = 'quitpage')
      )
    })

    ### Load Menu ###

    ## Load new data page ##

    # Environmental variable loading
    load.var <- reactiveValues(factors = c(), formats = c(), norm = TRUE,  vars = list())
    # working.directory <- get("working.directory", envir = .GlobalEnv)
    example = system.file("extdata", package = "SSDM")
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyFileChoose(input, 'envfiles', session=session,
                      roots=c(wd = working.directory,
                              example = example,
                              home = '/home',
                              root = '/'),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyFileChoose(input, 'envfiles', session=session,
                      roots=c(wd = working.directory,
                              example = example,
                              disks),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    } else {
      shinyFileChoose(input, 'envfiles', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                home = '/user',
                                root = '/'),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    }
    observeEvent(input$envfiles,{
      if(!is.integer(input$envfiles)){
          load.var$vars = lapply(input$envfiles$files, function(x) x[[length(x)]])
          names(load.var$vars) <- unlist(load.var$vars)
      }
    })
    output$envnames <- renderTable({
      matrix(names(data$Env),dimnames=list(c(1:length(names(data$Env))),c("Selected data")))
      })

    output$factors <- renderUI({
      selectInput('factors', 'Categorical', load.var$vars, multiple = TRUE, selectize = TRUE)
    })
    observeEvent(input$load, {
      validate(
        need(length(load.var$vars) > 0, 'Choose environment variable files first !')
      )
      if(Sys.info()[['sysname']] == 'Linux') {
        path = switch(input$envfiles$root,
                      'wd' = working.directory,
                      'example' = example,
                      'home' = '/home',
                      'root' = '/')
      } else if (Sys.info()[['sysname']] == 'Windows') {
        path = switch(input$envfiles$root,
                      'wd' = working.directory,
                      'example' = example,
                      input$envfiles$root)
      } else {
        path = switch(input$envfiles$root,
                      'wd' = working.directory,
                      'example' = example,
                      'home' = '/home',
                      'root' = '/')
      }
      for(i in 2:(length(input$envfiles$files[[1]]))-1){
        path = paste0(path, '/', input$envfiles$files[[1]][i])
      }
      load.var$formats = c()
      for (i in seq_len(length(load.var$vars))) {
        format = paste0('.',strsplit(load.var$vars[[i]], '.', fixed = TRUE)[[1]][2])
        if (!(format %in% load.var$formats)) {load.var$formats = c(load.var$formats, format)}
      }
      if('Normalization' %in% input$load.var.options) {
        load.var$norm = TRUE
      } else {
        load.var$norm = FALSE
      }
      a = try(withProgress(message = 'Variables loading',
                           load_var(path,
                                    files = unlist(load.var$vars),
                                    format = load.var$formats,
                                    Norm = load.var$norm,
                                    tmp = FALSE,
                                    categorical = load.var$factors,
                                    verbose = FALSE,
                                    GUI = TRUE)))
      if(inherits(a, 'try-error')){
        output$Envbug <- renderUI(p('Environmental variables loading failed, please check your inputs and try again'))
      } else {
        output$Envbug <- renderUI(p())
        data$Env = a
        for (i in seq_len(length(load.var$vars))) {
          names(data$Env)[i] = strsplit(load.var$vars[[i]], '.', fixed = TRUE)[[1]][1]
        }
        output$layerchoice <- renderUI({
          selectInput('layer', 'Variable', as.list(names(data$Env)), multiple = FALSE, selectize = TRUE)
        })
        output$env <- leaflet::renderLeaflet({
          if(!is.null(input$layer)){
            i = as.numeric(which(as.list(names(data$Env)) == input$layer))
            if(data$Env[[i]]@data@isfactor) {
              map = !as.factor(data$Env[[i]])
            } else {
              map = data$Env[[i]]
            }
            pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                         values(map), na.color = "transparent")
            leaflet::leaflet() %>%
            leaflet::addTiles() %>%
            leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
            leaflet::addLegend(pal = pal, values = values(map), title = "Diversity")
          }
        })
      }
      updateTabItems(session, "actions", selected = "newdata")
    })

    # Occurrences loading
    load.occ <- reactiveValues(columns = c())
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyFileChoose(input, 'Occ', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                home = '/home',
                                root = '/'),
                      filetypes=c('',"csv", "txt"))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyFileChoose(input, 'Occ', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                disks),
                      filetypes=c('',"csv", "txt"))
    } else {
      shinyFileChoose(input, 'Occ', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                home = '/user',
                                root = '/'),
                      filetypes=c('',"csv", "txt"))
    }
    observeEvent(input$Occ, {
      if(!is.integer(input$Occ)) {
        file = paste0(switch(input$Occ$root,
                             'wd' = working.directory,
                             'example' = example,
                             'home' = '/home',
                             'root' = '/',
                             input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
        load.occ$columns = names(read.csv2(file))
      }
    })
    observeEvent(input$sep, {
       if(!is.integer(input$Occ)) {
        file = paste0(switch(input$Occ$root,
                             'wd' = working.directory,
                             'example' = example,
                             'home' = '/home',
                             'root' = '/',
                             input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
        load.occ$columns = names(read.csv2(file, sep = input$sep, nrows = 0))
      }
    })
    observeEvent(input$Occ, {
      if(!is.integer(input$Occ)) {
        file = paste0(switch(input$Occ$root,
                             'wd' = working.directory,
                             'example' = example,
                             'home' = '/home',
                             'root' = '/',
                             input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
        load.occ$columns = names(read.csv2(file, sep = input$sep, nrows = 0))
      }
    })
    output$Xcol <- renderUI({selectInput('Xcol', 'X column', load.occ$columns, multiple = FALSE)})
    output$Ycol <- renderUI({selectInput('Ycol', 'Y column', load.occ$columns, multiple = FALSE)})
    output$Pcol <- renderUI({selectInput('Pcol', 'Presence (1) / absence (0) column', c('None', load.occ$columns), multiple = FALSE)})
    output$Spcol <- renderUI({selectInput('Spcol', 'Species column', c('None', load.occ$columns), multiple = FALSE)})
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
                           'wd' = working.directory,
                           'example' = example,
                           'home' = '/home',
                           'root' = '/',
                           input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
      a = try(withProgress(message = 'Occurrences loading',
                           load_occ(path = {},
                                    data$Env,
                                    file,
                                    Xcol = input$Xcol,
                                    Ycol = input$Ycol,
                                    Spcol = Spcol,
                                    GeoRes = input$GeoRes,
                                    reso = max(res(data$Env@layers[[1]])) * as.numeric(input$reso),
                                    verbose = FALSE,
                                    GUI = TRUE,
                                    sep = sep,
                                    dec = dec)))
      if(inherits(a, 'try-error')){
        output$Occbug <- renderUI(p('Occurrences loading failed, please check your inputs and try again'))
      } else {
        output$Occbug <- renderUI(p(' '))
        data$Occ = a
      }
      updateTabItems(session, "actions", selected = "newdata")
    })
    output$occ <- renderDataTable({if(length(data$Occ) > 0) {data$Occ}})

    ## Load previous model page ##
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyDirChoose(input, 'prevmodel', session=session,
                     roots = c(wd = working.directory,
                               home = '/home',
                               root = '/'),
                     filetypes=c(''))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyDirChoose(input, 'prevmodel', session=session,
                     roots = c(wd = working.directory,
                               disks),
                     filetypes=c(''))
    } else {
      shinyDirChoose(input, 'prevmodel', session=session,
                     roots = c( wd= working.directory,
                                home = '/user',
                                root = '/'),
                     filetypes=c(''))
    }
    observeEvent(input$load.model, {
      validate(
        need(input$model.type != ' ', 'You need to choose the model type first !')
      )
      if(Sys.info()[['sysname']] == 'Linux') {
        path = switch(input$prevmodel$root,
                      'wd' = working.directory,
                      'home' = '/home/',
                      'root' = '/')
      } else if (Sys.info()[['sysname']] == 'Windows') {
        path = switch(input$prevmodel$root, 'wd' = working.directory, input$prevmodel$root)
      } else {
        path = switch(input$prevmodel$root,
                      'wd' = working.directory,
                      'home' = '/home',
                      'root' = '/')
      }
      for(i in 2:(length(input$prevmodel$path)-1)){
        path = paste0(path, '/', input$prevmodel$path[[i]][1])
      }
      name = input$prevmodel$path[[length(input$prevmodel$path)]][1]
      if (input$model.type == 'Ensemble SDM') {
        a = try(load_esdm(name, path))
      }
      if (input$model.type == 'SSDM') {
        a = try(withProgress(message = 'Model loading', load_stack(name, path, GUI = TRUE)))
      }
      if(inherits(a, 'try-error')){
        output$prevmodelbug <- renderText('Previous model loading failed, please check your inputs and try again')
      } else {
        output$prevmodelbug <- renderText(' ')
        if (input$model.type == 'Ensemble SDM') {
          data$ESDM = a
          if(!is.null(data$ESDM)){result$ESDM = data$ESDM}
          output$model.preview <- leaflet::renderLeaflet({
            diversity <- data$ESDM@projection
            pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                         values(diversity), na.color = "transparent")
            leaflet::leaflet() %>%
              leaflet::addTiles() %>%
              leaflet::addRasterImage(diversity, colors = pal, opacity = 0.8) %>%
              leaflet::addLegend(pal = pal, values = values(diversity), title = "Diversity")
          })
        }
        if (input$model.type == 'SSDM') {
          data$Stack = a
          for (i in seq_len(length(names(data$Stack@esdms)))) {data$esdms[[i]] = strsplit(names(data$Stack@esdms), '.Ensemble.SDM', fixed = TRUE)[[i]][1]}
          output$model.preview <- leaflet::renderLeaflet({
            diversity <- data$Stack@diversity.map
            pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                         values(diversity), na.color = "transparent")
            leaflet::leaflet() %>%
              leaflet::addTiles() %>%
              leaflet::addRasterImage(diversity, colors = pal, opacity = 0.8) %>%
              leaflet::addLegend(pal = pal, values = values(diversity), title = "Diversity")
          })
        }
      }
    })

    ## Basic modelling parameters ##
    output$species <-renderUI({
      if(input$modellingchoice != 'Stack modelling'){
        if (input$Spcol != 'None') {
          selectInput('species', 'Species', as.list(levels(data$Occ[,which(names(data$Occ) == input$Spcol)])))
        }
      }
    })
    output$algoUI <- renderUI({
      if(input$modellingchoice == 'Algorithm modelling'){
        title = 'Algorithm'
        multiple = FALSE
      } else {
        title = 'Algorithms'
        multiple = TRUE
      }
      selectInput('algo', title, c('GLM','GAM','MARS','GBM','CTA','RF','MAXENT','ANN','SVM'), multiple = multiple, selectize = TRUE)})
    output$repUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){sliderInput('rep','Repetitions',1,50,10, step = 1)}})
    output$nameUI <- renderUI({textInput('name','Name')})
    output$uncertUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){checkboxInput('uncert', 'Uncertainty mapping', value = TRUE)}})
    output$uncertinfoUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){p('Between-algorithm variance map')}})
    output$endemismUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){selectInput('endemism', 'Endemism mapping', c('None', 'WEI', 'CWEI'), selected = 'WEI')}})
    output$endemisminfoUI <- renderUI({
      if(input$modellingchoice == 'Stack modelling'){
        if(length(input$endemism) == 1) {
          p(switch(input$endemism,
                   'None' = 'No endemism map will be built',
                   'WEI' = 'Endemism map will be built by counting all species in each cell and weighting each by the inverse of its number of occurrences',
                   'CWEI' = 'Endemism map will be built by dividing the weighted endemism index by the total count of species in the cell')
          )
        }
      }
    })
    output$endemismrangeUI <- renderUI({
      if(input$modellingchoice == 'Stack modelling'){
        if(length(input$endemism) == 1) {
          if(input$endemism != 'None'){
            selectInput('endemismrange', 'Range in endemism mapping', c('NbOcc', 'Binary'), selected = 'Binary')
          }
        }
      }
    })
    output$endemismrangeinfoUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){
      if(length(input$endemismrange) == 1) {
        p(switch(input$endemismrange,
                 'NbOcc' = 'Range in endemism index computing is the total number of occurrences.',
                 'Binary' = 'Range in endemism index comuting is the surface of the binary map species distribution.'))}}})
    output$metricUI <- renderUI({selectInput('metric', 'Evaluation metric', c('Kappa','CCR','TSS','SES','LW','ROC'), selected = 'SES')})
    output$metricinfoUI <- renderUI({
      if(length(input$metric) == 1) {
        p(switch(input$metric,
                 'Kappa' = 'Maximizes the Kappa',
                 'CCR' = 'Maximizes the sum of sensitivity and specificity',
                 'TSS' = '(True Skill Statistic) maximizes the sum of sensitivity and specificity',
                 'SES' = 'Uses the sensitivity-specificity equality',
                 'LW' = 'Uses the lowest occurrence prediction probability',
                 'ROC' = 'Minimizes the distance between the ROC plot (receiving operating curve) and the upper left corner (1,1).'))
      }})
    output$methodUI <- renderUI({
      if(input$modellingchoice == 'Stack modelling'){
        selectInput('method', 'Diversity mapping method',
                    c('pSSDM','Bernoulli','bSSDM','MaximumLikelihood','PRR.MEM','PRR.pSSDM'),
                    selected = 'Probability')
      }
    })
    output$methodinfoUI <- renderUI({
      if(input$modellingchoice == 'Stack modelling'){
        if(length(input$method) == 1) {
          p(switch(input$method,
                   'pSSDM' = 'Sum probabilities of habitat suitability maps',
                   'Bernoulli' = 'Drawing repeatedly from a Bernoulli distribution',
                   'bSSDM' = 'Sum the binary map obtained with the thresholding (depending on the metric, see metric parameter)',
                   'MaximumLikelihood' = 'Adjust species richness of the model by linear regression',
                   'PRR.MEM' = 'Model richness with a macroecological model (MEM) and adjust each ESDM binary map by ranking habitat suitability and keeping as much as predicted richness of the MEM',
                   'PRR.pSSDM' = 'Model richness with a pSSDM and adjust each ESDM binary map by ranking habitat suitability and keeping as much as predicted richness of the pSSDM'))
        }
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
    output$ensemblemetricUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){selectInput('ensemblemetric', 'Ensemble selection metric', c('AUC','Kappa','sensitivity','specificity','prop.correct'), selected = 'AUC', multiple = TRUE, selectize = TRUE)}})
    output$ensembleinfoUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){p('Metric(s) and threshold(s) used to select the best SDMs that will be kept in the ensemble SDM')}})
    output$AUCUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('AUC' %in% input$ensemblemetric){sliderInput('AUC','AUC threshold',0.5,1,0.75, step = 0.05)}}})
    output$KappaUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('Kappa' %in% input$ensemblemetric){sliderInput('Kappa','Kappa threshold',0,1,0.5, step = 0.05)}}})
    output$sensitivityUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('sensitivity' %in% input$ensemblemetric){sliderInput('sensitivity','Sensitivity threshold',0.5,1,0.75, step = 0.05)}}})
    output$specificityUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('specificity' %in% input$ensemblemetric){sliderInput('specificity','Specificity threshold',0.5,1,0.75, step = 0.05)}}})
    output$propcorrectUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){if('prop.correct' %in% input$ensemblemetric){sliderInput('propcorrect','Correct proportion threshold',0.5,1,0.75, step = 0.05)}}})
    output$weightUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){checkboxInput('weight', 'Ensemble weighted', value = TRUE)}})
    output$rangeUI <- renderUI({if(input$modellingchoice == 'Stack modelling'){checkboxInput('range', 'Range restriction', value = FALSE)}})
    output$rangevalUI <- renderUI({
      if(is.logical(input$range)){
        if(input$modellingchoice == 'Stack modelling' && input$range){
          sliderInput('rangeval', 'Range restriction value', 1,100,5, step = 1)
        }
      }
    })
    output$rangeinfoUI <- renderUI({
      if(is.logical(input$range)){
        if(input$modellingchoice == 'Stack modelling' && input$range){
          p('Set a value of range restriction (in pixels) around presence occurrences on habitat suitability maps (all further points will have a null probability)')
        }
      }
    })

    ## Advanced modelling parameters ##
    output$tmpUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){checkboxInput('tmp', 'Temporary files', value = FALSE)}})
    output$tmpinfoUI <- renderUI({if(input$modellingchoice != 'Algorithm modelling'){p('Rasters are saved in temporary files to release memory')}})
    output$testUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam){selectInput('test','Test (GLM/GAM)',c('AIC', 'BIC'))}})
    output$epsilonUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'SVM' %in% input$algoparam){sliderInput('epsilon','Epsilon (GLM/GAM/SVM) = 1e-X',1,20,8, step = 1)}})
    output$maxitUI <- renderUI({if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'ANN' %in% input$algoparam){sliderInput('maxit','Maximum number of iterations (GLM/GAM/ANN)',100,1000,500, step = 100)}})
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
      data$ESDM = NULL
      data$Stack = NULL
      data$esdms = list()
      result$ESDM = NULL
      output$modelprev = renderPlot({plot.new()})
      output$modelfailed = renderText({''})
      if (input$Pcol == 'None') {Pcol = NULL} else {Pcol = input$Pcol}
      if (input$Spcol == 'None') {Spcol = NULL} else {Spcol = input$Spcol}
      if (input$name == "") {name = NULL} else {name = input$name}
      algo = c()
      for (i in seq_len(length(input$algo))) {algo = c(algo, input$algo[i])}
      if(!is.null(input$PA)) {PA = TRUE}
      if(input$PA) {PA = NULL} else {PA = list('nb' = input$PAnb, 'strat' = input$PAstrat)}
      if(is.null(input$ensemblemetric)){
        ensemble.metric = c('AUC')
        ensemble.thresh = c(0.75)
      } else {
        ensemble.metric = c()
        ensemblethresh = c(input$AUC, input$Kappa, input$sensitivity, input$specificity, input$propcorrect)
        ensemble.thresh = c()
        for (i in seq_len(length(input$ensemblemetric))) {
          ensemble.metric = c(ensemble.metric, input$ensemblemetric[i])
          ensemble.thresh = c(ensemble.thresh, as.numeric(ensemblethresh[i]))
        }
      }
      if(is.null(input$range) || !input$range) {range = NULL} else {range = input$rangeval}
      if(is.null(input$endemism) || input$endemism == 'None') {endemism = NULL} else {endemism = c(input$endemism, input$endemismrange)}
      algoparam = list()
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam){algoparam$test = input$test} else {algoparam$test = 'AIC'}
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'SVM' %in% input$algoparam){algoparam$epsilon = as.numeric(paste0('1e-', as.character(input$epsilon)))} else {algoparam$epsilon = 1e-08}
      if('GLM' %in% input$algoparam || 'GAM' %in% input$algoparam || 'ANN' %in% input$algoparam){algoparam$maxit = as.numeric(input$maxit)} else {algoparam$maxit = 500}
      if('MARS' %in% input$algoparam){algoparam$degree = as.numeric(input$degree)} else {algoparam$degree = 3}
      if('GBM' %in% input$algoparam){algoparam$thresh.shrink = as.numeric(input$threshshrink)} else {algoparam$threshshrink = 1e-03}
      if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam){algoparam$trees = as.numeric(input$trees)} else {algoparam$trees = 2500}
      if('GBM' %in% input$algoparam || 'RF' %in% input$algoparam || 'CTA' %in% input$algoparam){algoparam$final.leave = as.numeric(input$finalleave)}  else {algoparam$finalleave = 1}
      if('GBM' %in% input$algoparam || 'CTA' %in% input$algoparam || 'SVM' %in% input$algoparam){algoparam$cv = as.numeric(input$cv)}  else {algoparam$cv = 3}
      if(is.null(input$repB)) {rep.B = 1000} else {rep.B = as.numeric(input$repB)}
      if(is.null(input$cval)) {cval = 'holdout'} else {cval = input$cval}
      if(length(c(as.numeric(input$cvalparam1), as.numeric(input$cvalrep))) < 2) {cv.param = c(0.7, 1)} else {cv.param = c(as.numeric(input$cvalparam1), as.numeric(input$cvalrep))}
      if(!inherits(input$tmp,'logical')) {tmp = FALSE} else {tmp = input$tmp}
      if(!inherits(input$weight,'logical')) {weight = TRUE} else {weight = input$weight}
      if(input$modellingchoice == 'Algorithm modelling'){
        if (input$Spcol != 'None') {
          Occ = data$Occ[which(data$Occ[,which(names(data$Occ) == input$Spcol)] == input$species),]
        } else {
          Occ = data$Occ
        }
        data$ESDM = withProgress(message = 'SDM',
                                modelling(algo,
                                          Occ, data$Env,
                                          Xcol = input$Xcol,
                                          Ycol = input$Ycol,
                                          Pcol = Pcol,
                                          name = name,
                                          save = FALSE,
                                          path = 'nowhere',
                                          PA = PA,
                                          cv = cval,
                                          cv.param = cv.param,
                                          thresh = as.numeric(input$thresh),
                                          axes.metric = input$axesmetric,
                                          select.metric = c('AUC'),
                                          select.thresh = c(0),
                                          select = FALSE,
                                          metric = input$metric,
                                          verbose = FALSE,
                                          GUI = TRUE,
                                          test = algoparam$test,
                                          maxit = algoparam$maxit,
                                          epsilon = algoparam$epsilon,
                                          degree = algoparam$degree,
                                          threshshrink = algoparam$threshshrink,
                                          trees = algoparam$trees,
                                          finalleave = algoparam$finalleave,
                                          algocv = algoparam$cv))
        output$modelprev = renderPlot(spplot(data$ESDM@projection,
                                             main = 'Habitat suitability map',
                                             xlab = 'Longitude (\u02DA)',
                                             ylab = 'Latitude (\u02DA)',
                                             col.regions = rev(terrain.colors(10000))))
        result$ESDM = data$ESDM
      }
      if(!inherits(input$uncert,'logical')) {uncert = TRUE} else {uncert = input$uncert}
      if(input$modellingchoice == 'Ensemble modelling'){
        if (input$Spcol != 'None') {
          Occ = data$Occ[which(data$Occ[,which(names(data$Occ) == input$Spcol)] == input$species),]
        } else {
          Occ = data$Occ
        }
        data$ESDM = withProgress(message = 'Ensemble SDM',
                                ensemble_modelling(algo,
                                                   Occ, data$Env,
                                                   Xcol = input$Xcol,
                                                   Ycol = input$Ycol,
                                                   Pcol = Pcol,
                                                   rep = as.numeric(input$rep),
                                                   name = name,
                                                   save = FALSE,
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
                                                   verbose = FALSE,
                                                   GUI = TRUE,
                                                   test = algoparam$test,
                                                   maxit = algoparam$maxit,
                                                   epsilon = algoparam$epsilon,
                                                   degree = algoparam$degree,
                                                   threshshrink = algoparam$threshshrink,
                                                   trees = algoparam$trees,
                                                   finalleave = algoparam$finalleave,
                                                   algocv = algoparam$cv))
        if(!is.null(data$ESDM)){
          output$modelprev = leaflet::renderLeaflet({
            diversity <- data$ESDM@projection
            pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                         values(diversity), na.color = "transparent")
            leaflet::leaflet() %>%
              leaflet::addTiles() %>%
              leaflet::addRasterImage(diversity, colors = pal, opacity = 0.8) %>%
              leaflet::addLegend(pal = pal, values = values(diversity), title = "Habitat\nsuitability")
          })
          result$ESDM = data$ESDM
        } else {
          output$modelfailed = renderText('No ensemble SDM were kept, maybe you should try lower ensemble threshold(s) ?')
        }
      }
      if(requireNamespace("parallel", quietly = TRUE)){
        cores = parallel::detectCores() - 1
      } else {
        cores = 1
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
                                                  save = FALSE,
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
                                                  method = input$method,
                                                  metric = input$metric,
                                                  rep.B =  rep.B,
                                                  range = range,
                                                  endemism = endemism,
                                                  verbose = FALSE,
                                                  GUI = TRUE,
                                                  cores = 0,
                                                  test = algoparam$test,
                                                  maxit = algoparam$maxit,
                                                  epsilon = algoparam$epsilon,
                                                  degree = algoparam$degree,
                                                  threshshrink = algoparam$threshshrink,
                                                  trees = algoparam$trees,
                                                  finalleave = algoparam$finalleave,
                                                  algocv = algoparam$cv))
        if(!is.null(data$Stack)){
          output$modelprev = leaflet::renderLeaflet({
            diversity <- data$Stack@diversity.map
            pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                         values(diversity), na.color = "transparent")
            leaflet::leaflet() %>%
              leaflet::addTiles() %>%
              leaflet::addRasterImage(diversity, colors = pal, opacity = 0.8) %>%
              leaflet::addLegend(pal = pal, values = values(diversity), title = "Species\nrichness")
          })
          for (i in seq_len(length(names(data$Stack@esdms)))) {data$esdms[[i]] = strsplit(names(data$Stack@esdms), '.Ensemble.SDM', fixed = TRUE)[[i]][1]}
        } else {
          output$modelfailed = renderText('You have less than two remaining ensemble SDMs, maybe you should try lower ensemble threshold(s) ?')
        }
      }
    })

    ## Stack results ##
    output$Stackname <- renderUI(h2(data$Stack@name, align = 'center'))
    output$diversity.title <- renderText({paste0('<b>Mean species richness error: ', round(data$Stack@evaluation['mean','species.richness.error'], 3), "<br>")})
    output$Diversity <- leaflet::renderLeaflet({
      map <- data$Stack@diversity.map
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(map), title = "Diversity")
    })
    output$endemism <- leaflet::renderLeaflet({
      map<- data$Stack@endemism.map
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(map), title = "Endemism")
    })
    output$Uncertainity <- leaflet::renderLeaflet({
      map <- data$Stack@uncertainty
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(map), title = "Uncertainty")
    })
    # Evaluation
    output$evaluation.barplot <- renderPlot({
      evaluation <- as.matrix(data$Stack@evaluation['mean',])
      colnames(evaluation)[1:2] <- c('richness', 'success')
      barplot(evaluation, col = rainbow(length(evaluation)), beside=TRUE)
    })
    output$evaluation.table <- renderTable({
      evaluation <- data$Stack@evaluation
      names(evaluation)[1:2] <- c('richness', 'success')
      evaluation
    })
    # Algorithms correlation
    output$algo.corr.table <- renderTable({data$Stack@algorithm.correlation})
    output$algo.corr.heatmap <- renderPlot({
      m <- as.matrix(data$Stack@algorithm.correlation)
      ggplot(reshape2::melt(m, id.vars = NULL),
             aes(Var1, Var2, fill = value, label = round(value, 2))) +
        geom_tile() +
        geom_text(col = "white") +
        scale_fill_gradient(guide = "none", na.value = "white",
                            low = scales::muted("blue"), high = scales::muted("red")) +
        theme_minimal() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(angle = 90))
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
      for (i in seq_len(length(strsplit(data$Stack@parameters$algorithms, '.', fixed = TRUE)[[1]][-1]))) {
        algo.info = paste(algo.info, strsplit(data$Stack@parameters$algorithms, '.', fixed = TRUE)[[1]][-1][i])
      }
      if (data$Stack@parameters$PA) {PA = 'default'}
      if(data$Stack@parameters$cv == 'LOO') {cv.param = 'None'}
      if(data$Stack@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                                  strsplit(data$Stack@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                                  'rep =',
                                                                  strsplit(data$Stack@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
      if(data$Stack@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                                 strsplit(data$Stack@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                                 'rep =',
                                                                 strsplit(data$Stack@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
      summary$Summary = c(data$Stack@parameters$data, length(data$Stack@esdms), algo.info, data$Stack@parameters$rep, PA, data$Stack@parameters$cv, cv.param)
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
      for (i in seq_len(length(strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1]))) {
        if (i == 1) {
          evaluation.info = paste0(evaluation.info, ' ', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],' (>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        } else if (i == length(data$Stack@parameters$axes.metric) && i != 1) {
          evaluation.info = paste0(evaluation.info, ' and ', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],' (>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste0(evaluation.info, ', ', strsplit(data$Stack@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],' (>',strsplit(data$Stack@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        }
      }
      evaluation.info
    })

    ## ESDM Stack Result ##
    observeEvent(input$esdmchoice, {
      if(!is.null(data$Stack)){
        if(length(data$esdms) > 0){result$ESDM = data$Stack@esdms[[which(data$esdms == input$esdmchoice)]]}
      }
    })
    output$probability.title <- renderText({
      paste0('<b>AUC :', round(result$ESDM@evaluation$AUC,3),
             '  Kappa', round(result$ESDM@evaluation$Kappa,3), "<br>")
    })
    output$probability <- leaflet::renderLeaflet({
      map <- result$ESDM@projection
      obs <- data.frame(X = result$ESDM@data$X[which(result$ESDM@data$Presence == 1)],
                        Y = result$ESDM@data$Y[which(result$ESDM@data$Presence == 1)])
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                   values(map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
        leaflet::addCircles(data = obs, lng = ~ X, lat = ~ Y, radius = 1) %>%
        leaflet::addLegend(pal = pal, values = values(map), title = "Habitat\nsuitability")
    })
    output$niche <- leaflet::renderLeaflet({
      map<- reclassify(result$ESDM@projection, c(-Inf,result$ESDM@evaluation$threshold,0, result$ESDM@evaluation$threshold,Inf,1))
      pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                  values(map), na.color = "transparent")
      leaflet::leaflet() %>%
        leaflet::addTiles() %>%
        leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
        leaflet::addLegend(pal = pal, values = values(map), title = "Niche")
    })
    output$esdm.uncertainty <- leaflet::renderLeaflet({
      if(!inherits(result$ESDM,'Algorithm.SDM')) {
        map <- result$ESDM@uncertainty
        pal <- leaflet::colorNumeric(rev(terrain.colors(1000)),
                                     values(map), na.color = "transparent")
        leaflet::leaflet() %>%
          leaflet::addTiles() %>%
          leaflet::addRasterImage(map, colors = pal, opacity = 0.8) %>%
          leaflet::addLegend(pal = pal, values = values(map), title = "Uncertainty")
      }
    })
    # Evaluation
    output$esdm.evaluation.barplot <- renderPlot({
      evaluation = result$ESDM@algorithm.evaluation
      if(!(0 %in% dim(evaluation))){
        for (i in seq_len(length(row.names(result$ESDM@algorithm.evaluation)))) {row.names(evaluation)[i] = strsplit(as.character(row.names(result$ESDM@algorithm.evaluation)[i]), '.SDM')[[1]][1]}
        for (i in seq_len(length(row.names(evaluation)))) {row.names(evaluation)[i] = tail(strsplit(as.character(row.names(evaluation)[i]), '.', fixed = TRUE)[[1]], n = 1)}
        evaluation$kept.model = evaluation$kept.model / max(evaluation$kept.model)
        table <- t(cbind(evaluation$AUC, evaluation$Kappa, evaluation$kept.model))
        barplot(table, col=c("darkblue","red","green"), names.arg = row.names(evaluation), beside=TRUE)
        legend('bottomright', c('AUC', 'Kappa','Kept model'), fill = c("darkblue","red","green"))
      }
    })
    output$esdm.evaluation.table <- renderTable({
      algo.eval = result$ESDM@algorithm.evaluation
      if(!(0 %in% dim(algo.eval))){
        for (i in seq_len(length(row.names(result$ESDM@algorithm.evaluation)))) {row.names(algo.eval)[i] = strsplit(as.character(row.names(result$ESDM@algorithm.evaluation)[i]), '.SDM')[[1]][1]}
        for (i in seq_len(length(row.names(result$ESDM@algorithm.evaluation)))) {row.names(algo.eval)[i] = tail(strsplit(as.character(row.names(algo.eval)[i]), '.', fixed = TRUE)[[1]], n = 1)}
        algo.eval[c(2,4:8)]
      }
    })

    # Algorithms correlation
    output$esdm.algo.corr.table <- renderTable({
      correlation = result$ESDM@algorithm.correlation
      for (i in seq_len(length(row.names(result$ESDM@algorithm.correlation)))) {row.names(correlation)[i] = strsplit(as.character(row.names(result$ESDM@algorithm.correlation)[i]), '.SDM')[[1]][1]}
      names(correlation) = row.names(correlation)
      if (length(correlation) > 0) {
        correlation[upper.tri(correlation, diag = TRUE)] = NA
        correlation
      }
    })
    output$esdm.algo.corr.heatmap <- renderPlot({
      correlation = result$ESDM@algorithm.correlation
      for (i in seq_len(length(row.names(result$ESDM@algorithm.correlation)))) {row.names(correlation)[i] = strsplit(as.character(row.names(result$ESDM@algorithm.correlation)[i]), '.SDM')[[1]][1]}
      names(correlation) = row.names(correlation)
      if (length(correlation) > 0) {
        correlation[upper.tri(correlation, diag = TRUE)] = NA
        m <- as.matrix(correlation)
        ggplot(reshape2::melt(m, id.vars = NULL),
               aes(Var1, Var2, fill = value, label = round(value, 2))) +
          geom_tile() +
          geom_text(col = "white") +
          scale_fill_gradient(guide = "none", na.value = "white",
                              low = scales::muted("blue"), high = scales::muted("red")) +
          theme_minimal() +
          theme(axis.title = element_blank(),
                axis.text.x = element_text(angle = 90))
      }
    })
    # Algo Eval Corr UI
    output$algoevalcorr <- renderUI({
      if(!inherits(result$ESDM,'Algorithm.SDM')) {
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
      }
    })
    # Variable importance
    output$esdm.varimp.barplot <- renderPlot({
      varimp = as.data.frame(t(result$ESDM@variable.importance))
      names(varimp) = 'Axes.evaluation'
      barplot(varimp$Axes.evaluation, names.arg = abbreviate(row.names(varimp)), las = 2, ylab = 'Variable relative contribution (%)')
    })
    output$esdm.varimp.table <- renderTable({result$ESDM@variable.importance})
    output$esdmvarimplegend <- renderTable({data.frame('Abbreviation' = abbreviate(names(result$ESDM@variable.importance)), 'Variable' = names(result$ESDM@variable.importance))})
    # Parameters
    output$esdm.summary <- renderTable({
      if(!inherits(result$ESDM,'Algorithm.SDM')) {
        summary = data.frame(matrix(nrow = 7, ncol = 1))
        names(summary) = 'Summary'
        row.names(summary) = c('Type of occurrences', 'Number of occurrences', 'Originally selected algorithms', 'Number of repetitions',
                               'Pseudo-absence selection method', 'Cross-validation method', 'Cross-validation parameters')
        algo.info = character()
        for (i in seq_len(length(strsplit(result$ESDM@parameters$algorithms, '.', fixed = TRUE)[[1]][-1]))) {
          algo.info = paste(algo.info, strsplit(result$ESDM@parameters$algorithms, '.', fixed = TRUE)[[1]][-1][i])
        }
        if (result$ESDM@parameters$PA) {PA = 'default'}
        if (result$ESDM@parameters$data == "presence-only data set") {
          nb.occ =  length(as.factor(result$ESDM@data$Presence[which(result$ESDM@data$Presence==1)])) / sum(result$ESDM@algorithm.evaluation$kept.model)
        } else {
          nb.occ =  length(as.factor(result$ESDM@data$Presence)) / sum(result$ESDM@algorithm.evaluation$kept.model)
        }
        if(result$ESDM@parameters$cv == 'LOO') {cv.param = 'None'}
        if(result$ESDM@parameters$cv == 'holdout') {cv.param = paste('fraction =',
                                                                    strsplit(result$ESDM@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                                    'rep =',
                                                                    strsplit(result$ESDM@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
        if(result$ESDM@parameters$cv == 'k-fold') {cv.param = paste('k =',
                                                                   strsplit(result$ESDM@parameters$cv.param, '|', fixed = TRUE)[[1]][2],
                                                                   'rep =',
                                                                   strsplit(result$ESDM@parameters$cv.param, '|', fixed = TRUE)[[1]][3])}
        summary$Summary = c(result$ESDM@parameters$data, nb.occ, algo.info, result$ESDM@parameters$rep, PA, result$ESDM@parameters$cv, cv.param)
        if(!is.null(result$ESDM@parameters$sp.nb.origin)) {
          summary = rbind(summary,
                          data.frame(Summary = result$ESDM@parameters$sp.nb.origin, row.names = 'Original number of species'))
        }
        summary
      }
      else {
        data.frame('Not computed')
      }
    })
    output$esdm.binary.info <- renderText({
      metric = switch(result$ESDM@parameters$metric,
                      'Kappa' = 'maximizing the Kappa',
                      'CCR' = 'maximizing the proportion of correctly predicted occurrences (CCR)',
                      'TSS' = 'maximizing the sensitivity and specificity sum (TSS)',
                      'SES' = 'using the sensitivity and specificity equality',
                      'LW' = 'using the lowest occurrence prediction probability',
                      'ROC' = 'minimizing the distance between the ROC plot and the upper left corner')
      text = paste('Binary map realized by', metric,
                   'with a final threshold of',  round(result$ESDM@evaluation$threshold, digits = 3))
      text
    })
    output$esdm.varimp.info <- renderText({
      varimp.info = 'Variable relative contribution evaluated with '
      for (i in seq_len(length(result$ESDM@parameters$axes.metric))) {
        if (i == 1) {
          varimp.info = paste(varimp.info, result$ESDM@parameters$axes.metric[i])
        } else if (i == length(result$ESDM@parameters$axes.metric) && i != 1) {
          varimp.info = paste(varimp.info, 'and', result$ESDM@parameters$axes.metric[i], '.')
        } else {
          varimp.info = paste(varimp.info, ',', result$ESDM@parameters$axes.metric[i])
        }
      }
      varimp.info
    })
    output$esdm.evaluation.info <- renderText({
      evaluation.info = 'Models evaluated with'
      for (i in seq_len(length(strsplit(result$ESDM@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1]))) {
        if (i == 1) {
          evaluation.info = paste0(evaluation.info, ' ', strsplit(result$ESDM@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],' (>',strsplit(result$ESDM@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        } else if (i == length(result$ESDM@parameters$axes.metric) && i != 1) {
          evaluation.info = paste0(evaluation.info, ' and ', strsplit(result$ESDM@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],' (>',strsplit(result$ESDM@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')','.')
        } else {
          evaluation.info = paste0(evaluation.info, ' , ', strsplit(result$ESDM@parameters$ensemble.metric, '.', fixed = TRUE)[[1]][-1][i],' (>',strsplit(result$ESDM@parameters$ensemble.thresh, '|', fixed = TRUE)[[1]][-1][i],')')
        }
      }
      if (result$ESDM@parameters$weight) {evaluation.info = paste(evaluation.info, ', and then weighted with the previous metrics means')}
      evaluation.info
    })

    ### Save Menu ###
    ## Save maps
    # select map to save
    output$speciesSave.sel <- renderUI({
      if(!is.null(data$Stack)){
      specieschoices <- reactive({c("all",names(data$Stack@esdms))})
      }
      if(is.null(data$Stack) & !is.null(data$ESDM)){
        specieschoices <-reactive({data$ESDM@name})
      }
      selectInput('speciesSave', 'Species:', specieschoices())
      })

    observeEvent(input$speciesSave,{
      if(input$speciesSave == "all"){
        mapchoices <- list("Diversity map" = "diversity.map", "Endemism map"="endemism.map")}
      if(input$speciesSave != "all" & !is.null(input$speciesSave)){
        mapchoices <- list("Ensemble map" = "projection")}
      if(input$uncert == TRUE){
        mapchoices <- c(mapchoices,"Uncertainty map" = "uncertainty")
      }
      output$mapSave.sel <- renderUI({
        selectInput('mapSave','Map type:', mapchoices)
      })
    })

    if(Sys.info()[['sysname']] == 'Linux') {
      shinyDirChoose(input, 'savem', session=session, roots=c( wd='.', home = '/home', root = '/'), filetypes=c(''))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyDirChoose(input, 'savem', session=session, roots=c( wd='.', disks), filetypes=c(''))
    } else {
      shinyDirChoose(input, 'savem', session=session, roots=c( wd='.', home = '/user', root = '/'), filetypes=c(''))
    }
    observeEvent(input$savemaps, {
      path = switch(input$savem$root,
                    'wd' = working.directory,
                    'home' = '/home',
                    'root' = '/',
                    input$savem$root)
      for(i in 2:length(input$savem$path)){
        path = paste0(path, '/', input$savem$path[[i]][1])
      }
      if(grepl("Ensemble",input$speciesSave) & !is.null(data$ESDM)) {
        # path <- "/home/lukas/Bilder/output/"
        writeRaster(eval(parse(text=paste0("data$ESDM@",input$mapSave))), filename = paste0(path,"/",data$ESDM@name,"_",input$mapSave), format="GTiff",overwrite=TRUE)}
      if(grepl("Ensemble",input$speciesSave) & !is.null(data$Stack)) {
        #path <- "/home/lukas/Bilder/output"
        writeRaster(eval(parse(text=paste0("data$Stack@esdms$",input$speciesSave,"@",input$mapSave))), filename = paste0(path,"/",input$speciesSave,"_",input$mapSave), format="GTiff",overwrite=TRUE)
      }
      if(input$speciesSave == "all") {
        # path <- "/home/lukas/Bilder/output/"
        writeRaster(eval(parse(text=paste0("data$Stack@",input$mapSave))), filename = paste0(path,"/",data$Stack@name,"_",input$mapSave), format="GTiff")}
        }
    )

    ## Save model page ##
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyDirChoose(input, 'save', session=session, roots=c( wd='.', home = '/home', root = '/'), filetypes=c(''))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
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
      path = switch(input$save$root,
                    'wd' = working.directory,
                    'home' = '/home',
                    'root' = '/',
                    input$save$root)
      for(i in 2:length(input$save$path)){
        path = paste0(path, '/', input$save$path[[i]][1])
      }
      if(!is.null(data$ESDM) && is.null(data$Stack)) {save.esdm(data$ESDM, name = data$ESDM@name, path, verbose = FALSE, GUI = TRUE)}
      if(!is.null(data$Stack)) {save.stack(data$Stack, name = data$Stack@name, path, verbose = FALSE, GUI = TRUE)}
    })

    ### Model forecasting

    output$projcheck <- renderText({
      if(is.null(data$Env) |  (is.null(data$ESDM) & is.null(data$Stack))){"Environmental data and model must both be provided"}
      else if(!is.null(data$ESDM)){if(!all(colnames(data$ESDM@data) %in% names(data$Env))){"Provided environmental variables do not match variables used model training. Please check your raster selection."}}
      else if(!is.null(data$Stack)){if(!all(colnames(data$Stack@variable.importance) %in% names(data$Env))){"Provided environmental variables do not match variables used for model training. Please check your raster selection."}}
      else if(!is.null(data$ESDM)){if(!is.null(data$Stack)){"Several models loaded. Please restart with only one model."}}
     })
    observeEvent(input$project,{
      if(!is.null(data$ESDM) & is.null(data$Stack)){
        data$ESDM <- project(data$ESDM,data$Env)
        updateTabItems(session, "actions", selected="stackesdm")
      }
      if(is.null(data$ESDM) & !is.null(data$Stack)){
        data$Stack <- project(data$Stack,data$Env)
        updateTabItems(session, "actions", selected="stack")
      }
    })


    ### Quit Menu ###

    ## quit page ##
    observeEvent(input$quitgui, {
      stopApp()
    })
  }
}
