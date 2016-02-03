if (!requireNamespace("shinyFiles", quietly = TRUE)) {
  stop("shinyFiles package is needed for the Global User Interface. Please install it.")
}
library(shinydashboard)
library(shinyFiles)
ui <- dashboardPage(dashboardHeader(title = 'SSDM'),
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
                                        shinyFilesButton('envfiles', 'Raster selection', 'Please select rasters', FALSE, multiple = T),
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
                                             uiOutput('endemismrangeUI'),
                                             uiOutput('endemismrangeinfoUI'),
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
                                        shinyDirButton('save', 'Folder selection', 'Please select folder to save the model', FALSE),
                                        actionButton('savemodel', 'save', icon = icon('floppy-o'))
                                    )
                                  )
                                )
                        )
                      )
                    )
)
