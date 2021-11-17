Sys.setenv(RETICULATE_PYTHON = ".epmv/bin/python3")
#reticulate::use_virtualenv(".epmv", required = TRUE) 
library(reticulate)
library(ggplot2)
library(DT)

#options(shiny.trace=TRUE)
#reticulate::virtualenv_create("epmv")
#reticulate::virtualenv_install("epmv", packages = c("EpigeneticPacemaker", "matplotlib","numpy"))
#reticulate::use_virtualenv("epmv", required = TRUE)
#reticulate::use_virtualenv('C:/Users/XZzzZzz/Desktop/EPM/epmvenv', required = TRUE)
#library(car)
#conda_create("r-scrublet")
#conda_install(envname="r-scrublet", packages ="numpy","pip","git")
#conda_python(envname =  "r-scrublet")
#reticulate::conda_create("my-environment")
#reticulate::use_condaenv("epmvenv", required = TRUE)

#py_install('matplotlib',pip=TRUE)
#py_module_available('numpy')

source_python("data_get.py")
#source_python("epi_function.py")
source_python("epi_fortest.py")


### SHINY UI ###
ui <- fluidPage(
  # Application title
  navbarPage(title="Epigenetic Pacemaker Model",id="mainpage",
    tabPanel("Home",value="tabhome",
       h1("Welcome to Epigenetic Pacemaker",
       style="color:grey;front-size:35px;text-align:center;"),
       br(),
       fluidRow(
         column(10,offset=1,
         fluidRow(
           img(src="mainpage_fig.png",width="20%"),align="center")
         )
         #style = "margin-bottom:20%; margin-top:20%"
       ),
       br(),
       p("DNA methylation is widely used to model physiological phenotypes, 
          such as aging and type II diabetes. The Epigenetic Pacemaker, EPM, 
          is an implementation of a fast conditional expectation maximization algorithm that models epigenetic 
          states under and evolutionary framework. The EPM was first introduced by Snir et al. 
          as an extension of the Universal Pacemaker (UPM) model of genome evolution. 
          In contrast to regression bases approaches, the EPM does not assume a linear relationship 
          between the epigenetic state and a trait of interest. As a result the EPM can model non-linear 
          epigenetic trait associations directly without transformation of the phenotype of interest.",
          style = "text-align:justify;font-family:'times'; font-si16pt; width:80%; margin-left:10%; margin-right:10%", 
          align="center"),
       p("Reference: https://epigeneticpacemaker.readthedocs.io/en/latest/",
          style = "font-family:'times'; font-si16pt; width:80%; margin-left:10%; margin-right:10%")
       ),
    navbarMenu("Datasets",
          
       # tabPanel("Example Dataset",value="example_dataset",
       #    sidebarLayout(
       #      sidebarPanel(
       #        fileInput("file",label=h3("Choose CSV file"),
       #                  multiple = TRUE,
       #                  accept = c("text/csv", ".csv")),
       #        tags$hr(),
       #        fluidRow(column(4, verbatimTextOutput("value")))
       #      ),
       #      mainPanel(
       #        tabsetPanel(
       #          tabPanel("Simple",plotOutput("prediction"))
       #        )
       #      )
       #    )
       # ),
      
       tabPanel("Sample File",value="sample_file",
                sidebarPanel(h3("Before uploading your file..."),
                             hr(),
                             h5("NOTE: This sample on the right shows the format of the file you should 
                                upload to the test_dataset part (The first line with Col1, Col2... 
                                is not required in your file. For a correct output,
                                your file should have exactly the same format
                                starting from the row with Sample IDs.",
                                style = "text-align:justify;font-family:'times'")),
                mainPanel(dataTableOutput("sample_table"))
       ),
       
			 tabPanel("Test Dataset",value="test_dataset",
          sidebarLayout(
             sidebarPanel(
                tabsetPanel(id='tabset',
                   tabPanel('YOUR TURN NOW!!!',
                      fileInput("test_file",label=h3("File Input")),
                      
                      textInput('num_input',
                                h3('Define your PCC value'),
                                value="Enter text..."),
                      h6('Note: Your PCC value should be between 0~1 . 
                         If your value is too high (i.e. >0.9) then sometimes there might be an error
                         when generating the graph because there is no site with such high PCC value.',
                         style = "text-align:justify;font-family:'times'")
                   )
                ),
                actionButton('go','Plot'),
                hr(),
                h3("Download the results of the Epigenetic Pacemaker"),
                downloadButton("downloadData","Download")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel('Plot',plotOutput("test_plot")),
                 tabPanel('Summary',verbatimTextOutput("selected_sites")),
                 tabPanel('Table',dataTableOutput("table"))
               )
             )
          ) 
       )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$prediction <- renderPlot({
    plot(test_ages, test_predict, main="EpigeneticPacemaker", 
         xlab="Chronological Age",ylab="Epigenetic State",pch=19)
    abline(lm(log(test_predict)~log(test_ages)),col="red")
  })
  
  output$sample_table<-renderDataTable({
    df <- rbind(c("Sample IDs","ID1","ID2","ID3","..."),c("Age","Age1","Age2","Age3","..."),
                c("site1","value1","value2","value3","..."),c("site2","value1","value2","value3","..."),
                c("site3","value1","value2","value3","..."),c("site...","...","...","...","..."))
    colnames(df) <- c('Col1','Col2',"Col3","Col4","...")
    df
  })
  
  pcc_val <- reactiveValues(doPlot=FALSE)
  observeEvent(input$go,{
    pcc_val$doPlot <- input$go
  })

  output$test_plot <- renderImage({
    if(pcc_val$doPlot==FALSE)
      return()
    isolate({
      req(input$test_file)
      tryCatch({ 
        filename <- (input$test_file)$name
        print(filename)
        pcc <- as.numeric(input$num_input)
        #print(typeof(pcc),pcc)
      },
      error = function(e) {
        stop(safeError(e))
      })
      r <- plot_testdata(filename,pcc)
      print('***********************')
      #print(r[[3]])
      list(src='myplot.png',contentType='image/png',width="100%")
    })
  },
  deleteFile=TRUE
  )
  
  output$selected_sites <- renderText({
    if(pcc_val$doPlot==FALSE)
      return()
    isolate({
      filename <- (input$test_file)$name
      pcc <- as.numeric(input$num_input)
      r <- plot_testdata(filename,pcc)
      paste("No. of sites you have selected:",r[[1]],
            "No. of individuals:",r[[2]], sep='\n')
    })
  })
  
  output$table <- renderDataTable({
    filename <- (input$test_file)$name
    pcc <- as.numeric(input$num_input)
    r <- plot_testdata(filename,pcc)
    df <- data.frame(r[[3]],r[[4]])
    colnames(df) <- c('Chronological Ages','Epigenetic States')
    df
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$test_file$name, ".csv", sep = "")
    },
    content = function(file) {
      filename <- (input$test_file)$name
      pcc <- as.numeric(input$num_input)
      r <- plot_testdata(filename,pcc)
      df <- data.frame(r[[3]],r[[4]])
      colnames(df) <- c('Chronological Ages','Epigenetic States')
      write.csv(df, file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
