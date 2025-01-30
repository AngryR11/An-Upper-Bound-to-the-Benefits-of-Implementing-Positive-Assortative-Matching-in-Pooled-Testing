rm(list=ls())





library(poisbinom)
library(shiny)

ui <- fluidPage(
  
  # CSS to remove number input arrows
  tags$head(tags$style(HTML("
    /* Remove spinners from number inputs */
    input[type=number]::-webkit-inner-spin-button, 
    input[type=number]::-webkit-outer-spin-button { 
      -webkit-appearance: none; 
      margin: 0; 
    }
    input[type=number] { 
      -moz-appearance: textfield;
    }
  "))),
  
  titlePanel("Maximum benefits of matching together those with similar probabilities of infeciton when implementing Dorfman Testing"),
 
  # Adding reference to the working paper
  tags$p("Please cite the following paper when referencing this application: ",
         tags$a(href = "https://ssrn.com/abstract=4779050",
                "An Upper Bound to the Benefits of Implementing Positive Assortative Matching in Pooled Testing by Gustavo Quinderé Saraiva.", target = "_blank"), "This paper provides a detailed explanation of the methodology used to compute the upper bounds."),
  
   sidebarLayout(
    sidebarPanel(
      # Add a link to the documentation
      #tags$h4("Documentation"),
      tags$a("View Documentation", href = "README.html", target = "_blank"),
      # Step 1: Initial Inputs for Dilution Function
      #h4("Step 1: Specify the Dilution Function"),
      h4("Step 1: Specify the Dilution Function"),
      tags$p("The dilution function is defined as:"),
      withMathJax(
        tags$p("$$(1 - Sp) + (Sp+Se- 1)\\left(\\frac{I}{K}\\right)^\\delta$$"),
        tags$p("This represents the probability that infection is detected, conditional on 'I' of the 'K' specimens in the pool being infected.")
      ),
      tags$div(
        class = "form-group shiny-input-container",
        tags$label("Sensitivity of the test (Se)"),
        tags$input(
          id = "input7", 
          type = "number", 
          class = "form-control", 
          placeholder = "e.g., 0.9",
          min = 0, 
          step = 0.01,
          max = 1
        )
      ),
      tags$div(
        class = "form-group shiny-input-container",
        tags$label("Specificity of the test (Sp)"),
        tags$input(
          id = "input8", 
          type = "number", 
          class = "form-control", 
          placeholder = "e.g., 0.95",
          min = 0, 
          step = 0.01,
          max = 1
        )
      ),
      tags$div(
        class = "form-group shiny-input-container",
        tags$label(HTML("Dilution Parameter &delta;: higher values indicate greater dilution effects")),
        tags$input(
          id = "input4", 
          type = "number", 
          class = "form-control", 
          placeholder = "e.g., 0.1",
          min = 0, 
          step = 0.01
        )
      ),
      actionButton("plot_dilution", "Display Dilution Function"),
      
      # Step 2: Conditional Inputs
      conditionalPanel(
        condition = "output.dilutionReady",
        h4("Step 2: Proceed with Additional Inputs"),
        tags$div(
          class = "form-group shiny-input-container",
          tags$label("Prevalence:"),
          tags$input(
            id = "input1", 
            type = "number", 
            class = "form-control", 
            placeholder = "e.g., 0.05",
            min = 0, 
            max = 1,
            step = 0.01
          )
        ),
        tags$div(
          class = "form-group shiny-input-container",
          tags$label("Max. Prob. Infection:"),
          tags$input(
            id = "input2", 
            type = "number", 
            class = "form-control", 
            placeholder = "e.g., 0.1",
            min = 0, 
            max = 1,
            step = 0.01
          )
        ),
        tags$div(
          class = "form-group shiny-input-container",
          tags$label("Min. Prob. Infection:"),
          tags$input(
            id = "input3", 
            type = "number", 
            class = "form-control", 
            placeholder = "e.g., 0.01",
            min = 0, 
            max = 1,
            step = 0.01
          )
        ),
        tags$div(
          class = "form-group shiny-input-container",
          tags$label("Pool size"),
          tags$input(
            id = "input5", 
            type = "number", 
            class = "form-control", 
            placeholder = "e.g., 5",
            min = 2
          )
        ),
        tags$div(
          class = "form-group shiny-input-container",
          tags$label("Number of pools"),
          tags$input(
            id = "input6", 
            type = "number", 
            class = "form-control", 
            placeholder = "e.g., 20"
          )
        ),
        tags$div(
          class = "form-group shiny-input-container",
          tags$label("Cost of a test (USD):"),
          tags$input(
            id = "input9", 
            type = "number", 
            class = "form-control", 
            placeholder = "e.g., 25"
          )
        ),
        actionButton("calculate", "Calculate")
      )
    ),
    mainPanel(
      h4("Plot Dilution Function"),
      plotOutput("dilution_plot"),
      h4("Maximum benefits obtained by matching together those with similar probabilities of infection"),
      uiOutput("result")
    )
  )
)


server <- function(input, output, session) {
  # Reactive value to track if dilution is ready
  dilutionReady <- reactiveVal(FALSE)
  
 
  # Reactive expression to store the dilution function
  dilution_reactive <- reactive({
    req(input$input7, input$input8, input$input4) # Ensure inputs are not NULL
    Se <- as.numeric(input$input7)
    Sp <- as.numeric(input$input8)
    delta <- as.numeric(input$input4)
    
    function(I, k) {
      (1 - Sp) + (Sp - 1 + Se) * (I / k)^delta
    }
  })
  
  # Observe the button for the first step (Dilution Plot)
  observeEvent(input$plot_dilution, {
    Se <- as.numeric(input$input7)
    Sp <- as.numeric(input$input8)
    delta <- as.numeric(input$input4)
    
    if (is.na(Se) || is.na(Sp) || is.na(delta)) {
      showNotification("Please fill in all fields for Sensitivity, Specificity, and Delta", type = "error")
    } else if (Se < 0 || Se > 1 || Sp < 0 || Sp > 1 || delta < 0|| delta>1) {
      showNotification(HTML("Please ensure Se, Sp and &delta; lie between 0 and 1"), type = "error")
    } else if(Se+Sp<1){
      showNotification("Please ensure that Se+Sp>1", type = "error")
    }else{
      # Enable step 2
      dilutionReady(TRUE)
      
      # Render the dilution plot
      output$dilution_plot <- renderPlot({
        dilution <- function(I, k) {
          (1 - Sp) + (Sp - 1 + Se) * (I / k)^delta
        }
        
        k_max <- 100
        i_grid <- seq(0, k_max, by = 0.001)
        y <- sapply(i_grid, function(I) dilution(I, k_max))
        
        plot(
          i_grid / k_max, y,
          type = "l",
          xlab = "Proportion of Infected within the pool, I/K",
          ylab = "Probability infection is detected",
          cex.lab = 1.5,
          lwd = 2,
          main = "Dilution Function"
        )
      })
    }
  })
  
  # Calculate the final result
  observeEvent(input$calculate, {
    # Access the remaining inputs
    p <- as.numeric(input$input1)#prevalence
    b <- as.numeric(input$input2)#max probability of infection
    a <- as.numeric(input$input3)#minimum prob. of infection
    k <- as.numeric(input$input5)#pool size
    m <- as.numeric(input$input6)#number of pools
    n=m*k#batch size
    C <- as.numeric(input$input9)#cost of a test
    
    if (is.na(p) || is.na(b) || is.na(a) || is.na(k) || is.na(C)) {
      showNotification("Please fill in all fields for Step 2", type = "error")
    }else if(k>60){
      showNotification("Please select a smaller pool size", type = "error")
    } else if(a>p){
      showNotification("Please select a value lower than the prevalence for the minimum probability of infection", type = "error")
    }
    else if(b<p){
      showNotification("Please select a value higher than the prevalence for the maximum probability of infection", type = "error")
    }
    else if((a==p || b==p) && (a<p || b>p)){
      showNotification("The minimum and maximum probabilities of infection must be consistent with the prevalence", type = "error")
    }
    else if(n>1000 || is.na(m)){#if the batch size is too large, or the user did not provide the number of pools, use the asymptotic upper bound
      source('function_tight_upper_bound.R')
      out <- asymptotic_UB_ET(a,b,p,k,dilution=dilution_reactive())
      out_fp<-asymptotic_UB_FP(a,b,p,k,dilution=dilution_reactive())
      out_fn<-asymptotic_UB_FN(a,b,p,k,dilution=dilution_reactive())
      
      print(class(out))  # Check if it's numeric
      if (!is.numeric(out)) {
        stop("Error: Output 'out' is not numeric.")
      }
      
      # Calculate the additional result (e.g., cost adjustment)
      cost_adjusted_result <- out * input$input9  # Multiply result by cost (input$input9)
      
  
      output$result <- renderUI({
        HTML(paste0(
          "<strong style='font-size:16px;'>Maximum reduction in the expected number of tests per subject:</strong> ", 
          "<span style='color:blue; font-size:18px; font-weight:bold;'>", round(out, 4), "</span><br><br>",
          
          "<strong style='font-size:16px;'>Maximum cost reduction per subject:</strong> ", 
          "<span style='color:green; font-size:18px; font-weight:bold;'>$", round(cost_adjusted_result, 4), "</span><br><br>",
          
          "<strong style='font-size:16px;'>Maximum reduction in false positives per subject:</strong> ", 
          "<span style='color:red; font-size:18px; font-weight:bold;'>", round(out_fp, 5), "</span><br><br>",
          
          "<strong style='font-size:16px;'>Maximum reduction in false negatives per subject:</strong> ", 
          "<span style='color:purple; font-size:18px; font-weight:bold;'>", round(out_fn, 5), "</span><br><br>"
        ))
      })
      
      
      
      
    }
    else {
      source('function_tight_upper_bound.R')
      q=max_var_q_2(a,b,n,p,k)
      lambda=.5#only used for equity analysis
      out <- Random_ET_FN_FP_EQT(p,k,n,dilution=dilution_reactive())[1]-Ordered_ET_FN_FP_EQT(q,k,lambda,dilution=dilution_reactive())[1]
      out_fp<-Random_ET_FN_FP_EQT(p,k,n,dilution=dilution_reactive())[2]-Ordered_ET_FN_FP_EQT(q,k,lambda,dilution=dilution_reactive())[2]
      out_fn<-Random_ET_FN_FP_EQT(p,k,n,dilution=dilution_reactive())[3]-Ordered_ET_FN_FP_EQT(q,k,lambda,dilution=dilution_reactive())[3]
      
      print(class(out))  # Check if it's numeric
      if (!is.numeric(out)) {
        stop("Error: Output 'out' is not numeric.")
      }
      
      # Calculate the additional result (e.g., cost adjustment)
      cost_adjusted_result <- out * input$input9  # Multiply result by cost (input$input9)
      

    

      output$result <- renderUI({
        HTML(paste0(
          "<strong style='font-size:16px;'>Maximum reduction in the expected number of tests per subject:</strong> ", 
          "<span style='color:blue; font-size:18px; font-weight:bold;'>", round(out,5), "</span><br><br>",
          
          "<strong style='font-size:16px;'>Maximum cost reduction per subject:</strong> ", 
          "<span style='color:green; font-size:18px; font-weight:bold;'>$", round(cost_adjusted_result, 4), "</span><br><br>",
          
          "<strong style='font-size:16px;'>Maximum reduction in false positives per subject:</strong> ", 
          "<span style='color:red; font-size:18px; font-weight:bold;'>", round(out_fp, 5), "</span><br><br>",
          
          "<strong style='font-size:16px;'>Maximum reduction in false negatives per subject:</strong> ", 
          "<span style='color:purple; font-size:18px; font-weight:bold;'>", round(out_fn, 5), "</span><br><br>"
        ))
      })      
      #  output$result <- renderUI({
      #   HTML(paste(out_1,out_2,out_3,out_4,sep='<br>'))
      # })
    
    }
  })
  
  # Reactive expression for conditionalPanel
  output$dilutionReady <- reactive({
    dilutionReady()
  })
  outputOptions(output, "dilutionReady", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)

