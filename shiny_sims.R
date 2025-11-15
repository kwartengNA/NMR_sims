# libs ##############
# remove parameter grouping
library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(shinycssloaders)
# library(plotly)
# import data ###############
graph_res<-readRDS(file="shiny_results/graph_res.rds")
diag_res<-readRDS(file="shiny_results/diag_res.rds")
alldata2<-readRDS(file="shiny_results/input_data.rds")
N.sim<- max(as.numeric(diag_res$simulation))
gc()

# custom transform seed 545 #################
graph_res<-graph_res%>%
  mutate(density=factor(density,levels = levels(density),labels = paste0(levels(density), "%")),
         parameter=ifelse(parameter=="tau", "\u03C4", parameter))%>%#add to dataset
  rename(Source=src)

diag_res<-diag_res%>%
  mutate(density=factor(density,levels = levels(density),labels = paste0(levels(density), "%")),
         parameter=ifelse(parameter=="tau", "\u03C4", parameter))%>%#add to dataset
  rename(Source=src)

alldata2<-alldata2%>%
  mutate(density=factor(density,levels = levels(density),labels = paste0(levels(density), "%")))

# shiny part ##################

# scale transformation for Comparison Plots #############
stretch_trans<-scales::trans_new(
  name = "stretch_trans",
  transform = function(x) ifelse(x <= 0.8,
                                 x / 4,
                                 0.2 + 0.8 * (x - 0.8)/0.2),  # corrected forward
  inverse   = function(y) ifelse(y <= 0.2,
                                 y * 4,
                                 0.8 + (y - 0.2)/0.8 * 0.2)  # corrected inverse
)

# color palette for network  ################
direct_evidence_col<-data.frame(sumdir=1:4,
                                sumdircol=factor(1:4,levels=c(1:4),labels=c(  "#89CFF0","#627DEC","#281CA6","#090057")))#,"#090057"
# ui #########

ui <- dashboardPage(
  skin="blue",
  dashboardHeader(title = "NMR Simulation Summary Statistics",
                  # Insert tabs in header
                  tags$li(class = "dropdown",
                          tabsetPanel(
                            id = "main_tabs",
                            type = "tabs",
                            tabPanel("Coverage"),
                            tabPanel("Bias"),
                            tabPanel("MSE"),
                            tabPanel("CI Length"),
                            tabPanel("Diagnostics"),
                            tabPanel("Coverage Comparison"),
                            tabPanel("Bias Comparison"),
                            tabPanel("CI Length Comparison"),
                            tabPanel("Networks"),
                            tabPanel("Data")
                          )
                  )
  ),
  
  dashboardSidebar(
    conditionalPanel(
      condition = "input.main_tabs == 'Diagnostics'",
      selectInput(inputId="choose_diag", label="Choose Measure", 
                  choices=list("diagnostics"=unique(diag_res$parameter)), multiple=TRUE, selected="AIC")
    ),
    conditionalPanel(
      condition = "['Coverage', 'Bias', 'MSE', 'CI Length','Data','Coverage Comparison','Bias Comparison','CI Length Comparison'].includes(input.main_tabs)",
      selectInput(inputId="choose_coef", label="Choose Parameter", 
                  choices=list("coefficients"=unique(graph_res$parameter)), multiple=TRUE, selected="(Adjusted)\nTreatment Effect")
    ),
    conditionalPanel(
      condition = "['Data','Coverage Comparison','Bias Comparison','CI Length Comparison'].includes(input.main_tabs)",
      selectInput(inputId="choose_reference", label="Choose Reference Treatment", 
                  choices=list("coefficients"=unique(graph_res$reference_treat)), multiple=TRUE, selected=unique(graph_res$reference_treat))
    ),
    conditionalPanel(
      condition = "['Data','Coverage Comparison','Bias Comparison','CI Length Comparison'].includes(input.main_tabs)",
      selectInput(inputId="choose_comparator", label="Choose Comparator Treatment", 
                  choices=list("coefficients"=unique(graph_res$comparator_treat)), multiple=TRUE, selected=unique(graph_res$reference_treat))
    ),
    
    uiOutput('choose_arms'),
    uiOutput('choose_across_studies'),
    uiOutput('choose_consistency'),
    uiOutput('choose_tau'),
    uiOutput('choose_density'),
    uiOutput('choose_model'),
    uiOutput('choose_src'),
    uiOutput('choose_facets')
  ),
  
  dashboardBody(
    ##ui: aesthetics ###########
    tags$head(
      tags$style(HTML("
    /* Main header bar grey */
    .skin-blue .main-header .navbar {
      background-color: #4A4A4A;
    }

    /* Logo section on the left side grey */
    .skin-blue .main-header .logo {
      background-color: #4A4A4A;
      color: #ffffff !important;
      border-bottom: 0;
    }

    /* Logo hover: not an active button*/
    .skin-blue .main-header .logo:hover {
      background-color: #002878; 
    }
    
    /* Hover effect */
    .skin-blue .main-header .navbar .sidebar-toggle:hover {
      background-color: #002878;  /* hover background */
      color: #000000;              /* hover icon color */
    }
    
    /* Sidebar full background */
    .skin-blue .main-sidebar,
    .skin-blue .left-side {
      background-color: #002878;
    }

    /* Sidebar menu item text */
    .skin-blue .sidebar a {
      color: #ffffff;
    }

    /* Active menu item highlight */
    .skin-blue .sidebar .sidebar-menu > li.active > a {
      background-color: #FCB0B3;
      color: #ffffff;
    }

    /* Menu item hover */
    .skin-blue .sidebar .sidebar-menu > li > a:hover {
      background-color: #FCECC9;
      color: #ffffff;
    }
    
     /* Page background */
    .content-wrapper, .right-side, .main-body {
      background-color: #ffffff !important;
    }
    body, html {
      background-color: #ffffff !important;
    }
  "))
    ),
    ##ui: page content ##################
    
    
    conditionalPanel(
      condition = "input.main_tabs == 'Coverage'",
      shinycssloaders::withSpinner(imageOutput('coverageplot', height = "100vh"))
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Bias'",
      shinycssloaders::withSpinner(imageOutput('biasplot', height = "100vh"))
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'MSE'",
      shinycssloaders::withSpinner(imageOutput('mseplot', height = "100vh")),
      shinycssloaders::withSpinner(imageOutput('mseplottrend', height = "100vh"))
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'CI Length'",
      shinycssloaders::withSpinner(imageOutput('cilengthplot', height = "100vh"))
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Diagnostics'",
      shinycssloaders::withSpinner(imageOutput('aicplot', height = "80vh"))
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Coverage Comparison'",
      imageOutput('coverageComparison_es', height = "120vh"),
      imageOutput('coverageComparison_beta', height = "120vh")
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Bias Comparison'",
      imageOutput('biasComparison_es', height = "80vh"),
      imageOutput('biasComparison_beta', height = "80vh")
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'CI Length Comparison'",
      imageOutput('cilengthplotComparison_es', height = "80vh"),
      imageOutput('cilengthplotComparison_beta', height = "80vh")
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Networks'",
      shinycssloaders::withSpinner(imageOutput('networkplot', height = "80vh"))
    ),
    conditionalPanel(
      condition = "input.main_tabs == 'Data'",
      tableOutput('dataTable')
    )
  )
)

# server #########
server <- function(input, output, session) {
  # Create Filters ##############
  output$choose_arms<- renderUI({
    select_arms <- as.vector( unique(graph_res$true_arms) )
    radioButtons(inputId="true_arms",label="Choose Arms", choices=select_arms, selected="twoarm")
  })
  output$choose_across_studies <- renderUI({
    select_across_studies <- as.vector( unique(graph_res$true_assumption) )
    radioButtons(inputId="true_across_studies",label="Choose Across Study Assumption", choices=select_across_studies, selected="NMA")
  })
  output$choose_consistency <- renderUI({
    select_consistency <- as.vector( unique(graph_res$true_consistency) )
    radioButtons(inputId="true_consistency",label="Choose Consistency Assumptions", choices=select_consistency, selected=TRUE)
  })
  output$choose_tau <- renderUI({
    select_tau <- as.vector( unique(graph_res$tau_txt) )
    checkboxGroupInput(inputId="true_tau",label="Choose Heterogeneity", choices=select_tau, selected="moderate")
  })
  output$choose_density <- renderUI({
    select_density <- as.vector( unique(graph_res$density_txt) )
    checkboxGroupInput(inputId="true_density",label="Choose Density", choices=select_density, selected="full")
  })
  
  output$choose_model <- renderUI({
    select_model <- as.vector( unique(graph_res$implemented_model) )
    checkboxGroupInput(inputId="implemented_model",label="Choose Model", choices=select_model, selected="NMA")
  })
  
  output$choose_facets <- renderUI({
    select_facets <- c("tau", "density", "parameter", "none")
    selectInput(inputId="facet",label="Choose Faceting Variables", choices=select_facets, 
                multiple=TRUE, selected=c("tau", "density"))
  })
  
  output$choose_src <- renderUI({
    checkboxGroupInput(inputId="src",label="Choose Source", choices=unique(graph_res$Source), selected=c("Bayesian", "Frequentist"))
  })
  
  output$togglevalue_mse<-renderUI({
    radioButtons(inputId="togglevalue_mse",label="Show Value?", choices=c(TRUE,FALSE), selected=TRUE)
  })
  
  output$togglevalue_cov<-renderUI({
    radioButtons(inputId="togglevalue_cov",label="Show Value?", choices=c(TRUE,FALSE), selected=TRUE)
  })
  
  output$togglevalue_bias<-renderUI({
    radioButtons(inputId="togglevalue_bias",label="Show Value?", choices=c(TRUE,FALSE), selected=TRUE)
  })
  
  
  # Subset so that only the selected rows are in coef.data, unique_networks, and diag.data ####################
  unique_networks <- reactive({
    alldata2%>%
      filter(true_arms%in%c(input$true_arms))%>%
      filter(true_assumption%in%c(input$true_across_studies))%>%
      filter(true_consistency%in%c(input$true_consistency))%>%
      filter(tau_txt%in%c(input$true_tau))%>%
      filter(density_txt%in%c(input$true_density))%>%
      group_by(networkid)%>%# unique networks
      mutate(num_scenarios=length(unique(scenarionum)),
             which_scenarios=paste(unique(scenarionum), collapse=";"))%>%
      filter(scenario_sim==min(scenario_sim))%>%
      ungroup()
  })
  
  coef.data <- reactive({
    graph_res%>%
      filter(parameter%in%c(input$choose_coef))%>%
      filter(true_arms%in%c(input$true_arms))%>%
      filter(true_assumption%in%c(input$true_across_studies))%>%
      filter(true_consistency%in%c(input$true_consistency))%>%
      filter(tau_txt%in%c(input$true_tau))%>%
      filter(density_txt%in%c(input$true_density))%>%
      filter(implemented_model%in%input$implemented_model)%>%
      filter(Source%in%input$src)%>%
      filter(parameter%in%c(input$choose_coef))#%>%
      # filter(reference_treat %in%c(input$choose_reference))%>%
      # filter(comparator_treat %in%c(input$choose_comparator))
  })
  
  comparison.data <- reactive({
    coef.data()%>%
    filter(reference_treat %in%c(input$choose_reference))%>%
    filter(comparator_treat %in%c(input$choose_comparator))
  })
  
  diag.data <- reactive({
    diag_res%>%
      filter(parameter%in%c(input$choose_diag))%>%
      filter(true_arms%in%c(input$true_arms))%>%
      filter(true_assumption%in%c(input$true_across_studies))%>%
      filter(true_consistency%in%c(input$true_consistency))%>%
      filter(tau_txt%in%c(input$true_tau))%>%
      filter(density_txt%in%c(input$true_density))%>%
      filter(implemented_model%in%input$implemented_model)%>%
      filter(Source%in%input$src)
  })
  
  # outputs other ###########################
  output$networkplot <- renderPlot({
    
    n_facets <- length(unique(unique_networks()$networkid))+1
    facet_col <- ceiling(sqrt(n_facets))
    facet_row <- ceiling(n_facets / facet_col)
    par(mfrow=c(facet_row,facet_col),
        mar = c(1, 1, 1, 1))
    
    for(i in unique(unique_networks()$networkid)){
      tmp_pairwise<-meta::pairwise(treat=treatment,n = n,event = events,studlab=studlab,data=unique_networks()%>%filter(networkid==i),sm="OR")%>% # binary
        group_by(studlab)%>% # to get contributing direct comparisons in UMIE
        mutate(considered_direct=ifelse(arm_directionality1!=arm_directionality2, TRUE, FALSE))%>%
        ungroup()%>%
        group_by(treat1,treat2)%>%
        mutate(sumdir=sum(considered_direct))%>%
        ungroup()

      mynet<-netmeta::netmeta(TE = TE,seTE = seTE,treat1 = treat1,treat2 = treat2,studlab = studlab, data=tmp_pairwise, method.tau = "REML", reference.group = "A")
      if(nrow(tmp_pairwise%>%count(studlab)%>%filter(n>1))>0){
        netmeta::netgraph(mynet, plastic=FALSE, 
                          col=tmp_pairwise%>%select(treat1,treat2,sumdir)%>%distinct()%>%arrange(treat1, treat2)%>%mutate(sumdir=as.character(factor(sumdir, levels=direct_evidence_col$sumdir, labels=direct_evidence_col$sumdircol)))%>%pull(sumdir),
                          main=paste(unique(tmp_pairwise$num_scenarios), "Scenarios"), col.multiarm="#B8B8FF", multiarm=TRUE)
      } else {
        netmeta::netgraph(mynet, plastic=FALSE, 
                          col=tmp_pairwise%>%select(treat1,treat2,sumdir)%>%distinct()%>%arrange(treat1, treat2)%>%mutate(sumdir=as.character(factor(sumdir, levels=direct_evidence_col$sumdir, labels=direct_evidence_col$sumdircol)))%>%pull(sumdir),
                          main=paste(unique(tmp_pairwise$num_scenarios), "Scenarios"))
      }
    }
    legend
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), xaxt = "n", yaxt = "n", bty="n")
    title(main="Reference Containing\nComparisons (UMIE)", line=-1)
    legend(x="bottomleft",
           legend = direct_evidence_col$sumdir,
           # col = "black",
           fill = as.character(direct_evidence_col$sumdircol),
           # pch = 24,
           bty = "n",
           # pt.cex = 1,
           text.col = "black",
           cex=3,
           horiz = F ,
           inset = c(0.4,0))
    
  })

  output$dataTable <- renderTable({ coef.data() })
  
  output$coverageplot <- renderPlot({
    
    req(input$facet)
    req(nrow(coef.data()) > 0)
    facets<-paste(rev(input$facet), collapse="+")
    
    coverageplot<-
      coef.data()%>%
      filter(parameter %in%c("(Adjusted)\nTreatment Effect", "Treatment:Covariate\nInteraction Effect"))%>%#no tau
      select(Source,implemented_model,parameter,
             starts_with("mean"),prop_coverage,
             true_model,tau, density, true_arms)%>%
      distinct()%>%
      ggplot(aes(x=implemented_model, y=prop_coverage, fill=implemented_model, alpha=Source))+
      geom_bar(stat="identity",position=position_dodge(), color="black")+
      geom_hline(aes(yintercept = 0.95), lty=2, linewidth=1)+
      geom_hline(aes(yintercept = (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim))),lty=2,linewidth=1)+#ul
      geom_hline(aes(yintercept =(0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim))),lty=2, linewidth=1)+#ll
      scale_alpha_manual(values=c("Bayesian"=0.25,"Frequentist"=1))+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_y_continuous(name="Coverage Proportion",
                         trans=scales::trans_new(name="sq",
                                                 transform = function(x)x-0.95,
                                                 inverse=function(y) y+0.95),
                         breaks=c((0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)),
                                  (0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)),
                                  seq(0,1,0.05)),
                         labels = function(x) format(round(x, 2), nsmall = 2))+
      scale_x_discrete(name="Implemented Model")+
      labs(fill="Implemented Model")+
      facet_grid(as.formula(paste("parameter ~", facets)),scales = "free_y")+
      theme_bw()+
      theme(text = element_text(size = 16), axis.text.x = element_text(angle=90, vjust = 0), panel.spacing = unit(0, "lines"))
    
    print(coverageplot)
  })
  output$coverageComparison_es <- renderPlot({
    req(input$facet)
    req(nrow(comparison.data()) > 0)
    # print(unique(comparison.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    
    coverageComparison_es<-comparison.data()%>%
      filter(parameter%in%c("(Adjusted)\nTreatment Effect"))%>%
      mutate(modelsrc=paste(Source, implemented_model))%>%
      ggplot()+
      geom_hline(aes(yintercept =1),color = "black")+
      geom_hline(aes(yintercept =yint), data=data.frame(yint = seq((0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), #llimit, 
                                                                   (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), # ulimit, 
                                                                   0.0001)),color = "lightgrey")+# shaded confidence region
      geom_point(aes(x=comparison, y = coverage, color=implemented_model, shape = Source))+
      geom_line(aes(x=comparison, y = coverage, group=modelsrc, color=implemented_model, lty=Source), lwd=1)+
      scale_color_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_linetype_manual(values=c("Bayesian"=2, "Frequentist"=1))+
      scale_y_continuous(name = "Coverage Proportion", trans=stretch_trans, limits=c(0,1), breaks = seq(0.8,1, 0.05))+
      facet_grid(density~tau+parameter)+
      coord_polar()+
      theme_bw()
    print(coverageComparison_es)
    
  })
  output$coverageComparison_beta <- renderPlot({
    req(input$facet)
    req(nrow(comparison.data()) > 0)
    # print(unique(comparison.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    coverageComparison_beta<-comparison.data()%>%
      filter(parameter%in%c("Treatment:Covariate\nInteraction Effect"))%>%
      mutate(modelsrc=paste(Source, implemented_model))%>%
      ggplot()+
      # geom_hline(aes(yintercept =ulimit), size=1,color = "lightgrey")+# non shaded
      # geom_hline(aes(yintercept =llimit), color = "black")+
      geom_hline(aes(yintercept =1),color = "black")+
      geom_hline(aes(yintercept =yint), data=data.frame(yint = seq((0.95-qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), #llimit, 
                                                                   (0.95+qnorm(1-0.05/2)*sqrt((0.95*(1-0.95))/N.sim)), # ulimit, 
                                                                   0.0001)),color = "lightgrey")+# shaded confidence region
      geom_point(aes(x=comparison, y = coverage, color=implemented_model, shape = Source))+
      geom_line(aes(x=comparison, y = coverage, group=modelsrc, color=implemented_model, lty=Source), lwd=1)+
      scale_y_continuous(name = "Coverage Proportion",trans=stretch_trans, limits=c(0,1), breaks = seq(0.8,1, 0.05))+
      scale_color_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_linetype_manual(values=c("Bayesian"=2, "Frequentist"=1))+
      facet_grid(density~tau+parameter)+
      # facet_wrap(as.formula(paste("parameter ~", facets)))+
      coord_polar()+
      theme_bw()
    
    print(coverageComparison_beta)
    
  })
  
  output$biasComparison_es <- renderPlot({
    req(input$facet)
    req(nrow(comparison.data()) > 0)
    # print(unique(comparison.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    
    biasComparison_es<-comparison.data()%>%
      filter(parameter=="(Adjusted)\nTreatment Effect")%>%
      mutate(modelsrc=paste(Source, implemented_model))%>%
      ggplot()+
      geom_col(aes(x=comparison, y = bias, group=modelsrc, fill=implemented_model, alpha=Source), lwd=1, position = position_dodge())+
      geom_hline(aes(yintercept =0),color = "red",lwd=1)+
      facet_grid(density~tau+parameter)+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_alpha_manual(values=c("Bayesian"=0.5, "Frequentist"=1))+
      scale_y_continuous(name  = "Mean Bias", limits=c(-1*max(abs(comparison.data()$bias[comparison.data()$parameter=="(Adjusted)\nTreatment Effect"])),max(abs(comparison.data()$bias[comparison.data()$parameter=="(Adjusted)\nTreatment Effect"]))))+
      theme_bw()
    
    print(biasComparison_es)
    
  })
  output$biasComparison_beta <- renderPlot({
    req(input$facet)
    req(nrow(comparison.data()) > 0)
    # print(unique(comparison.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    biasComparison_beta<-comparison.data()%>%
      filter(parameter=="Treatment:Covariate\nInteraction Effect")%>%
      mutate(modelsrc=paste(Source, implemented_model))%>%
      ggplot()+
      geom_col(aes(x=comparison, y = bias, group=modelsrc, fill=implemented_model,  alpha=Source), lwd=1, position = position_dodge())+
      geom_hline(aes(yintercept =0),color = "red",lwd=1)+
      facet_grid(density~tau+parameter)+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_alpha_manual(values=c("Bayesian"=0.5, "Frequentist"=1))+
      scale_y_continuous(name  = "Mean Bias", limits=c(-1*max(abs(comparison.data()$bias[comparison.data()$parameter=="Treatment:Covariate\nInteraction Effect"])),max(abs(comparison.data()$bias[comparison.data()$parameter=="Treatment:Covariate\nInteraction Effect"]))))+
      theme_bw()
    
    print(biasComparison_beta)
    
  })
  output$biasplot <- renderPlot({
    req(input$facet)
    req(nrow(coef.data()) > 0)
    # print(unique(coef.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    
    biasplot<-coef.data()%>%
      select(Source,implemented_model,parameter,
             starts_with("mean"),prop_coverage,
             true_model,tau, density, true_arms)%>%
      distinct()%>%
      ggplot(aes(x = implemented_model, y = mean_bias, fill = implemented_model, alpha=Source)) +
      geom_bar(stat="identity",position=position_dodge(), color="black")+
      geom_hline(aes(yintercept = 0))+
      facet_grid(as.formula(paste("parameter ~", facets)),scales = "free_y")+ #prompt for user input
      scale_alpha_manual(values=c("Bayesian"=0.25,"Frequentist"=1))+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      theme_bw()+
      ylab(label = "Mean Bias")+
      scale_x_discrete(name="Implemented Model")+
      labs(fill="Implemented Model")+
      theme(text = element_text(size = 16), axis.text.x = element_text(angle=90, vjust = 0), panel.spacing = unit(0, "lines"))
    
    print(biasplot)
  })
  

  output$cilengthplot <- renderPlot({
    req(input$facet)
    req(nrow(coef.data()) > 0)
    # print(unique(coef.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    
    cilengthplot<-
      coef.data()%>%
      select(Source,implemented_model,parameter,
             starts_with("mean"),prop_coverage,
             true_model,tau, density, true_arms)%>%
      distinct()%>%
      ggplot(aes(x = implemented_model, y = mean_CI_length, fill = implemented_model, alpha=Source)) +
      geom_bar(stat="identity",position=position_dodge(), color="black")+
      geom_hline(aes(yintercept = 0))+
      # facet_grid(parameter~tau+density,scales = "free_y")+ #manual debugging
      facet_grid(as.formula(paste("parameter ~", facets)),scales = "free_y")+ #prompt for user input
      # scale_fill_manual(values=c("Bayesian"="#FAA613","Frequentist"="#000080"))+
      scale_alpha_manual(values=c("Bayesian"=0.25,"Frequentist"=1))+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_x_discrete(name="Implemented Model")+
      ylab(label = "Mean Confidence Interval Length")+
      labs(fill="Implemented Model")+
      theme_bw()+
      theme(text = element_text(size = 16), axis.text.x = element_text(angle=90, vjust = 0), panel.spacing = unit(0, "lines"))
    
    print(cilengthplot)
  })
  output$cilengthplotComparison_es <- renderPlot({
    req(input$facet)
    req(nrow(comparison.data()) > 0)
    # print(unique(comparison.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    
    cilengthplotComparison_es<-comparison.data()%>%
      filter(parameter%in%c("(Adjusted)\nTreatment Effect"))%>%
      mutate(modelsrc=paste(Source, implemented_model))%>%
      ggplot()+
      geom_col(aes(x=comparison, y = CI_length, fill=implemented_model, alpha=Source), position = position_dodge())+
      scale_alpha_manual(values=c("Bayesian"=0.5, "Frequentist"=1))+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      ylab(label = "Mean Confidence Interval Length")+
      facet_grid(density~tau+parameter)+
      theme_bw()
    
      print(cilengthplotComparison_es)
    
  })
  
  output$cilengthplotComparison_beta <- renderPlot({
    req(input$facet)
    req(nrow(comparison.data()) > 0)
    facets<-paste(rev(input$facet), collapse="+")
    
    cilengthplotComparison_beta<-comparison.data()%>%
      filter(parameter%in%c("Treatment:Covariate\nInteraction Effect"))%>%
      mutate(modelsrc=paste(Source, implemented_model))%>%
      ggplot()+
      geom_col(aes(x=comparison, y = CI_length, fill=implemented_model, alpha=Source), position = position_dodge())+
      scale_alpha_manual(values=c("Bayesian"=0.5, "Frequentist"=1))+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      ylab(label = "Mean Confidence Interval Length")+
      facet_grid(density~tau+parameter)+
      theme_bw()
    print(cilengthplotComparison_beta)
  })
  output$mseplottrend <- renderPlot({
    req(input$facet)
    req(nrow(coef.data()) > 0)
    facets<-paste(rev(input$facet), collapse="+")
    
    
    mseplottrend<-coef.data()%>%
      # filter(true_model=="NMA" & true_arms=="twoarm" & parameter=="(Adjusted)\nTreatment Effect")%>%
      select(Source,implemented_model,parameter,
             starts_with("mean"),prop_coverage,
             true_model,tau, density, true_arms)%>%
      mutate(modelsrc=paste(implemented_model,Source))%>%
      distinct()%>%
      ggplot(aes(x = interaction(density,tau, sep = "\n"), y = mean_mse, color = implemented_model, group = modelsrc, alpha=Source)) +
      geom_line(linewidth = 1) +
      geom_hline(aes(yintercept = 0))+
      geom_point(size = 2) +
      xlab(label = "Scenario")+
      ylab(label = "Mean MSE")+
      geom_point()+
      geom_hline(aes(yintercept = 0))+
      facet_grid(as.formula(paste("parameter ~", facets)),scales = "free_y")+
      scale_alpha_manual(values=c("Bayesian"=0.25,"Frequentist"=1))+
      scale_color_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      scale_x_discrete(name="Implemented Model")+
      labs(fill="Implemented Model")+
      theme_bw()+
      theme(text = element_text(size = 16), axis.text.x = element_text(angle=90, vjust = 0), panel.spacing = unit(0, "lines"))
    
    print(mseplottrend)
  })
  output$mseplot <- renderPlot({
    req(input$facet)
    req(nrow(coef.data()) > 0)
    facets<-paste(rev(input$facet), collapse="+")
    
    mseplot<-coef.data()%>%
      select(Source,implemented_model,parameter,
             starts_with("mean"),prop_coverage,
             true_model,tau, density, true_arms)%>%
      distinct()%>%
      ggplot(aes(x = implemented_model, y=mean_mse,
                 fill = implemented_model, alpha=Source)) +
      geom_bar(stat="identity",position=position_dodge(), color="black")+
      geom_hline(aes(yintercept = 0))+
      facet_grid(parameter~tau+density,scales = "free_y")+
      scale_alpha_manual(values=c("Bayesian"=0.25,"Frequentist"=1))+
      scale_fill_manual(values=c("Common"="#000080","Common\nUMIE"="#1C77C3", "Independent"="#FF33CC", "Independent\nUMIE"="#D7263D","Independent\nUMIE"="#D7263D", "NMA"="#FBCEB1"))+
      theme_bw()+
      theme(text = element_text(size = 16), axis.text.x = element_text(angle=90, vjust = 0), panel.spacing = unit(0, "lines"))
    
    print(mseplot)
  })
  
  

  
  output$aicplot <- renderPlot({
    req(input$facet)
    req(nrow(diag.data()) > 0)
    # # print(unique(diag.data()$parameter))
    facets<-paste(rev(input$facet), collapse="+")
    
    aic_res<-
      diag.data()%>%
      filter(parameter %in%c("AIC", "AICc", "BIC", "dev", "ll"))%>%
      select(fileloc, simulation,true_model, implemented_model, parameter, ML, REML, true_assumption, true_consistency)%>%
      mutate(model_filter=stringr::str_trim(gsub("\n"," ", implemented_model)))%>%
      filter(true_model==model_filter)%>%
      select(-c(true_consistency,true_assumption))%>%
      left_join(diag.data(), by=c("fileloc", "simulation", "parameter"), suffix=c("_truemodel", ""))%>%
      group_by(fileloc, simulation)%>%
      arrange(fileloc, simulation, (REML))%>%
      mutate(diff_truemodel_REML=REML-REML_truemodel,
             rank=row_number())%>%
      ungroup()
    
    
    aicplot<-
      aic_res%>%
      select(parameter, true_assumption, true_consistency,true_model, tau, 
             density,density_txt,
             implemented_model,diff_truemodel_REML, REML, REML_truemodel, rank)%>%
      group_by(parameter, tau, density,density_txt, true_model,
               implemented_model)%>%
      mutate(mean_diff_NMA_REML=mean(diff_truemodel_REML, na.rm=TRUE), median_rank=median(rank))%>%
      arrange(tau, density,diff_truemodel_REML )%>%
      ggplot()+
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2,fill = "grey", alpha = 0.5) +
      geom_boxplot(aes(y=diff_truemodel_REML, x=interaction(density, tau, sep="\n"), fill=implemented_model), alpha=0.85)+
      scale_fill_manual(values=c("NMA"="#FAA613","Common"="#000080",  "Independent"="#D7263D","Independent\nUMIE"="#F57788","Common\nUMIE"="#1C77C3"))+
      facet_grid(true_model~parameter,scales = "free_y", switch="y")+ #prompt for user input
      theme_bw()+
      theme(text = element_text(size = 16), axis.text.x = element_text(angle=0, vjust = 0), panel.spacing = unit(0, "lines"))
    
    print(aicplot)
    
  })
  
gc()
  }

# run shiny ###################
shinyApp(ui, server)

# END SHINY ####################



