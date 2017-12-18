
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library("shiny")
shinyUI(fluidPage(
  # Header ---------------------------------------------
  titlePanel("Defensins Sequence Space"),
  wellPanel(
    div(class="header",
        p("There are",
          a(href="https://doi.org/10.1093/molbev/msw106","two evolutionarily independent superfamilies of defensins"),
          "(1) with",
          a(href="https://doi.org/10.1007/s00018-016-2344-5","extensive convergent similarities"),
          "(2) . The cis-defensins contain plant, fungal, and most invertebrate defensins. The trans-defensins contain
          the mammalian defensins and invertebrate big defensins."),
        p("This webtool allows you to interactively view the sequence space maps for these defensin superfamilies.
          It can predict whether a given query sequence is a cis-defensin or trans-defensin, map
          its sequence space location, and determine which cluster it falls into.")
        )
    ),

  tabsetPanel(
  # TAB 1 -----------------------------------
    tabPanel("View",fluidPage(
      sidebarPanel(
        h3("Sequence Space"),
        p("View a quantitative sequence space map for either defensin superfamily."),
        selectInput("view_type",
                    label   = "Defensin superfamily",
                    choices  = c("cis-Defensin","trans-Defensin")),
        p("Clicking \"show\" will open the interactive diagram in a new window."),
        actionButton("button.view", "Show"),
        tags$hr(style="border-color: grey;"),
        h4("Identify sequences"),
        p("To identify points in the sequence space, use the button below,
          then drag a box around the sequences you are interested in."),
        actionButton("button2.view", "Activate selection mode")
      ),
      
      # Image legend
      # mainPanel(
      #   textOutput("report_legend")
      # ),
      # Display selected points alignment
      mainPanel(
        h4("Selected sequences"),
        htmlOutput("fasta2.view",style = "font-family:Courier New ; min-width:500px")
      )
      
    )),

      
    # TAB 2 -----------------------------------
    tabPanel("Find",fluidPage(
  
    # Inputs section
    sidebarLayout(
  
      # Select superfamily
      sidebarPanel(
        h3("Locate a sequence"),
        p("Input a query sequence to locate in the defensin sequence space. Since the defensins consist
          of two superfamilies, you may indicate whether it is a cis-defensin of trans-defensin (if known).
          Accuracy is best when the query sequence includes only the defensin domain."),
        textAreaInput("query_sequence",resize="vertical",rows=3,
                  label       = "Find query sequence",
                  value       = "RECKTESNTFPGICITKPPCRKACISEKFTDGHCSKILRRCLCTKPC",
                  placeholder = "Amino acid sequence"),
        selectInput("query_type",
                    label   = "Within superfamily (if known)",
                    choices  = c("unknown","cis-Defensin","trans-Defensin"),
                    selected = "unknown"),
        p("Clicking \"calculate\" will open the interactive diagram in a new window."),
        actionButton("button", "Calculate"),
        tags$hr(style="border-color: grey;"),
        h4("Identify sequences"),
        p("To identify points in the sequence space, use the button below,
          then drag a box around the sequences you are interested in."),
        actionButton("button2", "Activate selection mode")
        ),
  
      
      # Outputs section 
      mainPanel(
        
        # Display summary text
        h4(textOutput("title")),
        textOutput("report"),
        tableOutput("table"),
        # plotOutput("histplot",height = "250px"),
        
        # Display fasta alignment
        tags$hr(style="border-color: grey;"),
        h4("Nearest neighbours"),
        p("These are the sequences that are closest to the query sequence within the sequence space."),
        sliderInput("return_nearest","Number of nearest neighbours",value=10,min=1,max=100),
        htmlOutput("fasta",style = "font-family:Courier New ; min-width:500px"),
        
        # Display selected points alignment
        tags$hr(style="border-color: grey;"),
        h4("Selected sequences"),
        htmlOutput("fasta2",style = "font-family:Courier New ; min-width:500px")
      )
    )
  ))),
  
  # Footer --------------------------
  wellPanel(
    div(class="footer",
        h3("References"),
        p("1) Shafee T, et al. \"The defensins consist of two independent, convergent protein superfamilies.\" ",
          tags$i("Molecular Biology and Evolution"),
          "(2016): msw106.",
          a(href="https://doi.org/10.1093/molbev/msw106","DOI:10.1093/molbev/msw106")),
        p("2) Shafee T, et al. \"Convergent evolution of defensin sequence, structure and function.\"",
          tags$i("Cellular and Molecular Life Sciences"),
          " (2017): 1-20.",
          a(href="https://doi.org/10.1007/s00018-016-2344-5","DOI:10.1007/s00018-016-2344-5"))
    )
  ),
  p("DefSpace was developed by Thomas Shafee at the ",
    a(href="http://www.latrobe.edu.au/lims","La Trobe Institute of Molecular Science (LIMS)")," and ",
    a(href="http://hexima.com.au/","Hexima")),
  img(src="http://cysbar.science.latrobe.edu.au/img/lims.png",width="200px"),
  img(src="http://cysbar.science.latrobe.edu.au/img/hexima.png",width="200px"),
  plotOutput("mainplot"),
  plotOutput("mainplot.view")
))
