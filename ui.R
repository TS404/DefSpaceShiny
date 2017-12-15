
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
    )
))
