#' msigdb_search
#'
#' You can search MsigDB gene sets on R from new shiny tab.
#'
#' @importFrom shiny jsonlite DT DBI RSQLite dplyr function
#' @param species A specie, Hu or Mm.
#' @return The gene set's name list that you selected.
#' @export

msigdb_search <- function(species = NULL) {
  return_env <- new.env()

  ui <- fluidPage(
    h3("MsigDB Browser", style = "font-size: 24px; font-weight: bold; color: #007BFF"),

    tags$script(HTML("
    $(document).on('click', 'a.dynamic-link', function(e) {
      e.preventDefault();
      var baseUrl = $(this).attr('href');
      var keyword = $('#keyword').val();
      var url = baseUrl + '?keywords=' + encodeURIComponent(keyword);
      window.open(url, '_blank');
      return false;
    });
  ")),


    sidebarLayout(
      sidebarPanel(
        width = 3,
        h3("⚙️ MENU", style = "font-weight: bold"),
        h4("Selected Gene Sets", style ="font-weight: bold; text-align: center; border-top: 1px solid #ccc;border-bottom: 1px solid #ccc;  padding: 10px 0;"),
        DTOutput("selected.df"),
        br(), br(),
        div(
          style = "display: flex; gap: 10px; justify-content: center;",
          actionButton("select", "Add", class ="btn btn-success btn-lg",
                       style = "width: 50%; font-size: 15x; text-align: center;"),
          actionButton("remove", "Remove", class="btn btn-danger btn-lg",
                       style = "width: 50%; font-size: 15px; text-align: center;")
        ),
        br(),
        actionButton("export", "Export", class="btn btn-primary", style = "width: 100%"),
        br(), br(),
        actionButton("reset", "Reset", class = "btn btn-danger btn-sm", style = "width: 100%;"),
        br(), br(),
        downloadButton("exportcsv", "Export to csv", class="btn btn-secondary btn-sm", style = "width: 100%"),
        br(), br(),
        actionButton("gotomsigdb", "Go to MsigDB.org", class = "btn btn-secondary btn-sm", style = "width: 100%;"),
        br(), br(),
        actionButton("exit", "Exit", class="btn btn-danger btn-lg", style = "width: 100%;")
      ),


      mainPanel(
        width = 9,

        br(),

        fluidRow(
          column(2, textInput("keyword", "Keyword: ")),
          column(2, selectInput("species", "Species: ",
                                choices = c("all", "HS", "MM", "RN", "RM", "DR"),
                                selected = "all")),
          column(4, selectInput("collection", "Collection: ",
                                choices = if(is.null(species)) {
                                  c("all", "C2:CGP", "C1", "C2:CP:BIOCARTA", "C2:CP:KEGG_LEGACY", "C3:MIR:MIRDB",
                                    "C3:MIR:MIR_LEGACY", "C3:TFT:GTRD", "C3:TFT:TFT_LEGACY",
                                    "C2:CP:REACTOME", "C2:CP:WIKIPATHWAYS", "C2:CP", "C4:CGN",
                                    "C4:CM", "C6", "C4:3CA", "C7:IMMUNESIGDB", "C7:VAX", "C8",
                                    "C5:GO:BP", "C5:GO:CC", "C5:GO:MF", "H", "C5:HPO", "C2:CP:KEGG_MEDICUS")
                                }  else if(species == "Mm") {
                                  c("all", "M2:CGP", "M5:GO:BP", "M3:MIRDB", "M5:MPT", "M8", "MH",
                                    "M2:CP:REACTOME", "M2:CP:WIKIPATHWAYS", "M5:GO:CC",
                                    "M5:GO:MF", "M1", "M2:CP:BIOCARTA", "M3:GTRD")
                                } else {
                                  c("all", "C2:CGP", "C1", "C2:CP:BIOCARTA", "C2:CP:KEGG_LEGACY", "C3:MIR:MIRDB",
                                    "C3:MIR:MIR_LEGACY", "C3:TFT:GTRD", "C3:TFT:TFT_LEGACY",
                                    "C2:CP:REACTOME", "C2:CP:WIKIPATHWAYS", "C2:CP", "C4:CGN",
                                    "C4:CM", "C6", "C4:3CA", "C7:IMMUNESIGDB", "C7:VAX", "C8",
                                    "C5:GO:BP", "C5:GO:CC", "C5:GO:MF", "H", "C5:HPO", "C2:CP:KEGG_MEDICUS")
                                },
                                selected = "all"))
        ),

        br(),

        fluidRow(
          column(12, DTOutput("Hu.Msigdf.filtered")))
      ))
  )




  server <- function(input, output, session) {
    dn.Msigdf <- if(species == "Hu" || is.null(species)) {
      read.csv("https://raw.githubusercontent.com/HonglabKNU/Main/refs/heads/CS/DB/msigdb_v2024.1.Hs.csv")
    } else if(species == "Mm") {
      read.csv("https://raw.githubusercontent.com/HonglabKNU/Main/refs/heads/CS/DB/msigdb_v2024.1.Mm.csv")
    } else {
      read.csv("https://raw.githubusercontent.com/HonglabKNU/Main/refs/heads/CS/DB/msigdb_v2024.1.Hs.csv")
    }

    Hu.Msigdf <- dn.Msigdf[, -1]
    Hu.Msigdf$Detail_link <- paste0(
      '<a href="https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/',
      Hu.Msigdf$standard_name,
      '.html" class="dynamic-link">Detail Link</a>'
    )

    #filtering by keyword
    Hu.Msigdf.filtered <- reactive({
      filtered.data <- Hu.Msigdf %>% select(1:3, "source_species_code", "Detail_link")
      if (input$keyword != "" && !is.null(input$keyword)) {
        keyword <- tolower(input$keyword)
        filtered.data <- Hu.Msigdf %>%
          filter(if_any(c("standard_name", "description_brief", "description_full"), ~ grepl(keyword, ., ignore.case = TRUE))) %>%
          select(1:3, "source_species_code", "Detail_link")
      }

      if (input$species != "all") {
        filtered.data <- filtered.data %>% filter(source_species_code == input$species)
      }

      if (input$collection != "all") {
        filtered.data <- filtered.data %>% filter(collection_name == input$collection)
      }

      return(filtered.data)
    })

    output$Hu.Msigdf.filtered <- renderDT({
      datatable(Hu.Msigdf.filtered(), escape = FALSE, selection = "multiple",
                rownames = F
      )
    }, server = T)


    ##Sidebar
    Hu.Msigdf.selected <- reactiveVal(data.frame(standard_name = character(0)))
    Final.output <- reactiveVal(list())

    #add genesets
    observeEvent(input$select, {
      geneset.selected <- input$Hu.Msigdf.filtered_rows_selected
      if (length(geneset.selected) != 0) {
        temp.df <- Hu.Msigdf.filtered()
        selected_names <- temp.df$standard_name[geneset.selected]
        current_df <- Hu.Msigdf.selected()
        new_names <- selected_names[!(selected_names %in% current_df$standard_name)]
        new_df <- data.frame(standard_name = new_names)
        updated_df <- rbind(current_df, new_df)
        Hu.Msigdf.selected(updated_df)
      }
    })

    #selected render
    output$selected.df <- renderDT({
      datatable(Hu.Msigdf.selected(), colnames = character(0), selection = "multiple",
                options = list(dom = "t",
                               scrollX = TRUE,
                               scrollY = '200px',
                               paging = FALSE,
                               columnDefs = list(list(
                                 targets = 1,
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "if(type === 'display' && data.length > 50) {",
                                   "return '<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>';",
                                   "}",
                                   "return data;",
                                   "}"
                                 )
                               ))))
    })

    #remove selected
    observeEvent(input$remove, {
      geneset.selected <- input$selected.df_rows_selected
      if (geneset.selected != 0) {
        current.df <- Hu.Msigdf.selected()
        removed.df <- current.df[-geneset.selected, , drop = F]
        Hu.Msigdf.selected(removed.df)
      }
    })

    #reset selected
    observeEvent(input$reset, {
      Hu.Msigdf.selected(data.frame(standard_name = character(0)))
    })

    #export
    observeEvent(input$export, {
      temp.df <- Hu.Msigdf.selected()
      temp.list <- as.vector(temp.df$standard_name)
      Final.output(temp.list)

      if (length(temp.list) != 0) {
        showModal(modalDialog(
          title = "✅ Export complete!",
          "Export Done.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      } else {
        showModal(modalDialog(
          title = "⚠️ No data.",
          "There is no data.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }
    })

    #download to csv
    output$exportcsv <- downloadHandler(
      filename = function() {
        paste0("geneset_name_", Sys.Date(), ".csv")
      },
      content = function(file) {
        temp.df <- Hu.Msigdf.selected()
        write.csv(temp.df, file, row.names = FALSE)
      }
    )

    #go to msigDB
    observeEvent(input$gotomsigdb, {
      browseURL("https://www.gsea-msigdb.org/gsea/index.jsp")
    })

    #exit
    observeEvent(input$exit, {
      showModal(modalDialog(
        title = "❗ Exit",
        "Are you sure you want to exit the app?",
        easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_exit", "Exit", class = "btn btn-danger")
        )
      ))
    })

    observeEvent(input$confirm_exit, {
      removeModal()
      stopApp(Final.output())
    })
  }

  result <- runApp(shinyApp(ui, server))
  return(result)
}
msigdb.browse <- function(species, name) {

  if(species == "Hu") {
    geneset_path <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName="
  } else if (species == "Mm") {
    geneset_path <- "https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName="
  } else {
    stop("It's a specie that doesn't support.")
  }

  temp.list <- c()

  temp.list <- lapply(name, function(n) {
    geneset_url <- paste0(geneset_path, n, "&fileType=json")
    temp <- fromJSON(geneset_url)
    temp_gene <- c(temp[[1]]$geneSymbols)
    geneset <- data.frame(gene = unlist(temp_gene))
    return(geneset)
  })




  # if you wanna delete 'gene' from 1st row,
  # geneset <- data.frame(gene = temp_df[[1]]$geneSymbols)

  names(temp.list) <- name
  return(temp.list)
}
