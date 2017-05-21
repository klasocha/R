library(shiny)
library(ggplot2)

setwd("/home/kacper/R/gene_atlas")
expr <- read.csv(file="expr.csv", header=TRUE, sep=",")
anno <- read.csv(file="anno.csv", header=TRUE,  sep="\t")

gene_symbol <- anno$Symbol
gene_symbol <- unique(gene_symbol)

anno = anno[match(gene_symbol, anno$Symbol),]

expr = expr[match(anno$ProbesetID, expr$X),]
expr$X <- NULL
gene_symbol <- anno$Symbol
gene_desc <- anno$Description
remove(anno)


mean_expr = rowMeans(expr, na.rm = FALSE, dims = 1)




ui <- fluidPage(
      
   
    
    fluidRow(
      column(7,
             plotOutput("plot1", click = "plot1_click")
      ),
      column(5,
             selectInput("variable", "Gene:",
                         gene_symbol),
             htmlOutput("x_value"),
             dataTableOutput('table1')
      )),
    fluidRow(
      
      htmlOutput("gene_info")
      
    )
  )
  
  
server <- function(input, output) {
  output$plot1 <- renderPlot({
    
    index = match(input$variable, gene_symbol)
    bar_data = as.numeric(expr[index,])
    
    values = -as.numeric(unlist(sort(-bar_data)))[1:5]
    labels = colnames(expr)[sort(-bar_data, index.return=TRUE)[[2]][1:5]]
    
    
    mydata <- data.frame(labels, values)
    
    p <-ggplot(mydata, aes(reorder(labels,-values), values))
    p +geom_bar(stat = "identity") +
      xlab("Tissue") + ylab("Expression [a.u.]") +
      ggtitle(paste(input$variable, "is located mostly in:"))
  })
  

output$x_value <- renderText({
    if (is.null(input$plot1_click$x)) return("")
    else {
      lvls <- levels(ToothGrowth$supp)
      
      index = match(input$variable, gene_symbol)
      bar_data = as.numeric(expr[index,])
      labels = colnames(expr)[sort(-bar_data, index.return=TRUE)[[2]][1:5]]
      
      name <- labels[round(input$plot1_click$x)]
      HTML("Selected  tissue: <code>", name, "</code>")
    }
  })
  
  
  output$table1 <- renderDataTable({
    if (is.null(input$plot1_click$x)) return()
    else {
      
      index = match(input$variable, gene_symbol)
      bar_data = as.numeric(expr[index,])
      
      labels = colnames(expr)[sort(-bar_data, index.return=TRUE)[[2]][1:5]]
      
      name <- labels[round(input$plot1_click$x)]
      
      tissue_index = match(name, colnames(expr)) 
      
      
      ratio = expr[,tissue_index]/mean_expr
      
      top_genes = gene_symbol[sort(-ratio, index.return=TRUE)[[2]]]
      top_genes = as.data.frame(head(top_genes, 6))
      names(top_genes) <- c("Other genes highly present:")
    return(top_genes)
    }
  }, escape = FALSE, options = list(paging = FALSE, searching = FALSE, info = FALSE, colnames = FALSE))
  
  
output$gene_info <- renderText({
    
  
  index = match(input$variable, gene_symbol)
  
  gene_corr = cor(t(expr)[,index], t(expr)[,-index])
  
  corr_genes_names = (gene_symbol[-index])[sort(-gene_corr, index.return=TRUE)[[2]][1:5]]
      
  HTML("Gene description: <code>", as.character(gene_desc[index]), "</code> <br>",
       "Correlated genes: <code>", paste(corr_genes_names, collapse=", "), "</code>")
  })
  
}

shinyApp(ui, server)