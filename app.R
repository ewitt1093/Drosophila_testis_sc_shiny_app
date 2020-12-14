library(shiny)
library(ggplot2)
library(plyr)
library(data.table)


#this is the data for the marker enrichment table
markers<-read.table("markers.txt")
names(markers)<-c("p_val", "avg_logFC(cell type vs all)", "Percent expressed (in cell type)", "Percent expressed (all other cells)", "p.adj", "Cell type", "gene")

#This helps us convert between gene symbol and gene id later:
genestsv<-read.table("genes.tsv")
names(genestsv)<-c("ID", "Gene")


#######Make separate files for each gene- already run on my end, this is just for documentation on how I handled the seurat object

#markers<-FindAllMarkers(testis.integrated)
#write.table(markers, "markers.txt")
#tmp<-FetchData(testis.integrated, vars= rownames(GetAssayData(testis.integrated)))
#write.csv(tmp, "expressionvalues.csv")

#data for t-SNE embeddings
#embeds<-Embeddings(testis.integrated[["tsne"]])
#idents<-as.data.frame(Idents(testis.integrated))
#embeds<-cbind(embeds, idents)
#names(embeds)<-c("tSNE_1", "tSNE_2", "Celltype")
#write.table(embeds, "embeds.txt")


##########end of setup
embeds<-read.table("embeds.txt")


ui <- fluidPage(
  tags$head(tags$style(".shiny-output-error{color: white;}")),
  titlePanel("Gene expression in Drosophila testes"),
  sidebarLayout(
    sidebarPanel(
      textInput("geneInput", "Enter a gene symbol or FBgn id",value = ""),
      sliderInput("Featureplotslider",label = "Cutoffs for T-SNE expression (arbitrary units)", min = 0, max=10, value=c(0,3)  ),
      selectInput(inputId = "minColor", label="Min color", choices = c("red", "darkred","lightblue", "darkblue", "black", "yellow", "green", "lightgreen","cyan","purple", "orange", "black", "grey"), selected = "lightblue"),
      selectInput(inputId = "maxColor", label="Max color", choices = c("red", "darkred","lightblue", "darkblue", "black", "yellow", "green", "lightgreen","cyan","purple", "orange", "black", "grey"), selected="red"),
      downloadButton(outputId = 'downloadReport')
    ),
    mainPanel(
      plotOutput("coolplot2", height=435,width=500),
      br(),
      plotOutput("coolplot3", height=400,width=600),
      br(),
      plotOutput("coolplot4", height=400,width=600),
      h3("Cell types with enrichment for this gene"), ###this doesn't work
      tableOutput('enrichment')
    )
  )
)



server <- function(input, output) {

  state <- reactiveValues()
  #If name is FBgn, convert to gene symbol
  observe({
    state$x<-input$geneInput
    if( state$x %in% genestsv$ID){
      state$y<- levels(droplevels(subset(genestsv, ID==input$geneInput )$Gene))
    }else{
      state$y<-input$geneInput
      #state$y is the gene symbol for plotting
    }

  })
  
  reac<-reactiveValues()
  observe({ if (state$y %in% genestsv$Gene){ exp<-fread(file="expressionvalues.csv", select=state$y, sep=",")
     merged<-cbind(embeds, exp)
      names(merged)<-c("tSNE_1", "tSNE_2", "Celltype", "Expression")
      merged$Celltype<-factor(merged$Celltype, levels=c("Hub cells", "Cyst cells", "Epithelial cells", "GSC, Early spermatogonia", "Late spermatogonia", "Early spermatocytes", "Late spermatocytes", "Early spermatids", "Late spermatids"))
      reac$merged<-merged
  }
  })

  

  
  #fix this part

  #enrichment table
  output$enrichment<-renderTable(subset(markers,gene==state$y ))
  
  #gene dotplot
  
  #output$coolplot <-renderPlot({
  #  DotPlot(testis.integrated, assay="RNA", features=state$y,col.min=input$Dotplotslider[1],cols = c(input$minColor, input$maxColor), col.max=input$Dotplotslider[2])})

 
  
  #t-sne plot of all cell types, and expression of target gene within all cells 
observe({ if (state$y %in% genestsv$Gene){
  output$coolplot2 <- renderPlot({ggplot(reac$merged, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Expression))+theme_classic()+scale_color_gradientn(limits=c(input$Featureplotslider[1], input$Featureplotslider[2]),oob = scales::squish,colors=c(input$minColor, input$maxColor))+ggtitle(state$y)+theme(text=element_text(size=15),plot.title=element_text(size=40, face="italic"))})
  output$coolplot3<-renderPlot({ggplot(reac$merged,aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Celltype))+theme_classic()+theme(text=element_text(size=15))+scale_color_manual(values=c("#e714bf","#9b1c9f", "#14339d", "#305ed2", "#0bf1ae", "#1a903f", "#e0ff3d", "#f3992d", "#e44f00"))})
  output$coolplot4<-renderPlot({ggplot(reac$merged, aes(x=Celltype,y=Expression))+geom_boxplot(aes(fill=Celltype))+theme_classic()+theme(text=element_text(size=15), axis.text=element_text(angle=90))+
      scale_fill_manual(values=c("#e714bf","#9b1c9f", "#14339d", "#305ed2", "#0bf1ae", "#1a903f", "#e0ff3d", "#f3992d", "#e44f00"))})
}
  else{
    output$coolplot2<-renderPlot(ggplot(reac$merged)+theme_void()+ggtitle("Please enter a valid gene (i.e p-cup)"))
    output$coolplot3<-renderPlot(ggplot(reac$merged)+theme_void())
    output$coolplot4<-renderPlot(ggplot(reac$merged)+theme_void())
  }}
)
  #options to dowload plots as PDF
  output$downloadReport = downloadHandler(
    filename = 'graphs.pdf',
    content = function(file) {
      pdf(file)
     print(ggplot(reac$merged, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Expression))+theme_classic()+
             scale_color_gradientn(limits=c(input$Featureplotslider[1], input$Featureplotslider[2]),oob = scales::squish,colors=c(input$minColor, input$maxColor))+
             ggtitle(state$y)+theme(text=element_text(size=15),plot.title=element_text(size=40, face="italic")))
     print(ggplot(reac$merged,aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=Celltype))+theme_classic()+theme(text=element_text(size=15))
           +scale_color_manual(values=c("#e714bf","#9b1c9f", "#14339d", "#305ed2", "#0bf1ae", "#1a903f", "#e0ff3d", "#f3992d", "#e44f00")))
     print(ggplot(reac$merged, aes(x=Celltype,y=Expression))+geom_boxplot(aes(fill=Celltype))+
             theme_classic()+theme(text=element_text(size=15), axis.text=element_text(angle=90)))+
       scale_color_manual(values=c("#e714bf","#9b1c9f", "#14339d", "#305ed2", "#0bf1ae", "#1a903f", "#e0ff3d", "#f3992d", "#e44f00"))
       dev.off()
    })
  
}



shinyApp(ui = ui, server = server)
