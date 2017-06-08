library(shiny)
library(Biobase)
library(data.table)
library(rjson)
library(DT)
library(CBMRtools)
library(shinyBS)


#get top and bottom matches for given eset of connectivity scores
summarize_eset<-function(mat, 
  summarize.func = c("mean", "median", "max", "min"),
  do.scorecutoff = TRUE, scorecutoff = c(-0.6, 0.6), 
  do.nmarkers = TRUE, nmarkers = c(100, 100)
  ){

  summarize.func<- match.arg(summarize.func)
  x<-apply(mat, 1, match.fun(summarize.func))
  x<-as.numeric(x)
  n<-length(x)
  
  if(do.nmarkers){
    ord<-order(x, decreasing = TRUE)
    x.ind.nmarkers<-c(ord[1:nmarkers[1]], ord[(n-nmarkers[2]+1):n])
  } else
    x.ind.nmarkers<-1:n

  if(do.scorecutoff)
    #TODO: rank by score here too
    x.ind.scorecutoff<-which(x > scorecutoff[2] | x < scorecutoff[1])
  else
    x.ind.scorecutoff<-1:n

  inds<-intersect(x.ind.nmarkers, x.ind.scorecutoff)
  inds<-inds[order(x[inds], decreasing = TRUE)]
  return(list(inds = inds, scores = x[inds]))
} 


get_gutc<-function(input, tab, header, datlist, sort.by){
  i<-get_BUID(input, tab)
  tab<-datlist[[header]][[i]]
  tab<-clean_gutc(tab, sort.by)
  tab<-data.table.round(tab)
  return(tab)
}

capitalize <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

clean_gutc<-function(tab, col){
  if("pert_idose" %in% colnames(tab))
    tab[, "pert_idose"]<-gsub("<fd><fd>", "u", tab[, "pert_idose"])
  

  colmatch<-paste("score_", capitalize(col), sep = "")
  scores<-tab[, colmatch]
  scores.dir<-sapply(scores, function(i){
    if(i>0) return("up")
    else return("down")
    })
  tab$summaryscore<-scores
  tab$direction<-as.character(scores.dir)
  cols_first<-c("id", "summaryscore", "direction")
  cols_others<-setdiff(colnames(tab), cols_first)
  tab[, c(cols_first, cols_others)]
}

get_chemical_description<-function(input, 
  tab,
  cols.keep = c("BUID", "Chemical.name", "CAS", "Broad_external_Id", "carc_liv_final")){
  i<-get_BUID(input, tab)
  res<-tab[tab$BUID %in% i, cols.keep]
  return(data.table.round(res[1,]))
}

subset_fdat<-function(fdat, keyword, x.split){
  if(keyword %in% "all")
    return(fdat)
  else 
    inds<-grep(keyword, x.split)
  return(fdat[inds,])
}

get_cellines<-function(fdat){
  x<-as.character(fdat$id)
  x.split<-unlist(lapply(x, function(i){strsplit(i, split = "_")[[1]][2]}))
  unique(as.character(x.split))
}

get_BUID<-function(input, tab){
  as.character(tab[which(apply(tab, 1, function(i) any(i %in% input)))[1], "BUID"])
}

get_ids_pdat<-function(pdat, cols = c("Chemical.name", "CAS", "BUID"), col.unique = "BUID",
  val.ignore = c("", " ", NA, "NA", "NOCAS")){
  tab<-unique(pdat[, cols])
  res<-lapply(cols, function(i){
  x<-as.character(tab[,i])
  x.uniq<- setdiff(unique(x), union(x[duplicated(x)], val.ignore))
  })
  names(res)<-cols
  return(res)
}

summarize_gsproj<-function(eset, order.col = "median"){
  res<-apply(exprs(eset), 1, summary)
  res<-t(res)
  colnames(res)<-c("min","Q1", "median", "mean","Q3","max")
  
  mat<-exprs(eset)
  colnames(mat)<-paste("gsscore_",pData(eset)$pert_idose, sep = "")

  res<-cbind(rowIDs = rownames(res), fData(eset), summaryscore = res[,order.col], mat)
  res<-res[order(res$summaryscore, decreasing = TRUE),, drop = FALSE]
  res<-data.table.round(res)
  return(res)
}

get_gsproj<-function(input, gslist, tab, gsname, gsmethod){
  i<-get_BUID(input, tab)
  res<-gslist[[gsname]][[gsmethod]]
  res<-res[, res$BUID %in% i]
  return(res)
}

get_gsproj_list<-function(gsnames, gsmethods, gsdir){
  res<-lapply(names(gsnames), function(i){
    gsnameval<-gsnames[[i]]
    res2<-lapply(gsmethods, function(j){
      readRDS(paste(gsdir, "/", gsnameval, "_",j, ".RDS", sep = ""))
      })
    names(res2)<-gsmethods
    return(res2)
    })
  names(res)<-names(gsnames)
  return(res)
}

get_de<-function(input, tab, 
  eset, landmark = FALSE, do.scorecutoff = TRUE, scorecutoff = c(-2, 2), 
  do.nmarkers = TRUE, nmarkers = c(100, 100),
  summarize.func = c("mean", "median", "max", "min")){
  i<-get_BUID(input, tab)
  eset<-eset[, eset$BUID %in% i]
  
  if(landmark)
    eset<-eset[fData(eset)$pr_is_lmark %in% "Y",]
  
  mat<-exprs(eset)
  fdat<-fData(eset)[, c("id", "pr_gene_symbol", "pr_is_lmark")]
  pdat<-pData(eset)

  colnames(mat)<-paste("modz_", pdat$pert_idose, sep = "")
  res<-summarize_eset(mat, summarize.func, do.scorecutoff, scorecutoff,
    do.nmarkers, nmarkers)

  res.ind<-res$inds
  res.scores<-res$scores
  tab<-cbind(fdat[res.ind,, drop = FALSE], summaryscore=res.scores, mat[res.ind,, drop = FALSE])
  tab<-data.table.round(tab)
}

get_de_by_gene_hist<-function(input, eset){
    rowid<-which(fData(eset)$pr_gene_symbol %in% input)[1]
    x<-as.numeric(exprs(eset)[rowid,])
    p.title<-paste("Distribution of mod Z-scores across samples for ", input, sep = "")
    df<-data.frame(x = x)
    p<-ggplot(df, aes(x = x))+ geom_histogram(binwidth = 0.2)+
    xlab("moderated Z-score") + 
    ylab("Count")+ 
    ggtitle(p.title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
    return(p)
}

get_de_by_gene_table<-function(input, eset){
    rowid<-which(fData(eset)$pr_gene_symbol %in% input)[1]
    x<-as.numeric(exprs(eset)[rowid,])
    pdat<-pData(eset)
    cols<-c("sig_id", "Chemical.name", "CAS", "pert_idose", "pert_dose_RANK", "distil_ss", "distil_ss_RANK")
    df<-cbind(pdat[, cols], modz = x)
    df<-df[order(df$modz, decreasing = TRUE),]
}

get_de_by_gene_lm<-function(input, eset){
   lms<-fData(eset)$pr_is_lmark
   rowid<-which(fData(eset)$pr_gene_symbol %in% input)[1]
   if(lms[rowid] %in% "Y") return(paste(input, " is a landmark gene", sep = ""))
   else paste("Warning: ", input, " is an inferred gene!", sep = "")
}

filtereset<-function(eset, filteropt, tab){
  if (filteropt %in% "all")
    eset <- eset
  else if (filteropt %in% "carcinogens") 
    eset<-eset[, eset$carc_liv_final %in% "POSITIVE"]
  else if (filteropt %in% "non-carcinogens") 
    eset<-eset[, eset$carc_liv_final %in% "NEGATIVE"]
  else if (filteropt %in% "genotoxic")
    eset<-eset[, eset$Genotoxicity %in% "POSITIVE"]
  else if (filteropt %in% "non-genotoxic")
    eset<-eset[, eset$Genotoxicity %in% "NEGATIVE"]
  else if (filteropt %in% "genotoxic carcinogens") 
    eset<-eset[, eset$carc_liv_final %in% "POSITIVE" & eset$Genotoxicity %in% "POSITIVE"]
  else if (filteropt %in% "nongenotoxic carcinogens") 
    eset<-eset[, eset$carc_liv_final %in% "POSITIVE" & eset$Genotoxicity %in% "NEGATIVE"]
  else if (filteropt %in% "genotoxic non-carcinogens") 
    eset<-eset[, eset$carc_liv_final %in% "NEGATIVE" & eset$Genotoxicity %in% "POSITIVE"]
  else if (filteropt %in% "non-genotoxic non-carcinogens") 
    eset<-eset[, eset$carc_liv_final %in% "NEGATIVE" & eset$Genotoxicity %in% "NEGATIVE"]
  else if (filteropt %in% "D. Sherr requested chemicals")
    eset<-eset[, grep("collaborator_David_Sherr", eset$source)]
  else if (filteropt %in% "J.Schlezinger requested chemicals")
    eset<-eset[, grep("collaborator_J_Schlezinger", eset$source)]
  else if (filteropt %in% "signal strength > 4")
    eset<-eset[, eset$distil_ss > 4]
  else if (filteropt %in% "signal strength > 5")
    eset<-eset[, eset$distil_ss > 5]
  else if (filteropt %in% "signal strength > 6")
    eset<-eset[, eset$distil_ss > 6]
  else if (filteropt %in% "q75 replicate correlation > 0.2")
    eset<-eset[, eset$q75rep > 0.2]
  else if (filteropt %in% "q75 replicate correation > 0.3")
    eset<-eset[, eset$q75rep > 0.3]
  else if (filteropt %in% "tas > 0.2")
    eset<-eset[, eset$tas > 0.2]
  else if (filteropt %in% "tas > 0.4")
    eset<-eset[, eset$tas > 0.4]
  else if (filteropt %in% "tas > 0.6")
    eset<-eset[, eset$tas > 0.6]
  else 
    eset<-eset[, eset$BUID %in% get_BUID(filteropt, tab)]
  return(eset)
}

subset_names<-function(x,n){
  if(nchar(x)> n) return(paste(substr(x, 1, n), "...", sep = ""))
  else return(x)
}

get_morpheus_link<-function(url, domain){
  url = paste(domain, "/", url, sep = "")
  url = paste("{\"dataset\":", "\"", url, "\"}", sep ="")
  url = URLencode(URL = url)
  url = paste("https://software.broadinstitute.org/morpheus/?json=", url, sep = "")
  return(url)
}

get_heatmap_gct<-function(ds, dsmap, method, domain){
  dsname<-dsmap[[ds]]
  url<-paste(dsname, "_", method, ".gct", sep = "")

}

get_heatmap_eset<-function(ds, dsmap, method){
  dsname<-dsmap[[ds]]
  res<-paste(dsname, "_", method, ".RDS", sep = "")
  return(res)
}

data.table.round<-function(dt, digits = 4){
  cols<-sapply(colnames(dt), function(i) is.numeric(dt[,i]))
  cols<-names(which(cols))

  for(i in cols){
    dt[,i]<-round(dt[,i], digits)
  }
  dt<-data.table(dt)
}

##load data dirs
dirs<-fromJSON(file = "datadirs.json")

##load data for chemical annotation tab
chemannot<-readRDS(file = dirs$chemannotation_filename)

##load data for Differential Expression tab
deeset<-readRDS(dirs$diffexp_filename)

#chemical selection dropdown
chemicals<-get_ids_pdat(chemannot)

#genes selection dropdown
genes<-unique(fData(deeset)$pr_gene_symbol)

#summarization function
summarizefuncs<-c("max", "median", "mean", "min")

#directory of heatmap data files
heatmapdir<-dirs$genesetenrich_dir
heatmapfiles<-gsub(".RDS", "",list.files(heatmapdir))


premade_sets<-c("carcinogens", "non-carcinogens", "genotoxic", "non-genotoxic",
  "genotoxic carcinogens", "nongenotoxic carcinogens",
  "genotoxic non-carcinogens", "non-genotoxic non-carcinogens",
  "D. Sherr requested chemicals", "J.Schlezinger requested chemicals",
  "signal strength > 4", "signal strength > 5", "signal strength > 6",
  "q75 replicate correlation > 0.2",
  "q75 replicate correation > 0.3",
  "tas > 0.2", "tas > 0.4", "tas > 0.6")

filteropts<-c(list(premade_sets = premade_sets), chemicals)

domain<-dirs$gct_dir

dsmap<-list(Hallmark="gsscores_h.all.v5.0",
    C2="gsscores_c2.cp.reactome.v5.0", 
    NURSA="gsscores_nursa_consensome_Cbyfdrvalue_0.01.gmt")

gctfiles<-names(dsmap)
gctmethods<-c("gsproj", "gsva", "ssgsea", "zscore")

##load data for GeneSetEnrichment tab
gsnames<-dsmap
gsmethods<-gctmethods
gsdir<-dirs$genesetenrich_dir
gslist<-get_gsproj_list(gsnames, gsmethods, gsdir)
gssort<-c("min","Q1", "median", "mean","Q3","max")

##gutc data
gutcdir<-dirs$gutc_dir
gutcfiles<-list.files(gutcdir)
gutcheaders<-gsub(".RDS", "", gutcfiles)
gutcobjects<-lapply(gutcfiles, function(i){
  readRDS(paste(gutcdir, "/", i, sep = ""))
  })
names(gutcobjects)<-gutcheaders

#tooltip texts
helptextgutc<-HTML(paste("cs: raw weighted connectivity scores",
              "ns: normalized scores, accounts for cell-line and perturbational type",
              "ps: percentile normalized scores[-100, 100]",
              "pcl: PCL (perturbational classes)", 
              "pert: perturbagen level",
              "cell: cell-line level",
              "summary: cell line-summarized level", sep="<br/>"))

helptextgsname<-HTML(paste("Hallmark: MSigDB Hallmark Pathways (v5.0)",
              "C2: MSigDB C2 reactome Pathways (v5.0)",
              "NURSA: Nuclear Receptor Signaling Atlas, consensome data for human", 
               sep="<br/>"))

helptextgsmethod<-HTML(paste("gsproj: GeneSetProjection for R package montilab:CBMRtools",
              "gsva, ssgea, zscore: from R Bioconductor package GSVA", 
               sep="<br/>"))


helptextgsfilter<-HTML(paste("filter columns by premade sets or by chemical",
              "premade sets:",
              "carcinogenicity/genotoxicity based on CPDB mouse and rats data",
              "D. Sherr suggested chemicals: mainly AHR ligands",
              "J. Schlezinger suggested chemicals: mainly PPAR ligands",
               sep="<br/>"))

selectInputWithTooltip<-function(inputId, label, choices, bId, helptext, ...){
  selectInput(inputId, tags$span(label,  
            tipify(bsButton(bId, "?", style = "inverse", size = "extra-small"), 
              helptext)),
             choices, ...)}

##define app
app<-shinyApp(

ui = shinyUI(
  fluidPage(
  tags$head(includeScript("google-analytics.js")),
  navbarPage("CRCGN Portal",
    tabPanel("About",
      titlePanel("CRCGN Liver Portal"),
      fluidRow(
        column(11,
          includeMarkdown("introduction.Rmd")
        ),
        img(src="logo.png", align = "left", width = 600)
      )
    ),

    tabPanel("Chemical Annotation",
      DT::dataTableOutput("chemannot_result")
      ),

    tabPanel("Differential Expression",
      fluidPage(
        #search by chemical
        radioButtons("detype", "Search Type", 
          choices = c("by chemical", "by gene"), 
          selected = "by chemical"),

        conditionalPanel(
          condition = "input.detype ==  'by chemical'",

          #dropdown menus
          fluidRow(
            column(3, selectInput("chemical_de", "CRCGN Chemical:", chemicals)),
            column(2, checkboxInput("landmark_de", "Landmark only", value =FALSE)),
            column(2, selectInput("summarizefunc_de", "Summarization:",summarizefuncs,
              selected = "median")),
            
            column(1,checkboxGroupInput("filterbyinput_de", "Filter by:",
                         c("score" = "score",
                           "number" = "number"),
                         selected = c("score", "number"))),
            column(2,sliderInput("range_de", "score threshold", min = -10, max = 10, 
              value = c(-2,2), step = 0.01)),
            column(1,sliderInput("numberthresleft_de", "Num +",
              min = 0, max = 1000, value = 10, ticks = FALSE, step = 10)),
            column(1,sliderInput("numberthresright_de", "Num -",
              min = 0, max = 1000, value = 10, ticks = FALSE, step = 10))
          ),
          # Table returned bshowing description for query chemical
          dataTableOutput("chemical_description_de"),
          tags$style(type="text/css", '#chemical_description_de tfoot {display:none;}'),

          ##insert empty space
          fluidRow(column(width = 1, offset = 0, style='padding:10px;')),

          # Table returned showing de markers
          dataTableOutput("result_de")
        ),

        #search by gene
        conditionalPanel(
          condition = "input.detype == 'by gene'",
          fluidRow(
            column(3, selectInput("gene_de", "Gene symbol:", genes)),
            column(3, textOutput("gene_de_lm"))
          ),
          tags$style(type='text/css', "#gene_de_lm { width:100%; margin-top: 25px;}"),
          plotOutput("gene_de_hist"),
          dataTableOutput("gene_de_table")
        )
      )
    ),
		
		tabPanel("Gene set enrichment",
      fluidPage(
        fluidRow(
          column(3, selectInput("chemical_gs", "CRCGN Chemical:", chemicals)),
          column(2, 
              selectInputWithTooltip(inputId = "gsname_gs", label = "Gene set name", 
              choices = names(gsnames), bId = "Bgsname", helptext =helptextgsname)
            ),
          column(2,           
            selectInputWithTooltip(inputId = "gsmethod_gs", label = "Projection method", 
              choices = gsmethods, bId = "Bgsmethod", helptext =helptextgsmethod)
          ),
          column(2, selectInput("summarize_gs", "Sort by:", gssort,selected = "median"))
          )
        ),
      DT::dataTableOutput("gsproj_result")
    ),

    tabPanel("Heatmap (interactive)",
      fluidRow(
        column(3, 
              selectInputWithTooltip(inputId = "gctfile", label = "Dataset", 
              choices = gctfiles, bId = "Bgsnamegct", helptext =helptextgsname)
          ),
        column(2,
              selectInputWithTooltip(inputId = "gctmethod", label = "Projection method", 
              choices = gctmethods, bId = "Bgsmethodgct", helptext =helptextgsmethod)
          )
        ),
      htmlOutput("morpheus_result_link"),
      htmlOutput("morpheus_result_embedded")
      ),

    tabPanel("Heatmap (static)",
      fluidPage(
        fluidRow(         
        column(3, 
              selectInputWithTooltip(inputId = "gsfile_heatmap", label = "Dataset", 
              choices = gctfiles, bId = "Bgsnamegctstatic", helptext =helptextgsname)
              ),
        column(2, 
              selectInputWithTooltip(inputId = "gsmethod_heatmap", label = "Projection method", 
              choices = gctmethods, bId = "Bgsmethodgctstatic", helptext =helptextgsmethod)
          ), 
        column(3, 
          selectInputWithTooltip(inputId = "filteropt", label = "Column Filter", 
          choices = filteropts, 
          selected = "signal strength > 6", bId = "Bgsfilterstatic", helptext = helptextgsfilter)
          )
        ),
        plotOutput("heatmap_result")

      )
     ),

     tabPanel("Connectivity",
       fluidPage(
         fluidRow(         
          column(3, selectInput("chemical_gutc", "CRCGN Chemical:", chemicals)),
          column(3,
            selectInputWithTooltip(inputId = "header_gutc", label = "Dataset", 
              choices = gutcheaders, bId = "Bgutc", helptext =helptextgutc)
            ),
          column(2, selectInput("summarize_gutc", "Sort by:", gssort, selected = "median"))
         ),
         dataTableOutput("chemical_description_gutc"),
         tags$style(type="text/css", '#chemical_description_de tfoot {display:none;}'),

         ##insert empty space
         fluidRow(column(width = 1, offset = 0, style='padding:10px;')),

         dataTableOutput("gutc_result")

       )
     )

	))),

server = shinyServer(function(input, output, session) {

  output$chemannot_result<-DT::renderDataTable({data.table.round(chemannot)},
    extensions = 'Buttons', 
    server = TRUE,
    options = list(dom = 'T<"clear">Blfrtip', 
    deferRender=FALSE, 
    scrollX = TRUE,
    scrollY = 400,
    scrollCollapse = TRUE,
    pageLength = 1000,
    lengthMenu = c(10,50,100,200,1000),
    buttons=c('copy','csv','print')))

  output$gene_de_table<-renderDataTable({
    data.table.round(get_de_by_gene_table(input$gene_de, deeset))
    },
    extensions = 'Buttons', 
    server = TRUE,
    options = list(dom = 'T<"clear">Blfrtip', 
    deferRender=FALSE, 
    scrollX = TRUE,
    scrollY = 400,
    scrollCollapse = TRUE,
    pageLength = 1000,
    lengthMenu = c(10,50,100,200,1000),
    buttons=c('copy','csv','print'))
   )

  output$gene_de_hist<-renderPlot({
    get_de_by_gene_hist(input$gene_de, deeset)
    }, width = 1000, height = 300)

  output$gene_de_lm<-renderText({
    get_de_by_gene_lm(input$gene_de, deeset)
    })

  output$gsproj_result <- DT::renderDataTable({
      eset<-get_gsproj(input$chemical_gs, gslist, chemannot, input$gsname_gs, input$gsmethod_gs)
      res<-summarize_gsproj(eset, order.col = input$summarize_gs)
      return(res)
    }, 
    extensions = 'Buttons',
    server = FALSE,
    options = list (dom = 'T<"clear">Blfrtip',
      autoWidth = FALSE,
      columnDefs = list(list(width = '100px', targets = list(2), className = 'dt-center',
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 30 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 30) + '..</span>' : data;",
            "}")), 
      list(targets = list(1), visible = FALSE)
      ),
      deferRender=TRUE,
      scrollCollapse=TRUE,
      pageLength = 10, lengthMenu = c(10,50,100,200,1000),
      buttons=c('copy','csv','print')
    )
  )

  output$chemical_description_de<-DT::renderDataTable({
    get_chemical_description(input = input$chemical_de, 
      tab = chemannot,
      cols.keep = c("BUID", "Chemical.name", "CAS", "Broad_external_Id", "carc_liv_final"))

    }, options = list(dom = ''))

  output$chemical_description_gutc<-DT::renderDataTable({
    get_chemical_description(input = input$chemical_gutc, 
      tab = chemannot,
      cols.keep = c("BUID", "Chemical.name", "CAS", "Broad_external_Id", "carc_liv_final"))

    }, options = list(dom = ''))


  output$result_de<-DT::renderDataTable({
      get_de(input$chemical_de, tab=chemannot, 
        eset = deeset, 
        landmark = input$landmark_de, 
        do.scorecutoff = "score" %in% input$filterbyinput_de, 
        scorecutoff = c(input$range_de[1], input$range_de[2]), 
        do.nmarkers = "number" %in% input$filterbyinput_de, 
        nmarkers = c(input$numberthresleft_de, input$numberthresright_de),
        summarize.func = input$summarizefunc_de)
    },
    extensions = 'Buttons', 
    server = TRUE,
    options = list(dom = 'T<"clear">Blfrtip', 
    deferRender=FALSE, 
    scrollX = TRUE,
    scrollY = 400,
    scrollCollapse = TRUE,
    pageLength = 10,
    lengthMenu = c(10,50,100,200,1000),
    buttons=c('copy','csv','print')))

  output$heatmap_result<-renderPlot({
      outfileheader<-get_heatmap_eset(ds = input$gsfile_heatmap, dsmap = dsmap, method = input$gsmethod_heatmap)
      outfile<-paste(heatmapdir, "/", outfileheader, sep = "")
      eset<-readRDS(outfile)
      eset<-filtereset(eset, input$filteropt, chemannot)

      cns<-as.character(sapply(as.character(eset$Chemical.name), 
        function(i) subset_names(i,25)))
      cns<-make.unique(cns) #in case truncation produces non-unique ids
      colnames(eset)<-paste(cns, "_", eset$pert_dose_RANK, sep = "")
      hc<-clust_eset(eset)
      col_legend<-list(carc_liv_final = list(col_breaks = c("POSITIVE", "NEGATIVE", ""), 
        col_values = sapply(c("orange", "green", "white"), to.hex),
        col_labels = c("POSITIVE", "NEGATIVE", "")),
      Genotoxicity = list(col_breaks = c("POSITIVE", "NEGATIVE", ""), 
        col_values = sapply(c("orange", "green", "white"), to.hex),
        col_labels = c("POSITIVE", "NEGATIVE", "")),
      pert_dose_RANK = list(col_breaks = 1:6, 
        col_values = rev(gray.colors(6)),
        col_labels = 1:6))
      hmcolors<-function(...) scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, ...)
    
      ncols<-ncol(eset)
      xsize<-3
      if (ncols<500) 
        xsize<-4
      if (ncols<100)
        xsize<-6
      if(ncols<10)
        xsize<-10
      if(ncols > 1000)
        xsize<-0 

      p<-ggheat.continuous.single(eset = eset, 
        hc = hc$hc, 
        hr = hc$hr, 
        hmcolors = hmcolors,
        hmtitle = "geneset score",
        col_lab = c("carc_liv_final", "Genotoxicity", "pert_dose_RANK"), 
        col_legend = col_legend,
        ylabstr = "",
        fout = NA, 
        p.heights = c(1.5, 0.5, 5),
        xsize = xsize,
        ysize = 3, 
        ysizelab = 7,
        xright = 0.18)
      return(p)
  }, width = 1000, height = 600)

  output$morpheus_result_link<-renderText({
      paste(c('<a target="_blank" href="',
      get_morpheus_link(url =get_heatmap_gct(ds=input$gctfile, dsmap = dsmap, 
          method = input$gctmethod), 
        domain = domain),
        '">', 'click here to open in new tab', '</a>'), sep = "")
      })

  output$morpheus_result_embedded<-renderText({
      paste(c('<iframe width ="1200" height ="680" src="',
      get_morpheus_link(url = get_heatmap_gct(ds=input$gctfile, dsmap = dsmap, 
          method = input$gctmethod), 
        domain = domain),
        '">', '</iframe>'), sep = "")
      })

  output$gutc_result<-DT::renderDataTable({
      get_gutc(input = input$chemical_gutc, 
        tab = chemannot, 
        header = input$header_gutc, 
        datlist = gutcobjects,
        sort.by = input$summarize_gutc)
    },    
    extensions = 'Buttons', 
    server = TRUE,
    options = list(dom = 'T<"clear">Blfrtip', 
    deferRender=FALSE, 
    scrollX = TRUE,
    scrollY = 400,
    scrollCollapse = TRUE,
    pageLength = 50,
    lengthMenu = c(10,50,100,200,1000),
    buttons=c('copy','csv','print')))

})

)

