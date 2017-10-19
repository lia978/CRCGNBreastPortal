library(shiny)
library(Biobase)
library(data.table)
library(rjson)
library(DT)
library(shinyBS)

source("ggheat.continuous.R")

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
  tab<-data.table.multibyte.string(tab)
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

get_ids_pdat<-function(pdat, cols = c("Chemical Name", "CAS", "BUID"), col.unique = "BUID",
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
  summarize.func = c("mean", "median", "max", "min"),
  fdat.id = c("id", "pr_gene_symbol", "pr_is_lmark"), 
  pdat.id = c("pert_idose"), landmark.values = c("Y", "N")){
  i<-get_BUID(input, tab)
  eset<-eset[, eset$BUID %in% i]
  
  if(landmark)
    eset<-eset[fData(eset)$pr_is_lmark %in% landmark.values[1],]
  
  mat<-exprs(eset)
  fdat<-fData(eset)[, fdat.id]
  pdat<-pData(eset)

  pdat.names<- apply(pdat[, pdat.id, drop = FALSE], 1, function(i) paste(i, collapse = "_"))
  colnames(mat)<-paste("modz_", pdat.names, sep = "")
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

get_de_by_gene_table<-function(input, eset, 
  cols){
    rowid<-which(fData(eset)$pr_gene_symbol %in% input)[1]
    x<-as.numeric(exprs(eset)[rowid,])
    pdat<-pData(eset)
    df<-cbind(pdat[, cols], modz = x)
    df<-df[order(df$modz, decreasing = TRUE),]
}

get_de_by_gene_lm<-function(input, eset, values = c("Y", "N")){
   lms<-fData(eset)$pr_is_lmark
   rowid<-which(fData(eset)$pr_gene_symbol %in% input)[1]
   if(lms[rowid] %in% values[1]) return(paste(input, " is a landmark gene", sep = ""))
   else paste("Warning: ", input, " is an inferred gene!", sep = "")
}

filtereset<-function(eset, filteropt, tab, 
  carc_label = "breast_carcinogen_any", 
  geno_label = "Genotoxicity"
  ){

  pdat<-pData(eset)
  carc<-pdat[, carc_label]
  geno<-pdat[, geno_label]
  source<-pdat[, "source"]
  distil_ss<-pdat[, "distil_ss"]
  tas<-pdat[, "tas"]

  if (filteropt %in% "all")
    eset <- eset
  else if (filteropt %in% "carcinogens") 
    eset<-eset[, carc %in% "POSITIVE"]
  else if (filteropt %in% "non-carcinogens") 
    eset<-eset[, carc %in% "NEGATIVE"]
  else if (filteropt %in% "genotoxic")
    eset<-eset[, geno %in% "POSITIVE"]
  else if (filteropt %in% "non-genotoxic")
    eset<-eset[, geno %in% "NEGATIVE"]
  else if (filteropt %in% "genotoxic carcinogens") 
    eset<-eset[, carc %in% "POSITIVE" & geno %in% "POSITIVE"]
  else if (filteropt %in% "nongenotoxic carcinogens") 
    eset<-eset[, carc %in% "POSITIVE" & geno %in% "NEGATIVE"]
  else if (filteropt %in% "genotoxic non-carcinogens") 
    eset<-eset[, carc %in% "NEGATIVE" & geno %in% "POSITIVE"]
  else if (filteropt %in% "non-genotoxic non-carcinogens") 
    eset<-eset[, carc %in% "NEGATIVE" & geno %in% "NEGATIVE"]
  else if (filteropt %in% "D. Sherr requested chemicals")
    eset<-eset[, grep("collaborator_David_Sherr", source)]
  else if (filteropt %in% "J.Schlezinger requested chemicals")
    eset<-eset[, grep("collaborator_J_Schlezinger", source)]
  else if (filteropt %in% "signal strength > 4")
    eset<-eset[, distil_ss > 4]
  else if (filteropt %in% "signal strength > 5")
    eset<-eset[, distil_ss > 5]
  else if (filteropt %in% "signal strength > 6")
    eset<-eset[, distil_ss > 6]
  #else if (filteropt %in% "q75 replicate correlation > 0.2")
  #  eset<-eset[, eset$q75rep > 0.2]
  #else if (filteropt %in% "q75 replicate correation > 0.3")
  #  eset<-eset[, eset$q75rep > 0.3]
  else if (filteropt %in% "tas > 0.2")
    eset<-eset[, tas > 0.2]
  else if (filteropt %in% "tas > 0.4")
    eset<-eset[, tas > 0.4]
  else if (filteropt %in% "tas > 0.6")
    eset<-eset[, tas > 0.6]
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

data.table.multibyte.string<-function(df){
  strcol<-which(sapply(df, function(i) !is.numeric(i)))
  for(i in strcol){
    df[, i]<-as.character(iconv(df[, i]))
  }
  return(df)

}

get_genecard_link<-function(genesymbol){
  sprintf('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s&keywords=%s" target="_blank" class="btn btn-primary">%s</a>',
    genesymbol, genesymbol, genesymbol)
}

get_geneset_link<-function(geneset){
  sprintf('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/%s" target="_blank" class="btn btn-primary">%s</a>',
    geneset, geneset)

}

#make a named list
List <- function(...) {
  names <- as.list(substitute(list(...)))[-1L]
  setNames(list(...), names)
}

selectInputWithTooltip<-function(inputId, label, choices, bId, helptext, ...){
  selectInput(inputId, tags$span(label,  
            tipify(bsButton(bId, "?", style = "inverse", size = "extra-small"), 
              helptext)),
             choices, ...)}

load_data<-function(config = "datadirs.json"){

  ##load data dirs
  dirs<-fromJSON(file = config)

  ##load data for chemical annotation tab
  chemannot<-readRDS(file = dirs$chemannotation_filename)

  #remove multibyte string
  chemannot$Chemical.name<-iconv(chemannot$Chemical.name)

  colnames(chemannot)<-c("Broad_external_Id", "Mean Signal Strength", 
    "BUID", "Chemical Name", "CAS", "Genotoxicity",
   "Breast Carcinogenicity", "Lung Carcinogenicity", 
   "Breast Carcinogenicity CPDB", "CPDB Lung Carcinogenicity CPDB", 
   "source")
  chemannot<-chemannot[,   c( "BUID", "Chemical Name", "CAS", "Genotoxicity",
   "Breast Carcinogenicity", "Lung Carcinogenicity", 
   "Breast Carcinogenicity CPDB", "CPDB Lung Carcinogenicity CPDB", 
   "Mean Signal Strength", "Broad_external_Id",
   "source")]

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
    "signal strength > 4", "signal strength > 5", "signal strength > 6")

  filteropts<-c(list(premade_sets = premade_sets), chemicals)
  domain<-dirs$gct_dir

  dsmap<-list(Hallmark="gsscores_h.all.v5.0",
      C2="gsscores_c2.cp.reactome.v5.0", 
      NURSA="gsscores_nursa_consensome_Cbyfdrvalue_0.01.gmt")

  gctfiles<-names(dsmap)
  gctmethods<-c("gsva", "ssgsea", "zscore", "gsproj")

  ##load data for GeneSetEnrichment tab
  gsnames<-list(Hallmark="gsscores_h.all.v5.0",
      C2="gsscores_c2.cp.reactome.v5.0", 
      NURSA="gsscores_nursa_consensome_Cbyfdrvalue_0.01.gmt")

  gsmethods<-gctmethods
  gsdir<-dirs$genesetenrich_dir
  gslist<-get_gsproj_list(gsnames, gsmethods, gsdir)
  gssort<-c("min","Q1", "median", "mean","Q3","max")

  #gutc data
  gutcdir<-dirs$gutc_dir
  gutcfiles<-list.files(gutcdir)
  gutcheaders<-gsub(".RDS", "", gutcfiles)
  gutcobjects<-lapply(gutcfiles, function(i){
    readRDS(paste(gutcdir, "/", i, sep = ""))
    })
  names(gutcobjects)<-gutcheaders

  gutch2<-c("ps_pcl_summary", "ps_pcl_cell", "ps_pert_summary", "ps_pert_cell")
  gutch2<-c(gutch2, setdiff(gutcheaders, gutch2))

  gutcobjects<-gutcobjects[gutch2]
  gutcheaders<-names(gutcobjects)

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

  helptextgsmethod<-HTML(paste(
                "gsva, ssgea, zscore: from R Bioconductor package GSVA", 
                "gsproj: GeneSetProjection for R package montilab:CBMRtools",
                 sep="<br/>"))


  helptextgsfilter<-HTML(paste("filter columns by premade sets or by chemical",
                "premade sets:",
                "carcinogenicity/genotoxicity based on CPDB mouse and rats data",
                "D. Sherr suggested chemicals: mainly AHR ligands",
                "J. Schlezinger suggested chemicals: mainly PPAR ligands",
                 sep="<br/>"))

  fdat.id<-c("pr_id", "pr_gene_symbol", "pr_is_lmark", "reactome.category")
  pdat.id<-c("pert_idose", "genotype")

  dat <- List(chemannot, deeset, chemicals, genes, summarizefuncs, 
      heatmapdir, heatmapfiles, premade_sets, filteropts, domain, 
      dsmap, gctfiles, gctmethods, 
      gsnames, gsmethods, gsdir, gslist, gssort,
      helptextgutc, helptextgsname, helptextgsmethod, helptextgsfilter,
      fdat.id, pdat.id,
      gutcobjects, gutcheaders
    )

  return(dat)


}

dat<-load_data()

##define app
app<-shinyApp(

ui = shinyUI(
  fluidPage(
 # tags$head(includeScript("google-analytics.js")),
  navbarPage("MCF10A Portal",
    tabPanel("About",
      titlePanel("MCF10A Portal"),
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
            column(3, selectInput("chemical_de", "CRCGN Chemical:", dat$chemicals)),
            column(2, checkboxInput("landmark_de", "Landmark only", value =FALSE)),
            column(2, selectInput("summarizefunc_de", "Summarization:",dat$summarizefuncs,
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
            column(3, selectInput("gene_de", "Gene symbol:", dat$genes)),
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
          column(3, selectInput("chemical_gs", "CRCGN Chemical:", dat$chemicals)),
          column(2, 
              selectInputWithTooltip(inputId = "gsname_gs", label = "Gene set name", 
              choices = names(dat$gsnames), bId = "Bgsname", helptext =dat$helptextgsname)
            ),
          column(2,           
            selectInputWithTooltip(inputId = "gsmethod_gs", label = "Projection method", 
              choices = dat$gsmethods, bId = "Bgsmethod", helptext =dat$helptextgsmethod)
          ),
          column(2, selectInput("summarize_gs", "Sort by:", dat$gssort,selected = "median"))
          )
        ),
      DT::dataTableOutput("gsproj_result")
    ),

     tabPanel("Heatmap (interactive)",
       fluidRow(
         column(3, 
               selectInputWithTooltip(inputId = "gctfile", label = "Dataset", 
               choices = dat$gctfiles, bId = "Bgsnamegct", helptext =dat$helptextgsname)
           ),
         column(2,
               selectInputWithTooltip(inputId = "gctmethod", label = "Projection method", 
               choices = dat$gctmethods, bId = "Bgsmethodgct", helptext =dat$helptextgsmethod)
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
               choices = dat$gctfiles, bId = "Bgsnamegctstatic", helptext =dat$helptextgsname)
               ),
         column(2, 
               selectInputWithTooltip(inputId = "gsmethod_heatmap", label = "Projection method", 
               choices = dat$gctmethods, bId = "Bgsmethodgctstatic", helptext =dat$helptextgsmethod)
           ), 
         column(3, 
           selectInputWithTooltip(inputId = "filteropt", label = "Column Filter", 
           choices = dat$filteropts, 
           selected = "signal strength > 6", bId = "Bgsfilterstatic", helptext = dat$helptextgsfilter)
           )
         ),
         plotOutput("heatmap_result")

       )
      ) ,

     tabPanel("Connectivity",
       fluidPage(
         fluidRow(         
          column(3, selectInput("chemical_gutc", "CRCGN Chemical:", dat$chemicals)),
          column(3,
            selectInputWithTooltip(inputId = "header_gutc", label = "Dataset", 
              choices = dat$gutcheaders, bId = "Bgutc", helptext =dat$helptextgutc)
            ),
          column(2, selectInput("summarize_gutc", "Sort by:", dat$gssort, selected = "median"))
         ),
         dataTableOutput("chemical_description_gutc"),
         tags$style(type="text/css", '#chemical_description_de tfoot {display:none;}'),

         ##insert empty space
         fluidRow(column(width = 1, offset = 0, style='padding:10px;')),

         dataTableOutput("gutc_result")

       )
     ), 
    tabPanel(HTML("</a></li><li><a href=\"../\">Home"))


  ))),

server = shinyServer(function(input, output, session) {

  output$chemannot_result<-DT::renderDataTable({data.table.round(dat$chemannot)},
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
    data.table.round(get_de_by_gene_table(input$gene_de, dat$deeset,
      cols<-c("sig_id", "Chemical.name", "CAS", "pert_idose", 
        "dose_rank", "distil_ss", "ss_rank")))
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
    get_de_by_gene_hist(input$gene_de, dat$deeset)
    }, width = 1000, height = 300)

  output$gene_de_lm<-renderText({
    get_de_by_gene_lm(input$gene_de, dat$deeset, values = c(1, 0))
    })

  output$gsproj_result <- DT::renderDataTable({
      currgsname<-input$gsname_gs
      eset<-get_gsproj(input$chemical_gs, dat$gslist, dat$chemannot, input$gsname_gs, input$gsmethod_gs)
      res<-summarize_gsproj(eset, order.col = input$summarize_gs)
      res<-data.frame(res)
      res<-res[, setdiff(colnames(res), "rowIDs")]
      #return hyperlink to MSigDB genesets
      if(currgsname %in% c("Hallmark", "C2"))
        res$genesets<-sapply(as.character(res$genesets), get_geneset_link)
      return(res)
    }, 
    extensions = 'Buttons',
    server = FALSE,
    escape = FALSE#,
    # options = list (dom = 'T<"clear">Blfrtip',
    #   autoWidth = FALSE,
    #   columnDefs = list(list(width = '100px', targets = list(2), className = 'dt-center',
    #       render = JS(
    #         "function(data, type, row, meta) {",
    #         "return type === 'display' && data.length > 30 ?",
    #         "'<span title=\"' + data + '\">' + data.substr(0, 30) + '..</span>' : data;",
    #         "}")), 
    #   list(targets = list(1), visible = FALSE)
    #   ),
    #   deferRender=TRUE,
    #   scrollCollapse=TRUE,
    #   pageLength = 10, lengthMenu = c(10,50,100,200,1000),
    #   buttons=c('copy','csv','print')
    # )
  )

  output$chemical_description_de<-DT::renderDataTable({
    get_chemical_description(input = input$chemical_de, 
      tab = dat$chemannot,
      cols.keep = c("BUID", "Chemical Name", "CAS", "Broad_external_Id", "Breast Carcinogenicity"))

    }, options = list(dom = ''))

  output$chemical_description_gutc<-DT::renderDataTable({
    chem_desc_table<-get_chemical_description(input = input$chemical_gutc, 
      tab = dat$chemannot,
      cols.keep = c("BUID", "Chemical Name", "CAS", "Broad_external_Id", "Breast Carcinogenicity"))
   # colnames(chem_desc_table)<-c("BUID", "Chemical Name", "CAS", "Broad_external_Id", "Breast Carcinogenicity")
    return(chem_desc_table)
    }, options = list(dom = ''))


  output$result_de<-DT::renderDataTable({
      de_table<-get_de(input$chemical_de, tab=dat$chemannot, 
        eset = dat$deeset, 
        landmark = input$landmark_de, 
        do.scorecutoff = "score" %in% input$filterbyinput_de, 
        scorecutoff = c(input$range_de[1], input$range_de[2]), 
        do.nmarkers = "number" %in% input$filterbyinput_de, 
        nmarkers = c(input$numberthresleft_de, input$numberthresright_de),
        summarize.func = input$summarizefunc_de,
        fdat.id = dat$fdat.id, 
        pdat.id = dat$pdat.id,
        landmark.values = c(1, 0))
      de_table$pr_gene_symbol<-sapply(as.character(de_table$pr_gene_symbol), get_genecard_link)

      cn<-colnames(de_table)
      cn[cn %in% "pr_id"]<-"Affy ID"
      cn[cn %in% "pr_gene_symbol"]<-"Gene Symbol"
      cn[cn %in% "pr_is_lmark"]<-"Is Landmark Gene"
      cn[cn %in% "summaryscore"]<-"Summary Mod-Zscore"
      cn[cn %in% "reactome.category"]<-"Reactome Category"
      colnames(de_table)<-cn

      return(de_table)
    },
    escape = FALSE,
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
      outfileheader<-get_heatmap_eset(ds = input$gsfile_heatmap, dsmap = dat$dsmap, method = input$gsmethod_heatmap)
      outfile<-paste(dat$heatmapdir, "/", outfileheader, sep = "")
      eset<-readRDS(outfile)
      eset<-filtereset(eset, input$filteropt, dat$chemannot)

      cns<-as.character(sapply(as.character(eset$Chemical.name), 
        function(i) subset_names(i,25)))
      cns<-make.unique(cns) #in case truncation produces non-unique ids
      colnames(eset)<-paste(cns, "_", eset$pert_dose_RANK, sep = "")
      hc<-clust_eset(eset)

      col_legend_names<-c("breast_carcinogen_any", "Genotoxicity")

      col_legend<-list(list(col_breaks = c("POSITIVE", "NEGATIVE"), 
        col_values = sapply(c("orange", "green"), to.hex),
        col_labels = c("POSITIVE", "NEGATIVE")),
       list(col_breaks = c("POSITIVE", "NEGATIVE"), 
        col_values = sapply(c("orange", "green"), to.hex),
        col_labels = c("POSITIVE", "NEGATIVE"))
      )
      names(col_legend)<-col_legend_names

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
        col_lab = col_legend_names, 
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
      get_morpheus_link(url =get_heatmap_gct(ds=input$gctfile, dsmap = dat$dsmap, 
          method = input$gctmethod), 
        domain = dat$domain),
        '">', 'click here to open in new tab', '</a>'), sep = "")
      })

  output$morpheus_result_embedded<-renderText({
      paste(c('<iframe width ="1200" height ="680" src="',
      get_morpheus_link(url = get_heatmap_gct(ds=input$gctfile, dsmap = dat$dsmap, 
          method = input$gctmethod), 
        domain = dat$domain),
        '">', '</iframe>'), sep = "")
      })

  output$gutc_result<-DT::renderDataTable({
      get_gutc(input = input$chemical_gutc, 
        tab = dat$chemannot, 
        header = input$header_gutc, 
        datlist = dat$gutcobjects,
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


