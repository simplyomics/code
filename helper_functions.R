
# date <- Sys.Date()
# summary.file <- paste0('data/Summary/GEO_Series_Summary_Table.', date, '.RDS')
# 
# if (!file.exists(summary.file)) {
#   series.table <- downloadSeriesSummary()
#   saveRDS(series.table, file=paste0('~/Documents/Publications/webGEO/data/Summary/GEO_Series_Summary_Table.', date, '.RDS'))
# } else {
#   series.table <- readRDS(summary.file)
# }


downloadSeriesSummary <- function(num.series=NULL) { # , dest.dir=NULL
  
  if (is.null(num.series)) {
    geo.homepage = readLines('https://www.ncbi.nlm.nih.gov/geo/')
    idx <- grep('Series.*counts.*\\d+', geo.homepage)
    num.series <- as.numeric(str_extract(geo.homepage[idx], '\\d+'))
  }
  
  num.pages <- seq_len(as.integer(num.series/5000+1))
  
  start.acc <- seq(1, num.series - (num.series %% 5000) + 1, 5000)
  end.acc <- start.acc + 5000 - 1
  end.acc[num.pages[length(num.pages)]] <- num.series
  
  series.list <- list()
  
  for (i in num.pages) {
    
    url <- paste0('https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&zsort=date&mode=csv&page=', i, '&display=5000')
    
    if (i==max(num.pages)) {
      txt <- paste0('Summary information for all the ', formatC(num.series, format="d", big.mark=","), ' datasets has been downloaded !')
    } else {
      txt <- paste0('Summary information for ', formatC(i*5000, format="d", big.mark=","), ' datasets has been downloaded !')
    }
    
    series.list[[i]] <- read.csv(url, header = T, stringsAsFactors = F)
    
    message (txt)
    
  }
  
  # for (i in rev(num.page)) {
  #   
  #   url <- paste0('https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&zsort=date&mode=csv&page=', i, '&display=5000')
  #   target.file <- file.path(dest.dir, paste0('Series.', i, '.', Sys.Date(), '.csv'))
  #   download.file(url, destfile = target.file)
  #   
  # }
  
  series.table <- do.call(rbind, series.list)
  series.table <- data.frame(series.table, stringsAsFactors = F)
  
  return (series.table)
  
}



statsFun <- function(x) {
  stats <- data.frame(sort(table(x), decreasing = T),
                      stringsAsFactors = F)
  
  colnames(stats) <- c('Group', 'Count')
  stats$Group <- as.character(stats$Group)
  stats$Count <- as.numeric(stats$Count)
  
  return (stats)
  
}

splitFun <- function(group, count) {
  groups <- strsplit(group, '\\s*\\.+\\s*|;')[[1]]
  counts <- rep(count, length(groups))
  
  stats <- data.frame(Group=groups, Count=counts)
  
  return (stats)
}


GEOStatsFun <- function(series.table, group='Series') {
  
  if (group=='Series') {
    series.stats <- series.table %>% count(Series.Type, sort = TRUE)
  } else if (group=='Organisms') {
    series.stats <- series.table %>% count(Taxonomy, sort = TRUE)
  }
  
  colnames(series.stats) <- c('Group','Count')
  idx <- grep('\\s*\\.+\\s*|;', series.stats$Group)
  
  tmp.stats <- c()
  for (i in idx) {
    tmp.stats <- rbind(tmp.stats, splitFun(series.stats[i,1], series.stats[i,2]))
  }
  
  series.stats <- rbind(series.stats[-idx,], tmp.stats)
  
  series.stats <- series.stats %>% group_by(Group) %>% 
    summarise(Count=sum(Count)) %>% 
    arrange(desc(Count))
  
}


# Read SOFT files

readSOFTFileFun <- function(soft.file, method='readLines') { # fread
  
  if (method=='readLines') {
    txt <- readLines(soft.file)
    txt <- as.character(txt)
    txt <- sub('=' ,'`', txt)
    
    txt <- strsplit(txt, ' ` ', fixed = T)
    txt <- lapply(txt, function(x) {
      V1 = x[1];
      V2=ifelse(length(x)==2, x[2],paste(x[-1], collapse=' = '))
      return (data.frame(V1, V2, stringsAsFactors = F))
    })
    
    txt <- data.frame(do.call(rbind, txt), stringsAsFactors = F)
    
  } else if (method=='fread') {
    
    cmd <- paste("sed 's/=/`/'", soft.file, sep = ' ')
    txt <- fread(cmd=cmd, sep = '`', header = FALSE, fill = TRUE)
    
  }
  
  return (txt)
  
}


# Parse Series metadata from SOFT files
parseSOFTMetaFun <- function(accession, dest.dir, scope = 'Samples', amount = 'brief', forceReadLines=FALSE) {
  
  special.accessions <- c('GSE104650', 'GSE113401', 'GSE117271','GSE117403','GSE117748',
                          'GSE12071','GSE120716','GSE129721','GSE150170','GSE154855',
                          'GSE160341','GSE165203','GSE17312','GSE174829','GSE19465',
                          'GSE2187','GSE37625','GSE4270','GSE43963','GSE54179',
                          'GSE5662','GSE6256','GSE7692','GSE92449','GSE95664','GSE99011')
  
  multi.column.accessions <- c("GSE152295","GSE173083","GSE160310","GSE163475","GSE112050",
                               "GSE162511","GSE147847","GSE141706","GSE141330","GSE137305",
                               "GSE156161","GSE153779","GSE124638","GSE130334","GSE130335",
                               "GSE138569","GSE139384","GSE120752","GSE129236","GSE114198",
                               "GSE117383","GSE115505","GSE111000","GSE79555","GSE91040",
                               "GSE91041","GSE89624","GSE60683","GSE97801","GSE77521",
                               "GSE80494","GSE79297","GSE86566","GSE63464","GSE63013","GSE47917",
                               "GSE68387","GSE59298","GSE63048","GSE44647","GSE44675","GSE40355",
                               "GSE40820","GSE19000","GSE21731","GSE25933","GSE18017","GSE18884",
                               "GSE14843","GSE16369","GSE16370","GSE16374","GSE9807","GSE5617",
                               "GSE3697")
  
  # GSE165203, not publicly available
  
  
  soft.file <- file.path(dest.dir, paste0(accession, '.', scope, '.txt'))
  
  if (! file.exists(soft.file)) {
    downloadAccessionInfo(accession = accession, dest.dir = dest.dir, scope, amount)
  }
  
  ### may stop early
  # txt <- fread(soft.file, sep = '=', header = FALSE, fill = TRUE)
  # 
  # if (ncol(txt)==1) {
  #   txt <- fread(soft.file, sep = ' ', header = FALSE, fill = TRUE)
  #   txt$V2 <- apply(txt[, -c(1:2)], 1, function(x) paste(x[which(x!='')], collapse=' '))
  # } else {
  #   txt$V2 <- apply(txt[, -1], 1, function(x) paste(x[which(x!='')], collapse='='))
  # }
  # 
  # txt <- txt[,1:2]
  
  ### works, a little slow to combine columns
  # txt <- readLines(soft.file)
  # txt <- as.character(txt)
  # 
  # txt <- strsplit(txt, ' = ', fixed = T)
  # txt <- lapply(txt, function(x) {
  #   V1 = x[1];
  #   V2=ifelse(length(x)==2, x[2],paste(x[-1], collapse=' = '))
  #   return (data.frame(V1, V2, stringsAsFactors = F))
  # })
  # 
  # txt <- data.frame(do.call(rbind, txt), stringsAsFactors = F)
  
  ### great !!!
  
  if (accession %in% c(special.accessions,multi.column.accessions) | forceReadLines) {
    txt <- readSOFTFileFun(soft.file, method = 'readLines')
  } else {
    txt <- readSOFTFileFun(soft.file, method = 'fread')
  }
  
  if (ncol(txt)==1) {
    txt <- readSOFTFileFun(soft.file, method = 'readLines')
  }
  
  if (scope=='Samples') {
    
    #levels.fields <- c("title","geo_accession","status","submission_date","last_update_date","type")
    
    #TODO: for large files, try series matrix ?
    
    idx <- which(txt$V1=='!Sample_series_id' & txt$V2 != accession)
    if (length(idx)>0) {
      txt <- txt[-idx,]
    }
    
    filter <- which(startsWith(txt$V1, '^SAMPLE'))
    nsam <- length(filter)
    
    nfield.per.sample <- diff(c(filter,nrow(txt)+1))
    txt$V3 <- rep(1:nsam, nfield.per.sample)
    
    txt <- txt[-filter,]
    txt$V1 <- gsub('!Sample_', '', txt$V1)
    
    txt <- txt %>% group_by(V3) %>% mutate(V1=make.names(V1, unique = TRUE))
    fields <- unique(txt$V1)
    
    # idx <- grep('description', fields)
    # if (length(idx)>1) {
    #   fields <- insert(fields[-idx], idx[1], fields[idx])
    # }
    # 
    # idx <- grep('relation', fields)
    # if (length(idx)>1) {
    #   fields <- insert(fields[-idx], idx[1], fields[idx])
    # }
    
    idx <- grep(':ch1|data_row_count', fields)
    idx <- tail(idx, n=1)
    
    if (idx < length(fields)) {
      
      for (i in (idx+1):length(fields)) {
        
        field <- strsplit(fields[i], '.', fixed = T)[[1]][1]
        
        pos <- grep(field, fields)
        fields <- insert(fields[-pos], pos[1], fields[pos])
      }
    }
    
    nfield <- length(fields)
    metadata <- data.frame(V1=rep(fields, nsam), V2=NA, V3=rep(1:nsam, each=nfield))
    
    idx <- match(paste0(metadata$V3, '_', metadata$V1), paste0(txt$V3, '_', txt$V1))
    metadata$V2 <- txt$V2[idx]
    metadata <- spread(metadata, V1, V2) %>% select(all_of(fields))
    metadata <- addMetadataFun(metadata)
    
    fields.all <- colnames(metadata)
    
    metadata <- data.frame(metadata, stringsAsFactors = F)
    
    rownames(metadata) <- metadata$geo_accession
    colnames(metadata) <- fields.all
    
    # idx <- which(startsWith(txt$V1, 'description'))
    # idx
    # 
    # txt$V1 <- factor(txt$V1, levels = fields)
    # txt$V1
    # 
    # txt$V1[idx] <- gsub('\\.\\d+', '', txt$V1[idx])
    # 
    # txt <- txt %>% group_by(V3, V1) %>% 
    #   summarise(V2 = paste0(V2, collapse = "; "))
    
  } else if (scope=='Series') {
    
    nsam <- sum(txt$V1=='!Series_sample_id')
    
    if (sum(duplicated(txt$V1))>0) {
      fields <- txt$V1[-which(duplicated(txt$V1))]
    } else {
      fields <- txt$V1
    }
    
    txt$V1 <- factor(txt$V1, levels=fields)
    
    txt <- txt %>% group_by(V1) %>% 
      summarise(V2 = paste0(V2, collapse = "; "))
    
    idx <- which(txt$V1=='!Series_sample_id')
    
    txt <- txt %>% add_row(V1='!Series_sample_size', V2=as.character(nsam), .before = idx) %>% 
      filter(V1 !='^SERIES')
    
    txt$V1 <- gsub('!Series_', '', txt$V1)
    
    txt <- txt %>% column_to_rownames('V1')
    fields.all <- rownames(txt)
    accession <- txt['geo_accession',]
    
    metadata <- data.frame(t(txt), stringsAsFactors = F)
    rownames(metadata) <- accession
    colnames(metadata) <- fields.all
    
  }
  
  return (metadata)
}


addMetadataFun <- function(sampledat) {
  
  rownames(sampledat) <- sampledat$geo_accession
  
  if(length(grep('characteristics_ch',colnames(sampledat)))>0) {
    
    pd = sampledat %>%
      dplyr::select(dplyr::contains('characteristics_ch')) %>%
      dplyr::mutate(accession = rownames(.)) %>%
      # these next two lines avoid warnings due
      # to columns having different factor levels
      # (attributes).
      mutate_all(as.character) %>%
      tidyr::gather(characteristics, kvpair, -accession) %>%
      dplyr::filter(grepl(':',kvpair) & !is.na(kvpair))
    # Thx to Mike Smith (@grimbough) for this code
    # sometimes the "characteristics_ch1" fields are empty and contain no 
    # key:value pairs. spread() will fail when called on an
    # empty data_frame.  We catch this case and remove the 
    # "charactics_ch1" column instead
    if(nrow(pd)) {
      pd = dplyr::mutate(pd, characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
        tidyr::separate(kvpair, into= c('k','v'), sep=":", fill = 'right', extra = "merge") %>%
        dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
        dplyr::select(-characteristics) %>%
        dplyr::filter(!is.na(v)) %>%
        dplyr::group_by(accession,k) %>%
        dplyr::mutate(v = paste0(trimws(v), collapse = ";")) %>%
        unique() %>%
        tidyr::spread(k,v)
    } else {
      pd = pd %>% 
        dplyr::select(accession)
    }
    
    sampledat = sampledat %>% dplyr::left_join(pd,by=c('geo_accession'='accession'))
    
    ##     dplyr::mutate(characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
    ##     dplyr::filter(grepl(':',kvpair)) %>% 
    ##     tidyr::separate(kvpair, into= c('k','v'), sep=":")
    ## if(nrow(pd)>0) {
    ##     pd = pd %>% dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
    ##         dplyr::select(-characteristics) %>%
    ##         tidyr::spread(k,v)
    
  } 
  
  return (sampledat)
  
}



###
listSeriesMatrixFilesFun <- function(accession) {
  
  folder <- gsub('\\d\\d\\d$', 'nnn', accession)
  url = paste0('https://ftp.ncbi.nlm.nih.gov/geo/series/', folder, '/', 
               accession, '/matrix/')
  
  a <- read_html(url)
  files <- xml_text(xml_find_all(a,'//a/@href'))
  files <- files[which(startsWith(files, 'GSE'))]
  return(files)
}


parseSeriesMetaFun <- function(fname,AnnotGPL=FALSE,destdir=tempdir(),getGPL=TRUE,parseCharacteristics=TRUE) {
  dat <- fread(fname, sep='\t', fill = TRUE, nrows = 200, header = F, tmpdir = destdir)
  
  series_header_row_count <- sum(startsWith(dat$V1, '!Series_'))
  # In the case of ^M in the metadata (GSE781, for example), the line counts
  # for "skip" and read.table are different.
  # This next line gets the "skip" line count for use
  # below in the tmpdat reading "skip"
  sample_header_start <- which(startsWith(dat$V1, '!Sample_'))[1]
  
  samples_header_row_count <- sum(startsWith(dat$V1, '!Sample_'))
  series_table_begin_line = which(startsWith(dat$V1, '!series_matrix_table_begin'))
  
  if(length(series_table_begin_line) != 1) {
    stop("parsing failed--expected only one '!series_data_table_begin'")
  }
  #con <- fileOpen(fname)
  
  ## Read the !Series_ and !Sample_ lines
  header <- dat[1:series_header_row_count,1:2]
  header <- data.frame(header, stringsAsFactors = F)
  
  # header <- data.frame(do.call(rbind, strsplit(dat[1:series_header_row_count], "\t")),
  #                      stringsAsFactors = F)
  
  tmpdat <- dat[sample_header_start:(sample_header_start+samples_header_row_count-1),]
  tmpdat <- data.frame(tmpdat, stringsAsFactors = F)
  
  #tmpdat <- data.frame(do.call(rbind, strsplit(dat[sample_header_start:(sample_header_start+samples_header_row_count-1)], "\t")),
  #                     stringsAsFactors = F)
  
  headertmp <- t(header)
  headerdata <- rbind(data.frame(), headertmp[-1,])
  colnames(headerdata) <- sub('!Series_','',as.character(header[,1]))
  
  # headerlist <- lapply(split.default(headerdata, names(headerdata)),
  #                      function(x) {
  #                        as.character(Reduce(function (a,b) {paste(a,b,sep = "\n")}, x))
  #                      })
  # 
  # link = "https://www.ncbi.nlm.nih.gov/geo/"
  # if (!is.null(headerlist$web_link)) {
  #   link <- headerlist$web_link
  # } else if (!is.null(headerlist$geo_accession)) {
  #   link <- paste(link, 'query/acc.cgi?acc=', headerlist$geo_accession, sep="")
  # }
  # 
  # ed <- new ("MIAME",
  #            name = ifelse(is.null(headerlist$contact_name), '', headerlist$contact_name),
  #            title = headerlist$title,
  #            contact = ifelse(is.null(headerlist$contact_email), '', headerlist$contact_email),
  #            pubMedIds = ifelse(is.null(headerlist$pubmed_id), '', headerlist$pubmed_id), 
  #            abstract = ifelse(is.null(headerlist$summary), '', headerlist$summary),
  #            url = link,
  #            other = headerlist)
  
  tmptmp <- t(tmpdat)
  sampledat <- rbind(data.frame(),tmptmp[-1,])
  colnames(sampledat) <- make.unique(sub('!Sample_','',as.character(tmpdat[,1])))
  
  sampledat[['geo_accession']] <- as.character(sampledat[['geo_accession']])
  rownames(sampledat) <- sampledat[['geo_accession']]
  
  sampledat <- addMetadataFun(sampledat)
  
  ## Lots of GSEs now use "characteristics_ch1" and
  ## "characteristics_ch2" for key-value pairs of
  ## annotation. If that is the case, this simply
  ## cleans those up and transforms the keys to column
  ## names and the values to column values.
  # if(length(grep('characteristics_ch',colnames(sampledat)))>0 && parseCharacteristics) {
  #   pd = sampledat %>%
  #     dplyr::select(dplyr::contains('characteristics_ch')) %>%
  #     dplyr::mutate(accession = rownames(.)) %>%
  #     # these next two lines avoid warnings due
  #     # to columns having different factor levels
  #     # (attributes).
  #     mutate_all(as.character) %>%
  #     tidyr::gather(characteristics, kvpair, -accession) %>%
  #     dplyr::filter(grepl(':',kvpair) & !is.na(kvpair))
  #   # Thx to Mike Smith (@grimbough) for this code
  #   # sometimes the "characteristics_ch1" fields are empty and contain no 
  #   # key:value pairs. spread() will fail when called on an
  #   # empty data_frame.  We catch this case and remove the 
  #   # "charactics_ch1" column instead
  #   if(nrow(pd)) {
  #     pd = dplyr::mutate(pd, characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
  #       tidyr::separate(kvpair, into= c('k','v'), sep=":", fill = 'right', extra = "merge") %>%
  #       dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
  #       dplyr::select(-characteristics) %>%
  #       dplyr::filter(!is.na(v)) %>%
  #       dplyr::group_by(accession,k) %>%
  #       dplyr::mutate(v = paste0(trimws(v), collapse = ";")) %>%
  #       unique() %>%
  #       tidyr::spread(k,v)
  #   } else {
  #     pd = pd %>% 
  #       dplyr::select(accession)
  #   }
  #   ##     dplyr::mutate(characteristics=ifelse(grepl('_ch2',characteristics),'ch2','ch1')) %>%
  #   ##     dplyr::filter(grepl(':',kvpair)) %>% 
  #   ##     tidyr::separate(kvpair, into= c('k','v'), sep=":")
  #   ## if(nrow(pd)>0) {
  #   ##     pd = pd %>% dplyr::mutate(k = paste(k,characteristics,sep=":")) %>%
  #   ##         dplyr::select(-characteristics) %>%
  #   ##         tidyr::spread(k,v)
  #   sampledat = sampledat %>% dplyr::left_join(pd,by=c('geo_accession'='accession'))
  # }
  
  ## used to be able to use colclasses, but some SNP arrays provide only the
  ## genotypes in AA AB BB form, so need to switch it up....
  ##  colClasses <- c('character',rep('numeric',nrow(sampledat)))
  # datamat <- read_tsv(fname,quote='"',
  #                     na=c('NA','null','NULL','Null'), 
  #                     # somewhere in the past month or so, read_tsv changed
  #                     # the way it dealt with skip!!! Had to add the -1 to
  #                     # avoid the problem.
  #                     skip = series_table_begin_line,
  #                     comment = '!series_matrix_table_end',
  #                     skip_empty_rows = FALSE)
  # tmprownames = datamat[[1]]
  # # need the as.matrix for single-sample or empty GSE
  # datamat <- as.matrix(datamat[!is.na(tmprownames),-1])
  # rownames(datamat) <- tmprownames[!is.na(tmprownames)]
  # datamat <- as.matrix(datamat)
  # rownames(sampledat) <- colnames(datamat)
  # GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)])
  # ## if getGPL is FALSE, skip this and featureData is then a data.frame with no columns
  # fd = new("AnnotatedDataFrame",data=data.frame(row.names=rownames(datamat)))
  # if(getGPL) {
  #   gpl <- getGEO(GPL,AnnotGPL=AnnotGPL,destdir=destdir)
  #   vmd <- Columns(gpl)
  #   dat <- Table(gpl)
  #   ## GEO uses "TAG" instead of "ID" for SAGE GSE/GPL entries, but it is apparently
  #   ##     always the first column, so use dat[,1] instead of dat$ID
  #   ## The next line deals with the empty GSE
  #   tmpnames=character(0)
  #   if(ncol(dat)>0) {
  #     tmpnames=as.character(dat[,1])
  #   }
  #   ## Fixed bug caused by an ID being "NA" in GSE15197, for example
  #   tmpnames[is.na(tmpnames)]="NA"
  #   rownames(dat) <- make.unique(tmpnames)
  #   ## Apparently, NCBI GEO uses case-insensitive matching
  #   ## between platform IDs and series ID Refs ???
  #   dat <- dat[match(tolower(rownames(datamat)),tolower(rownames(dat))),]
  #   # Fix possibility of duplicate column names in the
  #   # GPL files; this is prevalent in the Annotation GPLs
  #   rownames(vmd) <- make.unique(colnames(Table(gpl)))
  #   colnames(dat) <- rownames(vmd)
  #   fd <- new('AnnotatedDataFrame',data=dat,varMetadata=vmd)
  # }
  # if(is.null(nrow(datamat))) {
  #   ## fix empty GSE datamatrix
  #   ## samplename stuff above does not work with
  #   ## empty GSEs, so fix here, also
  #   tmpnames <- names(datamat)
  #   rownames(sampledat) <- tmpnames
  #   datamat=matrix(nrow=0,ncol=nrow(sampledat))
  #   colnames(datamat) <- tmpnames
  # } else {
  #   ## This looks like a dangerous operation but is needed
  #   ## to deal with the fact that NCBI GEO allows case-insensitive
  #   ## matching and we need to pick one.
  #   rownames(datamat) <- rownames(dat)
  # }
  # eset <- new('ExpressionSet',
  #             phenoData=as(sampledat,'AnnotatedDataFrame'),
  #             annotation=GPL,
  #             featureData=fd,
  #             experimentData=ed,
  #             exprs=as.matrix(datamat))
  # return(list(GPL=as.character(sampledat[1,grep('platform_id',colnames(sampledat),ignore.case=TRUE)]),eset=eset))
  #sampledat <- data.frame(sampledat, stringsAsFactors = F)
  return(sampledat)
  
}


### for shiny app


parseSampleMetaShinyFun <- function(accession, series.table, dest.dir = 'data/SOFT/Samples') {
  
  scope <- 'Samples'
  soft.file <- file.path(dest.dir, paste0(accession, '.', scope, '.txt'))
  
  if (! accession %in% series.table$Accession) {
    samples.meta <- 'Invalid'
    message ('Invalid GSE accession !')
    
  } else if (paste0('Sample_Metadata.',accession,'.RDS') %in% list.files('data/Special/')) {
    message ('Large SOFT file! Return preprocessed data !')
    samples.meta <- readRDS(file.path('data/Special/', paste0('Sample_Metadata.',accession,'.RDS')))
    
  } else {
    if (! file.exists(soft.file)) {
      message ('Download SOFT file !')
      downloadAccessionInfo(scope = scope, amount = 'brief', accession = accession, dest.dir = dest.dir)
    }
    
    if (! file.exists(soft.file)) {
      message ('Failed to download soft file! Parse sample metadata from series matrix file !')
      filenames <- listSeriesMatrixFilesFun(accession)
      
      #dest.dir.matrix <- 'data/SeriesMatrix/'
      
      folder <- gsub('\\d\\d\\d$', 'nnn', accession)
      url = paste0('https://ftp.ncbi.nlm.nih.gov/geo/series/', folder, '/', 
                   accession, '/matrix/')
      
      meta.list.accession <- c()
      for (fname in filenames) {
        
        series.file.url <- paste0(url, fname)
        #series.matrix.file <- file.path(dest.dir.matrix, fname)
        #download.file(url = series.file.url, destfile = series.matrix.file)
        #series.meta <- parseSeriesMetaFun(series.matrix.file)
        series.meta <- parseSeriesMetaFun(series.file.url)
        
        meta.list.accession[[fname]] <- series.meta
      }
      
      meta.summary <- bind_rows(meta.list.accession)
      rownames(meta.summary) <- meta.summary$geo_accession
      
      fields <- colnames(meta.summary)
      
      idx <- grep(':ch1|data_row_count', fields)
      idx <- tail(idx, n=1)
      
      if (idx < length(fields)) {
        
        for (i in (idx+1):length(fields)) {
          
          field <- strsplit(fields[i], '.', fixed = T)[[1]][1]
          
          pos <- grep(field, fields)
          fields <- insert(fields[-pos], pos[1], fields[pos])
        }
      }
      
      samples.meta <- meta.summary[sort(meta.summary$geo_accession),fields]
      
    } else {
      
      if (file.info(soft.file)$size == 0) {
        message ('No sample metadata available ! SOFT file is empty !')
        samples.meta <- 'Empty'
        
      } else {
        message ('Parse metadata from SOFT file ! Read the SOFT file using "fread" !')
        samples.meta <- parseSOFTMetaFun(accession, dest.dir = dest.dir, scope = 'Samples', amount = 'brief')
        
        idx <- which(series.table$Accession==accession)
        
        if (nrow(samples.meta) != series.table$Sample.Count[idx]) {
          message ('Special character "`" detected! Read the SOFT file using "readLines" !')
          samples.meta <- parseSOFTMetaFun(accession, dest.dir = dest.dir, scope = 'Samples', amount = 'brief', forceReadLines = TRUE)
        }
      }
    }
  }
  
  return (samples.meta)
  
}





### Download accession information
# SOFT file

downloadAccessionInfo <- function(scope='Series', amount='full', accession=NULL, dest.dir=NULL) {
  ### scope: series, platform, samples, family
  ### amount: brief, quick, full, data
  ### format: SOFT (HTML, MINiML) -- Only SOFT is used
  
  if (scope=='Series') {
    targ <- 'gse'
  } else if (scope=='Platform') {
    targ <- 'gpl'
  } else if (scope=='Samples') {
    targ <- 'gsm'
  } else if (scope=='Family') {
    targ <- 'all'
  }
  
  geo.url <- paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', accession, '&targ=', targ, '&form=text&view=', amount)
  
  target.file <- file.path(dest.dir, paste0(accession, '.', scope, '.txt'))
  download.file(url = geo.url, destfile = target.file)
  
}
