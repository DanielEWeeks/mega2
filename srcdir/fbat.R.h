const char *fbatR = "\
# R --slave --vanilla < test.R --args f1\n\
# Rscript   --vanilla   test.R        f1\n\
\n\
get_log = function(str) {\n\
  once = 0\n\
  show = 0\n\
  lines = NULL\n\
  con = file(str, 'r', blocking=F)\n\
  while (T) {\n\
    readLine = readLines(con, n=1, ok=T, encoding='latin1')\n\
    if (length(readLine) == 0) break\n\
    if (show == 2) {\n\
      if (readLine == '') {\n\
        show = 0;\n\
      } else {\n\
        lines = c(lines, readLine)\n\
      }\n\
    } else if (show == 1 && substring(readLine, 1, 12) == '------------') {\n\
      show = 2;\n\
    } else {\n\
      if (substring(readLine, 1, 6) == 'Marker') {\n\
        if (!once) {\n\
          header = paste('Marker', substring(readLine, 7), sep=' ')\n\
          once = 1\n\
        }\n\
        show = 1\n\
      }\n\
    }\n\
  }\n\
  close(con)\n\
  return (c(list(lines, str, header)))\n\
}\n\
\n\
log_loop = function () {\n\
#  for (dir in args) {\n\
#    files = list.files(dir, '\\\\.sh\\\\.log$', full.names=T)\n\
    files = c(paste(prefix, '.sh.log', sep=''))\n\
    for (f in files) {\n\
      ans = get_log(f)\n\
      if (is.null(header)) header = ans[[3]]\n\
\n\
      lines = ans[[1]]\n\
      if (is.null(lines)) {\n\
        cat('NO FBAT data found in file', f, 'for trait', trait, sep=' ', '\\n')\n\
        cat('The FBAT results described below were not created:', '\\n')\n\
        return (lines)\n\
      }\n\
      lines_no = length(lines)\n\
\n\
      flines = strsplit(lines, '\\\\s+', perl=T)\n\
      df = data.frame(matrix(c(flines, recursive=T), nrow=lines_no, byrow=T))\n\
\n\
      fpaths = rep(c(trait, fix), lines_no)\n\
      fp = matrix(fpaths, nrow=lines_no, byrow=T)\n\
\n\
      nl = cbind(df, fp)\n\
      log = rbind(log, nl)\n\
    }\n\
#  }\n\
  for (i in 2:8) log[,i] = as.numeric(as.character(log[,i]))\n\
  hfields = c(unlist(strsplit(header, '\\\\s+', perl=T)), 'Trait', 'Chromosome')\n\
  colnames(log) = make.names(hfields)\n\
  rownames(log) = NULL\n\
  return (log)\n\
}\n\
\n\
map_loop = function () {\n\
  map = NULL\n\
#  for (dir in args) {\n\
#    files = list.files(dir, '.*\\\\.map$', full.names=T)\n\
    files = mapfl\n\
    for (f in files) {\n\
      mf = read.table(f, header=F)\n\
      map = rbind(map, mf)\n\
    }\n\
#  }\n\
  hfields = c('Marker', 'Chr', 'GenPos', 'BPPos', 'SexLink')\n\
  colnames(map) = make.names(hfields)\n\
  rownames(map) = NULL\n\
  return (map)\n\
}\n\
\n\
args = commandArgs(TRUE)\n\
##args = c('fbat.01.map', 'dup', 'fbat', '01')\n\
header = NULL\n\
log    = NULL\n\
jpg    = FALSE\n\
mapfl  = args[1]\n\
trait  = args[2]\n\
fix    = args[4]\n\
prefix = paste(args[3], fix, sep='.')\n\
\n\
go = function ()\n\
{\n\
  log = log_loop()\n\
  if (is.null(log)) {\n\
    return(invisible(log))\n\
  }\n\
  assign('foo', log, envir = .GlobalEnv)\n\
\n\
  map = map_loop()\n\
  assign('map', map, envir = .GlobalEnv)\n\
\n\
  bar = merge(log, map, by='Marker')\n\
  bar = bar[order(bar$Chr, bar$BPPos, bar$GenPos), ]\n\
  assign('bar', bar, envir = .GlobalEnv)\n\
  write.table(bar, file=paste(prefix, '.tbl', sep=''), sep='\t', row.names=F, quote=F)\n\
\n\
  Data = paste('R', prefix, sep='')\n\
  Hdr = paste('R', prefix, '.hdr', sep='')\n\
  Map = paste('R', prefix, '.map', sep='')\n\
  Out = paste('R', prefix, '.pdf', sep='')\n\
\n\
  bar$Log10p = -log10(bar$P)\n\
\n\
  if (sum(bar$GenPos) == 0) {\n\
    bb = bar[bar$Allele==1, c('Marker', 'BPPos','Log10p')]\n\
  } else {\n\
    bb = bar[bar$Allele==1, c('Marker', 'GenPos','Log10p')]\n\
  }\n\
  colnames(bb) = c('Marker', 'Position', trait)\n\
  write.table(bb, file=Data, row.names=F, col.names=T, sep=' ', quote=F)\n\
\n\
  cc = bar[bar$Allele==1, c('Chr', 'Marker', 'BPPos')]\n\
  colnames(cc) = c('Chromosome', 'Name', 'BP.p')\n\
  write.table(cc, file=Map, row.names=F, col.names=T, sep=' ', quote=F)\n\
  plotdata=split(bb, cc$Chromosome)\n\
\n\
  library(nplplot)\n\
  nocustomtracks = all(cc$BP.p == 0)\n\
  if (nocustomtracks) {\n\
    nplplot.multi(plotdata=plotdata, filename=c(Data), col=2, row=2, mode='l', output=Out, headerfiles=Hdr, lgnd='page', customtracks=FALSE)\n \
  } else {\n\
    nplplot.multi(plotdata=plotdata, filename=c(Data), col=2, row=2, mode='l', output=Out, headerfiles=Hdr, lgnd='page', customtracks=TRUE, mega2mapfile=cc)\n \
  }\n\
\n\
\n\
  if (jpg) {\n\
    zz = split(bar, bar$Chr)\n\
    zl = length(unique(bar$Chr))\n\
    zsqrt = sqrt(zl)\n\
    zx = ceiling(zsqrt)\n\
    zy = ceiling(zl/zx)\n\
    jpeg(paste('R', prefix, '.jpg', sep=''))\n\
    par(mfrow=c(zx, zy))\n\
#   par(mfrow=c(5,5), mai=c(.3,.3,.3,.15))\n\
    for (mx in zz) plot(x=mx$GenPos, y=-log10(mx$P), main=mx$Chr[1], ylim=c(0,10))\n\
    dev.off()\n\
  }\n\
}\n\
\n\
go()\n\
\n";
