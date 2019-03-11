#' Load a collection of 10X data-sets
#'
#' Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#'
#' @export
#' @param dataDirs Vector of top level cellranger output directories (the directory that contains the "raw_gene_bc_matrices" folder).
#' @param channelNames To make droplet barcodes unique across experiment, each channel needs its own unique label.  If not given, this is set numerically.
#' @param ... Extra parameters passed to \code{SoupChannel} construction function.
#' @return A SoupChannelList object containing the count tables for each 10X dataset.
#' @seealso SoupChannel SoupChannelList estimateSoup
#' @importFrom Seurat Read10X
load10X = function(dataDirs,channelNames=NULL,...){
  if(is.null(channelNames))
    channelNames = sprintf('Channel%d',seq_along(dataDirs))
  channels = list()
  for(i in seq_along(dataDirs)){
    message(sprintf("Loading data for 10X channel %s from %s",channelNames[i],dataDirs[i]))
    dataDir = dataDirs[i]
    #Get reference
    ref = "ref"
    #Load the 10X data
    tod = Read10X(file.path(dataDir,'raw_feature_bc_matrix')) #remove ref
    #Get the barcodes that cell ranger considers to contain cells
    cells = read.delim(gzfile(file.path(dataDir,'filtered_feature_bc_matrix','barcodes.tsv.gz')),sep='\t',header=FALSE) #Unzip the compressed barcode file
    cells = gsub('-1','',cells[,1])
    #Get the index in the big table
    cellIdxs = match(cells,colnames(tod))
    channels[[channelNames[i]]] = SoupChannel(tod,tod[,cellIdxs,drop=FALSE],channelName=channelNames[i],ref=ref,path=dataDir,dataType='10X',...)
  }
  channels = SoupChannelList(channels)
  return(channels)
}
