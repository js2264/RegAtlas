get.ids <- function() {
    library("biomaRt")
    listMarts()
    ensembl = useDataset("celegans_gene_ensembl", mart = useMart("ensembl"))
    attributes <- listAttributes(ensembl)
    ids <- getBM(attributes = attributes[match(c('Gene name', 'WormBase Gene ID', 'NCBI gene ID', 'WikiGene name'), attributes$description),]$name, mart = ensembl)
    return(ids)
}