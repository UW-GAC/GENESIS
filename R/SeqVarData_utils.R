.hasSex <- function(genoData){
    if(class(genoData) == "GenotypeData"){
        hasSex(genoData)
    }else if(class(genoData) == "SeqVarData"){
        "sex" %in% varLabels(sampleData(genoData))
    }
}

.XchromCode <- function(genoData){
    if(class(genoData) == "GenotypeData"){
        XchromCode(genoData)
    }else if(class(genoData) == "SeqVarData"){
        "X"
    }
}

.YchromCode <- function(genoData){
    if(class(genoData) == "GenotypeData"){
        YchromCode(genoData)
    }else if(class(genoData) == "SeqVarData"){
        "Y"
    }
}
