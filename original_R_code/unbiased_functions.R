
##############################################################################################################
feature_selection <- function(dat,lab){
    require(limma)
    
    design <- model.matrix(~1+lab)
    fit <- lmFit(dat, design)
    fit <- eBayes(fit)
    top_features <- topTable(fit, number=nrow(dat), adjust = "fdr",coef='labPR')
    
    #using a wilcoxon test
    res <- sapply(1:nrow(dat),
                  function(x,dat,lab)wilcox.test(dat[x,lab==unique(lab)[1]],
                                                 dat[x,lab==unique(lab)[2]])$p.value,
                  dat,
                  lab)
    wilcox_res <- data.frame(Feature=rownames(dat),
                             p.value=res,
                             adj.p.Val=p.adjust(res,
                                                method = 'BH',
                                                n = length(res)))
    wilcox_res <- wilcox_res[order(wilcox_res$p.value),]
    rownames(wilcox_res) <- NULL
    return(list(t_test=top_features,
                wilcox=wilcox_res))
}


#extracts the values for NN and interaction analysis 
extractSimpleValues <- function(mat,remove_self=T){
    big_mat <- array(unlist(mat), dim = c(dim(mat[[1]]), length(mat)))
    phenos <- colnames(mat[[1]])
    combinations <- expand.grid(seq(length(phenos)),seq(length(phenos)))
    if (remove_self){
        #remove combinations where both values are the sames (nearest neighbor with itself doesn't make sense)
        combinations <- combinations[combinations[,1]!=combinations[,2],]
    }
    collapsed_measurements <- t(apply(combinations,1,function(x,bm)bm[x[2],x[1],],big_mat))
    combinations <- apply(combinations,2,function(x,p)p[x],phenos)
    rownames(collapsed_measurements) <- paste(combinations[,1],' -> ',combinations[,2])
    colnames(collapsed_measurements) <- names(mat)
    collapsed_measurements[is.nan(collapsed_measurements)] <- 0
    return(collapsed_measurements)
}

getPaired <- function(nams){
    tab <- table(nams)
    tab <- names(tab)[tab==2]
    return(nams%in%tab)
}

extractRatios <- function(mat,nam){
    if (length(grep(' -> ',rownames(mat),fixed=T)) == nrow(mat)){
        nams <- t(sapply(strsplit(sapply(strsplit(rownames(mat),' - '),function(x)x[2]),' -> '),function(x)x))
        nams[,1] <- sub('. $','',nams[,1])
        paired <- getPaired(paste(nams[,1],nams[,2]))
        nams <- nams[paired,]
        COUNTS <- F
    }else{
        COUNTS <- T
        nams <- sub('[+-]$','',sapply(strsplit(rownames(mat),' - '),function(x)x[2]))
        paired <- getPaired(nams)
        nams <- nams[paired]
    }
    if (sum(paired)==0){
        ratios <- matrix(nrow=0,ncol=ncol(mat))
        colnames(ratios) <- colnames(mat)
    }else{
        mat <- mat[paired,]
        
        #get the ratios
        num_pairs <- seq(nrow(mat)/2)
        indices <- sort(rep(num_pairs,2))
        ratios <- t(sapply(num_pairs,
                           function(x,indices,mat)
                               mat[grep(x,indices)[1],]/mat[grep(x,indices)[2],],
                           indices,
                           mat))
        #if a ratio is NA that means both were 0 so the ratio is 1
        ratios[is.nan(ratios)] <- 1
        #if a ratio is inf it means that the denominator was 0 -> in this case we just set it to the maximum value of the dataset
        ratios[is.infinite(ratios)] <- max(ratios[!is.infinite(ratios)])
        
        #log2 to get a nicer behavior
        ratios <- log2(ratios)
        if (COUNTS){
            rownames(ratios) <- paste('ratio -',unique(paste0(nams,'+/-')))
        }else{
            rownames(ratios) <- paste('ratio -',unique(paste0(nams[,1],'+/- -> ',nams[,2])))
        }    
        rownames(ratios) <- paste(nam,'-',rownames(ratios))
    }
    return(ratios)
}

extract_interaction_features <- function(interactions,nam){
    f_interactions <- extractSimpleValues(interactions)
    rownames(f_interactions) <- paste(nam,'-',rownames(f_interactions))
    f_int_ratios <- extractRatios(mat=f_interactions,nam)
    dat <- rbind(f_interactions,
                 f_int_ratios)
    return(dat)
}



collapse_marker_inter <- function(x,marker){
    coords <- grep(marker,colnames(x$total),fixed = T)
    #collapse
    x$total <- cbind(x$total[,-coords],rowSums(x$total[,coords]))
    x$total <- rbind(x$total[-coords,],colSums(x$total[coords,]))
    x$nums <- c(x$nums[-coords],sum(x$nums[coords]))
    
    #clean up the name
    rownames(x$total)[nrow(x$total)] <- colnames(x$total)[ncol(x$total)] <- names(x$nums)[length(x$nums)] <- sub(' [^ ]+$','',marker)
    
    #fix the ordering so they are ordered aplhabetically again
    x$total <- x$total[order(rownames(x$total)),order(colnames(x$total))]
    x$nums <- x$nums[order(names(x$nums))]
    return(x)
}

get_collapsed_interactions <- function(selector,marker_combos,interactions){
    current_markers <- marker_combos[selector]
    inter <- interactions
    for (marker in current_markers){
        inter <- lapply(inter,collapse_marker_inter,marker)
    }
    collapsed_interactions <- lapply(inter,function(x)sweep(x$total,2,x$nums,'/'))
    return(collapsed_interactions)
}

extract_interaction_combinations <- function(interactions){
    all_markers <- colnames(interactions[[1]]$total)
    marker_combos <- table(sub('.$','',all_markers))
    marker_combos <- names(marker_combos)[marker_combos>1]
    grid <- as.matrix(expand.grid(lapply(1:length(marker_combos),function(x)c(T,F))))
    colnames(grid) <- marker_combos
    inter <- apply(grid,1,get_collapsed_interactions,marker_combos,interactions)
    return(inter)
}


collapse_marker_nn <- function(set,marker){
    for (idx in 1:length(set)){
        class <- as.character(set[[idx]]$ppp$marks)
        class[grep(marker,class,fixed=T)]<- sub(' [^ ]+$','',marker)
        set[[idx]]$ppp$marks <- as.factor(class)
    }
    return(set)
}

get_collapsed_nn <- function(selector,marker_combos,ppp){
    current_markers <- marker_combos[selector]
    set <- ppp
    for (marker in current_markers){
        set <- lapply(set,collapse_marker_nn,marker)
    }
    nn <- generate_NN(set)
    nn <- lapply(nn,function(x)x$means)
    return(nn)
}

extract_nn_combinations <- function(ppp){
    all_markers <- levels(ppp[[1]][[1]]$ppp$marks)
    marker_combos <- table(sub('.$','',all_markers))
    marker_combos <- names(marker_combos)[marker_combos>1]
    grid <- as.matrix(expand.grid(lapply(1:length(marker_combos),function(x)c(T,F))))
    colnames(grid) <- marker_combos
    nn <- apply(grid,1,get_collapsed_nn,marker_combos,ppp)
    return(nn)
}


collapse_combo_marker <- function(x,cnts){
    mrk <- cnts[rownames(cnts)==x,]
    if (is.matrix(mrk)){
        mrk <- colSums(mrk)
    }
    return(mrk)
}

#removes a marker and get all the coutns based on the collapsed marker
count_combo <- function(drop_marker,cnts){
    drop_marker <- names(drop_marker)[drop_marker]
    for (mark in drop_marker){
        rownames(cnts) <- sub(mark,'',rownames(cnts))
    }
    cnts <- t(sapply(unique(rownames(cnts)),collapse_combo_marker,cnts))
    return(cnts)
}

extract_count_combinations <- function(cnts,markers=c(' PD1.',' PDL1.')){
    all_markers <- rownames(cnts)
    grid <- as.matrix(expand.grid(lapply(1:length(markers),function(x)c(T,F))))
    colnames(grid) <- markers
    all_cnts <- apply(grid,1,count_combo,cnts)
    all_cnts <- do.call(rbind, all_cnts)
    all_cnts <- all_cnts[!duplicated(rownames(all_cnts)),]
    all_cnts <- all_cnts[order(rownames(all_cnts),decreasing = T),]
    return(all_cnts)
}

extract_count_features <- function(f_counts,nam){
    rownames(f_counts) <- paste(nam,'-',rownames(f_counts))
    f_count_ratios <- extractRatios(mat=f_counts,nam)
    dat <- rbind(f_counts,f_count_ratios)
    return(dat)
}





