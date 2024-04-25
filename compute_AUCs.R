
# Code for computing the AUPR score and the AUROC score
# The 'trapz' function is part of the R package 'pracma'.


compute_AUPR_NEW <- function(True_Matrix, Aposteriori_Matrix){

n_nodes = dim(True_Matrix)[1]

for (i in 1:n_nodes){True_Matrix[i,i]=-1}

n_non_edges   = length(which(True_Matrix==0))
n_edges       = length(which(True_Matrix==1))

Aposterioris   = NULL

for (i in 1:n_nodes){
for (j in 1:n_nodes){
if (i!=j){Aposterioris = c(Aposterioris, Aposteriori_Matrix[i,j])}
}
}

Aposterioris = sort(Aposterioris,decreasing=TRUE)

APO_values = Aposterioris[1]

for (i in 2:length(Aposterioris)){
if (Aposterioris[i]==APO_values[length(APO_values)]){}
else
  {APO_values = c(APO_values,Aposterioris[i])}
}

TPS = NULL
FPS = NULL

MATRIX = matrix(0,n_nodes,n_nodes)

for (i in 1:length(APO_values)){

indicis = which(Aposteriori_Matrix>=APO_values[i])
MATRIX[indicis] = 1

TP = length(which(MATRIX==1 & True_Matrix==1))
FP = length(which(MATRIX==1 & True_Matrix==0)) 

TPS = c(TPS,TP)
FPS = c(FPS,FP)

}

TPS = c(TPS,n_edges)
FPS = c(FPS,n_nodes^2-n_nodes -n_edges )

i_lag = 0

for (j in 2:length(TPS)){

    i = j + i_lag

    if ((TPS[i]-TPS[i-1])>1){

        NEW_TPS = NULL
        NEW_FPS = NULL
        
        for (x in 1:(TPS[i]-TPS[i-1]-1)){
            skew    = (FPS[i]-FPS[i-1])/(TPS[i]-TPS[i-1])

            NEW_TPS = c(NEW_TPS,TPS[i-1]+x)
            NEW_FPS = c(NEW_FPS,FPS[i-1]+ skew*x)       
	}
        
        TPS = c(TPS[1:(i-1)],NEW_TPS,TPS[i:length(TPS)])
        FPS = c(FPS[1:(i-1)],NEW_FPS,FPS[i:length(FPS)])

        i_lag = i_lag + length(NEW_TPS)

}
}

PRECISION = TPS[1:length(TPS)]/(TPS[1:length(TPS)]+FPS[1:length(TPS)])

RECALL    =  TPS[1:length(TPS)]/n_edges

PRECISION = c(PRECISION[1],PRECISION)
RECALL    = c(0,RECALL)

auroc = trapz(RECALL,PRECISION)

return(auroc)

}

###########################################################
###########################################################

compute_AUROC_NEW <- function(True_Matrix, Aposteriori_Matrix){

n_nodes = dim(True_Matrix)[1]

for (i in 1:n_nodes){True_Matrix[i,i]=-1}

n_non_edges   = length(which(True_Matrix==0))
n_edges       = length(which(True_Matrix==1))

Aposterioris   = NULL

for (i in 1:n_nodes){
for (j in 1:n_nodes){
if (i!=j){Aposterioris = c(Aposterioris, Aposteriori_Matrix[i,j])}
}
}

Aposterioris = sort(Aposterioris,decreasing=TRUE)

APO_values = Aposterioris[1]

for (i in 2:length(Aposterioris)){
if (Aposterioris[i]==APO_values[length(APO_values)]){}
else
  {APO_values = c(APO_values,Aposterioris[i])}
}

ROC_x = 0
ROC_y = 0

MATRIX = matrix(0,n_nodes,n_nodes)

for (i in 1:length(APO_values)){

indicis = which(Aposteriori_Matrix>=APO_values[i])
MATRIX[indicis] = 1

TP = length(which(MATRIX==1 & True_Matrix==1))
TN = length(which(MATRIX==0 & True_Matrix==0)) 

Sensitivity = TP/n_edges;
inv_Specif  = 1 - (TN/n_non_edges)

ROC_y = c(ROC_y, Sensitivity)
ROC_x = c(ROC_x, inv_Specif)
}

ROC_x = c(ROC_x,1)
ROC_y = c(ROC_y,1)

auroc = trapz(ROC_x,ROC_y)

return(auroc)
}
