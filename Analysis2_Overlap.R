n_compare = length(q)
names(q)<-c('Sub1_Sub2','Sub3_Sub4','Sub1_Sub3')
q_overlap<-matrix(rep(0,n_compare^2),nrow=n_compare)
tmp_fdr = 0.1
for (i in 1:n_compare){
  for (j in 1:n_compare){
    q_overlap[i,j] = sum(q[[i]]<tmp_fdr & q[[j]]<tmp_fdr & !is.na(q[[i]]) & !is.na(q[[j]]))
  }
}
print(q_overlap)

path_overlap<-matrix(rep(0,n_compare^2),nrow=n_compare)
for (i in 1:n_compare){
  for (j in 1:n_compare){
    q_overlap[i,j] = sum(c2_path_msig[[i]][[1]] %in% c2_path_msig[[j]][[1]])
  }
}
print(q_overlap)

