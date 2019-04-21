library(igraph)

best<-0.0
gerSemMelhora<-0
init<-function(g,npop){
  numVertex<-length(g[1])
  ger<-matrix((1:npop),npop,numVertex+1,TRUE)
  for(i in 1:npop){
    lpa<-label.propagation.community(g)
    mem<-lpa$membership
    mod<-lpa$modularity
    ger[i,]=c(mem,mod)
  }
  print(ger)
  return(ger)
  
  
}
createGraph<-function(name){
  if(name=='zachary')
    return(make_graph(name))
  return(read_graph(name,format = 'edgelist',directed=FALSE))
}
avalia<-function(g,geracao){
  numver<-length(g[1])
  ac<-0.0
  for(g in geracao[,numver+1]){
    ac <- ac+g
  }
  return(ac)
}
selecao<-function(g,geracao){
  numver<-length(g[1])
  tfit<-c()
  ac<-0.0
  for(g in geracao[,numver+1]){
    ac <- ac+g
    tfit<-c(tfit,g)
  }
  pt<-c()
  for(g in 1:length(geracao[,numver+1])){
    pt<-c(pt,tfit[g]/ac)
  }
  
  pk<-c()
  pt<-sort(pt,decreasing = TRUE)
  
  pac<-0.0
  for(p in pt){
    pac<-pac+p
    pk<-c(pk,pac)
  }
  
  prob<-c(runif(1),runif(1))
  
  sel1<-geracao[which(pk>=prob[1]),]
  sel2<-geracao[which(pk>=prob[2]),]
  if(is.null(dim(sel1))||dim(sel1)==1){
    selecao1<-sel1
  }else{
    selecao1<-sel1[1,]
  }
  if(is.null(dim(sel2))||dim(sel2)==1){
    selecao2<-sel2
  }else{
    selecao2<-sel2[1,]
    
  }
  if(selecao1[numver+1]==selecao2[numver+1]){
    if(!is.null(dim(sel2))){
      selecao2 <- sel2[2,]
    }else if(!is.null(dim(sel1))){
      selecao2 <- sel1[2,]
    }else{
      selecao2 <- geracao[2,]
    }
  }
  
  return(rbind(selecao1,selecao2))
}
push <- function(l, x) {
  assign(l, append(eval(as.name(l)), x), envir=parent.frame())
}
crossover<-function(g,ind1,ind2,numvert){
  vizitados<-rep(-1,numvert)
  novo <- ind2
  k<-1
  while(length(which(vizitados==-1))>0){
    sel <- sample(which(vizitados==-1),size = 1)
    ri1<-ind1[sel]
    ri2<-ind2[sel]
    caminho <- dfs(g,sel)$order
    vizitados[sel]<-1
    novo[sel]<-k
    for(c in caminho){
      
      if(ind1[c+1]==ri1){
        vizitados[c+1]<-1
        novo[c+1] <- k
      }else{
        
        break
      }
      
    }
    
    k<-k+1
    
    
    
  }
  
  
  #fase de refinamento.
  for(i in 1:numvert){
    viz<-neighbors(g,i)
    if(length(viz)!=0){
      escolhido<-sample(viz,size=5)
      
      count<-c()
      names<-c()
      escolhido = sort(escolhido)
      for(e in escolhido){
        if(!(novo[e] %in% names)){
          if(length(count)>=1){
            count<-push(count,1)
            names<-push(names,novo[e])
          }else{
            count<-c(1)
            names<-c(novo[e])
          }
          
        }else{
          print(names)
          count[[str(novo[e])]]<-count[[str(novo[e])]]+1
        }
      }
      
      
      novo[i]=max(count)
    }
  }
  com<-novo[-(numvert+1)]
  novo[numvert+1]<-modularity(g,com)
  
  
  
  
  return(novo)
  
}
chooseBestModel <- function(x) {
  tabulatedOutcomes <- table(x)
  sortedOutcomes <- sort(tabulatedOutcomes, decreasing=TRUE)
  mostCommonLabel <- names(sortedOutcomes)[1]
  mostCommonLabel
}
muta<-function(g,ind1,numvert){
  sel = sample(1:numvert,size = 1)
  
  
  viz<-neighbors(g,sel)
  l<-ind1[viz]
  if(length(l)>0){
    
    s<-sel
  }else{
    s<-sel
  }
  ind1[sel] = s
  
  for(i in viz){
    ind1[i] = s
  }
  com<-ind1[-(numvert+1)]
  ind1[numvert+1]<-modularity(g,com)
  return(ind1)
  
}
runger<-function(g,npop,geracao_antiga,txaEli=0.1,txaCross=0.5,txaCC=0.5,txaMut=0.2,txaMutMult=0.4){
  numver<-length(g[1])
  elite<-as.integer(ceiling(npop*txaEli))
  geracao_antiga<-geracao_antiga[order(-geracao_antiga[,numver+1]),]
  ger<-geracao_antiga
  for(p in 1:npop){
    sel<-selecao(g,geracao_antiga)
    #sel<-selecaoTorneio(g,geracao_antiga)
    if(txaCross>runif(1)){
      ger[p,]<-crossover(g,sel[1,],sel[2,],numver)
    }else{
      ger[p,]<-muta(g,sel[1,],numver)
      
    }
    
    
  }
  ger<-ger[order(ger[,numver+1]),]
  for(i in 1:elite){
    ger[i,]<-geracao_antiga[i,]
  }
  ger<-ger[order(-ger[,numver+1]),]
  av<-avalia(g,ger)
  if(av>best){
    best<-av
    gerSemMelhora<-0
  }else{
    gerSemMelhora<-gerSemMelhora+1
  }
  return(ger)
}
run<-function(name_graph,name,npop,nger){
  g<-createGraph(name_graph)
  numver<-length(g[1])
  ger<-init(g,npop)
  for(i in 1:nger){
    if(gerSemMelhora>2) break
    ger<-runger(g,npop,ger)
    print(ger[,numver+1])
  }
  cmp<-ger[1,]
  print(cmp)
  f<-file(paste("test",name,i,".txt",sep = "_"),"w")
  write(ger[1,],f)
  close(f)
  plot(g,vertex.color=cmp)
}
lista<-c("C:/Users/lucao/Documents/codigos/datasets/polbooks.csv","C:/Users/lucao/Documents/codigos/datasets/football.csv",'zachary','/home/lucas/Downloads/binary_networks/datas/200/5c.dat',"/home/lucas/Área de Trabalho/datasets/dolphins.paj")
pop<-c(100,100,50,50,50)
#pop<-c(20,30,100)

for(dt in 1:length(lista)){
  print(system.time(run(lista[dt],dt,pop[dt],100)))
}