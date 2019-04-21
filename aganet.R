library(sparklyr)
library(dplyr)
library(DBI)
library(igraph)
library(parallel)
cores<-detectCores()-1
best<-0.0
gerSemMelhora<-0
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
selecaoTorneio<-function(g,geracao){
  t<-0.05
  numver<-length(g[1])
  #inicialmente irei criar um torneio entre dois cromossomos
  sel<-sample(1:length(geracao[,1]),5,replace = FALSE)
  
  if(geracao[sel[2],numver+1]>geracao[sel[1],numver+1]){
    sel1 = geracao[sel[2],]
  }else{
    sel1 = geracao[sel[1],]
  }
  sel<-sample(1:length(geracao[,1]),5,replace = FALSE)
  if(geracao[sel[2],numver+1]>geracao[sel[1],numver+1]){
    sel2 = geracao[sel[2],]
  }else{
    sel2 = geracao[sel[1],]
  }
  return(rbind(sel1,sel2))
  
  
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
vertViz<-function(g,a=-1,ant=-1){
  numVertex<-length(g[1])
  if(a>0 & a<=numVertex){
    x = which(g[a]==1)
    x<-x[which(x!=ant)]
    if(length(x)>=1){
      v1 = sample(x,1)
    }else{
      v1 = a
    }
    return(v1)
  }
  v1 <- sample(1:numVertex)
  #gera um vetor com vertices vizinhos
  for (a in 1:numVertex){
    x = which(g[a]==1)
    
    if(length(x)>=1){
      v1[a] = sample(x,1)
    }else{
      v1[a] = a
    }
    
    
  }
  return(v1)
  
}
cmpcon<-function(g,v){
  numVertex<-length(g[1])
  com<-1
  ncom<-1
  v2<-rep(-1,numVertex)
  for(i in 1:numVertex){
    pilha<-c(i)
    k<-v[i]
    while(!(k %in% pilha)){
      pilha<-c(pilha,k)
      if(v2[k]!=-1){
        ncom<-min(com,v2[k])
      }
      k<-v[k]
    }
    while(length(pilha)>0){
      v2[pilha[length(pilha)]]<-ncom
      pilha<-pilha[-length(pilha)]
    }
    
    if(ncom!=com){
      com<-com+1
      ncom<-com
    }else{
      com<-com+1
    }
    
  }
  return(v2)
}
initpop<-function(g,npop){
  numVertex<-length(g[1])
  ger<-matrix((1:npop),npop,numVertex+1,TRUE)
  # for(pop in 1:npop){
  #   v1 <- c()
  #   #v2 <- cmpcon(g,v1)
  #   #m <- modularity(g,v2)
  #   lpa<-label.propagation.community(g)
  #   lpaM<-lpa$membership
  #   for (a in 1:numVertex){
  #     x = lpaM[a]
  #     xx<- paste(x)
  #     eval(parse(text=paste("k<-lpa[",x,"]$'",xx,"'",sep = "")))
  #     if(length(k)>=1){
  #       #v1[a] = x
  #       v1[a] = sample(k,1)
  #     }else{
  #       v1[a] = a
  #     }
  #   }
  #   v2 <- cmpcon(g,v1)
  #   #m <- modularity(g,v2)
  #   m<-lpa$modularity
  #   print(m)
  #   ger[pop,]=c(v1,m)
  # }
  ger<-ger[order(-ger[,numVertex+1]),]
  for(pop in 1:npop){
    v1 <- vertViz(g)
    v2 <- cmpcon(g,v1)
    m <- modularity(g,v2)
    for (a in 1:numVertex){
      x = which(g[a]==1)
      
      if(length(x)>=1){
        v1[a] = sample(x,1)
      }else{
        v1[a] = a
      }
      
    }
    ger[pop,]=c(v1,m)
  }
  print(ger)
  
  return(ger)
}
cruza<-function(g,e1,e2,txaCC){
  numver<-length(g[1])
  for(i in 1:numver){
    if(txaCC>=runif(1)){
      e1[i]<-e2[i]
    }
  }
  com<-cmpcon(g,e1)
  e1[numver+1]<-modularity(g,com)
  return(e1)
}
muta<-function(g,e1,txaMut,txaMutMult){
  
  numver<-length(g[1])
  if(txaMutMult>runif(1)){
    for(i in 1:numver){
      if(0.5>=runif(1)){
        e1[i]<-vertViz(g,i,e1[i])
      }
    }
  }else{
    ver<-sample(1:numver,1)
    e1[ver]=vertViz(g,ver,e1[ver])
  }
  com<-cmpcon(g,e1)
  e1[numver+1]<-modularity(g,com)
  return(e1)
}
runger<-function(g,npop,geracao_antiga,txaEli=0.05,txaCross=0.8,txaCC=0.5,txaMut=0.2,txaMutMult=0.4){
  numver<-length(g[1])
  elite<-as.integer(ceiling(npop*txaEli))
  geracao_antiga<-geracao_antiga[order(-geracao_antiga[,numver+1]),]
  ger<-geracao_antiga
  for(p in 1:npop){
    sel<-selecao(g,geracao_antiga)
    #sel<-selecaoTorneio(g,geracao_antiga)
    if(txaCross>runif(1)){
      ger[p,]<-cruza(g,sel[1,],sel[2,],txaCC)
    }else{
      ger[p,]<-muta(g,sel[1,],txaMut = txaMut,txaMutMult = txaMutMult)
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
  ger<-initpop(g,npop)
  for(i in 1:nger){
    if(gerSemMelhora>2) break
    ger<-runger(g,npop,ger)
    print(ger[,numver+1])
  }
  cmp<-cmpcon(g,ger[1,])
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
