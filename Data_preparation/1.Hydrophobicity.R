setwd("/home/baptiste/Desktop/PhD/smorfs functionalization")

# calculation of hydrophobicity
library(Peptides)
library(stringr)

Annotation=read.csv("data/Annotation.csv")

avg_hydrophobicity_per_position <-function(seqs,L_max,L_min,hydrophobicity_scale){
  H <- rep(0,L_max) # hydrophobicity at each position
  N <- rep(0,L_max) # count of sequence at each position
  i=0
  for(seq in seqs){
    # print progress for every 1000
    i=i+1
    if (i%%1000==0) {
      print(i)
    }
    
    if(nchar(seq) >= L_min){
      # extend the N-terminal with x
      seq = paste(paste(rep('x',L_max),collapse=''),seq,sep='')
      seq = substr(seq,nchar(seq)-L_max+1,nchar(seq))
      
      hi = hydrophobicity(unlist(strsplit(seq,'')),hydrophobicity_scale) 
      H = H + hi 
      N = as.numeric(N) + as.numeric(hi!=0)
    }
  }
  H = H/N
  return(list(H,N))
}

L_max=100
L_min=10
hydrophobicity_scale = 'Miyazawa' 

smorfs_seq=Annotation$Peptide.seq
smorfs_seq=substr(smorfs_seq,1,nchar(smorfs_seq)-1)

smorfs_hydro=avg_hydrophobicity_per_position(smorfs_seq,L_max,L_min,hydrophobicity_scale)



L=100
skip_last=2
x=-L:-1
y=smorfs_hydro[[1]][(L_max-L+1):L_max]
lo <- loess(y[1:(L-skip_last)]~x[1:(L-skip_last)])
ylim=c(min(lo$fitted)-0.1,max(lo$fitted)+0.1)
plot(x,y,ylim=ylim,pch=16,cex=0.2,col='red',bty='l',ylab='\nAverage hydrophobicity',xlab='Position relative to C-termini\n')
lines(x, predict(lo,x), col='red', lwd=2)


# Add Uniprot

uniprot=fread("")







































