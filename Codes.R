################################################################
##                                                            ##
##  The codes demonstrate the "Example 1" in the paper "The   ##
##  Markov Chain Monte Carlo Revolution" by Persi Diaconis.   ##
##                                                            ##
##  To show the whole process, no external R packages are     ##
##  used, but it can be slower than using packages for MCMC.  ##
##  It may take a while to finish the whole computation. You  ##
##  can change the following 2 variables to change the load   ##
##  of computation.                                           ##
##                                                            ##
################################################################

tests=10 # Set the total number of tests. Each test is an independent MCMC process.
steps=30000 # Set the steps of each MCMC process.

set.seed(0) # Set the random seed. You can change it to what you want.

##
## Input files:
##   Hamlet.txt
##   War_and_Peace.txt
##
## Ouput files:
##   Results.txt
##   History_LOG.txt
##

##########################################
##
## Read "War and Peace" and "Hamlet"
##
##########################################

lns=readLines("War_and_Peace.txt")
sl=strsplit(lns,"")  # Decomposition

htxt=readLines("Hamlet.txt")
ht=strsplit(htxt,"")[[1]]  # Decomposition

long=3207934 # The number of characters in War_and_Peace

txt=matrix(nrow =long)
n=0  # Position pointer
for(i in 1:length(sl)){
  if(length(sl[[i]]!=0)){
    for(j in 1:length(sl[[i]])){
      n=n+1
      txt[n]=sl[[i]][j]
    }
    n=n+1
    txt[n]="\n"
  }
}
## Now, "Hamlet.txt" is stored in 'ht', 
## "War_and_Peace.txt" is stored in 'txt'.

##########################################
##
## Write the functions we need later
##
##########################################

as.n=function(c){
  ## This function is to transform characters to numbers
  ## Input: chr vector
  ## Output: numaric vector
  ## A,B,C,...,Z,a,b,c,...,z,Space  ---> 1:53
  n=matrix(nrow=length(c))
  for(i in 1:length(c)){
    ra=as.numeric(charToRaw(c[i]))
    if(ra>64 & ra<91){    # A-Z ASCII is 65:90
      ra=ra-64
    }
    else{
      if(ra>96 & ra<123){   #a-z ASCII is 97:122
        ra=ra-70
      }
      else{
        ra=ifelse(ra==32,53,0)  # space
      }
    }
    n[i]=ra
  }
  return(n)
}

as.c=function(n){
  ## This is the function to transform numbers back to characters,
  ## it is the inverse of function as.n().
  c=matrix(nrow=length(n))
  for(i in 1:length(n)){
    ra=n[i]
    if(ra>0 & ra<27){    # A-Z is 65:90
      ra=ra+64
    }
    else{
      if(ra>26 & ra<53){   #a-z is 97:122
        ra=ra+70
      }
      else{
        ra=ifelse(ra==53,32,0)
      }
    }
    ra=rawToChar(as.raw(ra))
    c[i]=ra
  }
  return(c)
}

as.str=function(c){
  ## To connect characters into a string
  str=""
  for(i in 1:length(c)){
    str=paste(str,c[i],sep="")
  }
  return(str)
}

de=function(t,k){
  ## To use a key to decrypt a text
  ## Input: text (number vector), key (number vector)
  ## Output: decrypted text (number vector)
  ## "Key" is an arangement of 1:53. For example, if key[1]=3, then all "1" in the original text will be replaced by "3".
  for(i in 1:length(t)){
    t[i]=k[t[i]]
  }
  return(t)
}

lk=function(t){
  ## To compute the likelihood 
  ## NOTICE: The output is a vector, whose product is the likelihood we need. 
  ##         However, we cannot compute the product here because it is too small 
  ##         and R will treat it as 0. So we output the vector and compute it 
  ##         later.
  ##
  lkl=matrix(nrow=length(t)-1)
  for(i in 1:(length(t)-1)){
    lkl[i]=tr[t[i],t[i+1]]
  }
  # not prod here because it will be 0
  return(lkl)
}

##########################################
##
## Compute the transition matrix and create
## the encrypted data.
##
##########################################

tn=as.n(txt)

## Compute the transition matrix
tr=matrix(1,nrow=53,ncol=53) #initail 1 is to avoid P=0
for(i in 1:(long-1)){
  l=tn[i]
  r=tn[i+1]
  if(l*r!=0){
    tr[l,r]=tr[l,r]+1
  }
}
for(i in 1:53){
  tr[i,]=tr[i,]/sum(tr[i,])
} 

## create encrypted data
key=sample(53)
rtkey=rep(NA,53) # This is the correct "key" vector we stored for computing correction rates.
for(i in 1:53){
  rtkey[i]=which(key==i)
}
hn=as.n(ht)
hne=de(hn,key)
as.str(as.c(hne))  ## Show the encryped data
## Write the encryped data in files 
cat(paste("The encryped data is:\n",as.str(as.c(hne)),"\n\n",sep=""),file="Results.txt")
cat(paste("The encryped data is:\n",as.str(as.c(hne)),"\n\n",sep=""),file="History_LOG.txt")

##########################################
##
## Use MCMC to decrypt
##
##########################################

## "tests" and "steps" have been set at the beginning of this file
for(i in 1:tests){  
  key=sample(53) # The initial key is a random arrangement.
  ## Write Log
  cat(paste("[Start the No.",i," test.]\n\n",sep=""),file="History_LOG.txt",append=T)
  ## Start Loop
  for(n in 1:steps){
    k0=key # Present State
    k1=key # Proposed State (before transposition)
    p=sample(53,2) # Choose the 2 random positions for transposition
    k1[p[1]]=k0[p[2]]  # Transposition
    k1[p[2]]=k0[p[1]]
    
    # compute pl(proposed)/pl(present)) :
    r=prod(lk(de(hne,k1))/lk(de(hne,k0))) 
    # hne is the number sequence representing the encrypted text
    # de(hne,k1) gives the number sequence representing the decrypted text using the key "k1" (Proposed key)
    # de(hne,k0) gives the number sequence representing the decrypted text using the key "k0" (Present key)
    # lk(de(hne,k1)) gives a sequence of values of possibility that a alphabet s_i is followed by s_(i+1)
    #               the values of possibility are from the matrix of transitions, i.e. M(x,y) in the reference paper
    #                by denotation in the reference paper it is the sequence of M(f(s_i), f(s_(i+1))) ,  i=1,2,...,n-1
    # prod(lk(de(hne,k1))/lk(de(hne,k0))) gives the product of all the possibility ratios, i.e. Pl(f*)/Pl(f) in the reference paper
    # So r is pl(proposed)/pl(present)
    
    # go if U < pl(proposed)/pl(present) :
    if(runif(1)<r){key=k1}else{key=k0}
    
    if(n%%1000 == 0){ # Write a Log every 1000 steps
      rst=as.str(as.c(de(hne,key)))
      cat(paste("In No. ",i," test, after ",n," steps, the result is:\n",rst,"\n\n",sep=""),file="History_LOG.txt",append=T)
    }
  }
  ## Transfrom the result into a character string
  rs=de(hne,key) 
  rst=as.str(as.c(rs))
  ## Write it to "Results.txt"
  cat(paste("No. ",i," test's result (after ",steps," steps):\n",rst,"\n\n",sep=""),file="Results.txt",append=T)
  ## Write Log
  cat(paste("No. ",i," test's result (after ",steps," steps):\n",rst,"\n\n",sep=""),file="History_LOG.txt",append=T)
  cat(paste("[End the No.",i," test.]\n\n\n",sep=""),file="History_LOG.txt",append=T)
  ## Also print the results in R.
  print(i)
  print(rst)
}

##############################################################
## Now, the final results of each tests are in "Results.txt",
## and the log of each 1000 steps are in "History_LOG.txt",
## in the current folder. You can check them, and you can 
## change the "total number of tests"(tests) and "steps of 
## MCMC"(steps) to see what happens. 
##############################################################
