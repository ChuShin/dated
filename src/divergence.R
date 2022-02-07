#R program to carry out Ks analysis. 
#The data are first trimmed to remove values less than 0.001,
#which could generate spurious frequency peaks due to rounding, and zeros.
#The values are "log"-ed and the frequency distribution is obtained in intervals of 0.1
#in log(Ks). Mclust is called to calculate gaussian mixture models. The number of components G1
#can be varied. This model with all its separate components is plotted. The chi-square statistic
#is calculated for the goodness of fit, for log(Ks) up to 3 i.e. Ks upto about 20 as in the plots.
#This avoids the spurious peak near Ks=70 often seen.

#Places where ************* is, have parameters that might need to be changed.
#1 filename
#2 G1 and G2
#3 Upper limit of frequency in plots

library(mclust)
rm(list=ls(all=TRUE))# remove all objects from old calculations

i=2 # number of column of data to be analysed *******************
# Read values from a tab-delimited text file
data <- read.table("./data/ab_ks.txt",fill=TRUE,sep='\t',header=FALSE)
list_data=as.list(data) #converts data to a 'list'
z=list_data[i]  # picks out the ith column of data including the header
d=z[[1]]        # picks out the data values
name=names(z)   # picks out the column header (name of data set)
subdata=log(subset(d,d>=0.0001))# removes unwanted values and takes logs
#sortsubdata=sort(subdata,decreasing=FALSE)
ilo=-200 # ********************* keep these fixed
ihi=59  
br=0.1*(-200:260)# steps of 0.1 for histogram # these limits have to include all the data after subsetting.
freq=hist(subdata,br,plot=FALSE)#get frequency distribution
#freq$mids                         #these are the midpoints
#freq$counts                       #these are the counts
graph_filename=paste("CS_AB",".pdf",sep="")# sets filename *********************
# (needs to incremented for each modification until old files are deleted)
G1=2 #***************** set these to 1 initially until mclust gives the value, it can be changed later
G2=2 #****************
mnbic<-mclustBIC(subdata,G1:G2,modelName="V")# does gaussian model fitting (G=# components)
datasummary <- summary(mnbic,subdata)#calculate new data summary (includes all the model parameters)
# Define colours
plot_colors <- c("forestgreen","red","green","blue","orange","cyan","magenta","brown","pink","black")

# Define plot line type
plot_lty <- c(1,1,1,1)
# Start PDF device driver to save output to figure.pdf
pdf(file=graph_filename, height=3.5, width=5)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(4.2, 3.8, 1,1),mgp=c(2,0.6,0))
y_upper=1000 #Upper limit for frequency in plot ************************* make it larger than max(freq$counts[valid]))
plot(exp(freq$mids),freq$counts,xlim=c(0.01,20),ylim=c(0,y_upper), log="x",type="h",lwd=1.0,lend=1,
     axes=F, ann=T, xlab="",xaxs="i",yaxs="i",
     ylab="", cex.lab=0.6,lty=plot_lty[1],col=plot_colors[1])
name="CS subgenome divergence analysis"
title(main=name,font.main=4)
#lend=1 is needed,lend=2 adds offset to column heights
# Make x axis with tick marks
#axis(1, las=1,cex.axis=0.4,c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,.02,.05,.1,.2,.5,1,2,5,10,20))
axis(1, las=1,cex.axis=0.4,c(0.01,.02,.05,.1,.2,.5,1,2,.5,1,2,5,10,20))
# Plot y axis with horizontal labels 
axis(2, las=2, cex.axis=0.4)

n=length(subdata) # non-zero data
m=datasummary$parameters$mean # pick out means
v=datasummary$parameters$variance$sigmasq # pick out variances
p=datasummary$parameters$pro #pick out probabilities

# Create box around plot
box()
j<- seq(ilo,ihi)
expected<-list()
for(i in 1:G1){
  expected[[i]] <- n*p[i]*(pnorm((((j+1)/10)-m[[i]])/sqrt(v[i]))-pnorm(((j/10)-m[[i]])/sqrt(v[i])))}
output <- matrix(unlist(expected),ncol=ihi-ilo+1, byrow = TRUE)
tot_expected=colSums(output)
valid=tot_expected>10 & j<30
frq=freq$counts[valid]
expd=tot_expected[valid]
chisq=sum(((frq-expd)^2)/expd)
tails=tot_expected<=10 & j<30
freq_tails=sum(freq$counts[tails])
exp_freq_tails=sum(tot_expected[tails])
message("tails expected frequency = ",exp_freq_tails)
chisq=chisq+(freq_tails-exp_freq_tails)^2/exp_freq_tails
df1=length(frq)+1-3*(G1-1)# assumes 1 component is out of range (i.e. near spurious peak)
df2=length(frq)+1-3*G1    # assumes no component is out of range 
p1=pchisq(chisq,df1,lower.tail=FALSE)
message("p1 = ",p1)
p2=pchisq(chisq,df2,lower.tail=FALSE)
message("If the variable m has a value greater than or near 4 use p2 else use p1
for the p-value for the chi-square goodness of fit test.")
message("p2 = ",p2)
#plot histogram

#This is left here so you can remove the # to check the expected frequencies
#lines(exp(freq$mids),tot_expected,type="h",lty=plot_lty[1],lwd=2,col=rgb(0,0,0),lend=1)
#lines(exp(freq$mids),freq$counts,type="h",lty=plot_lty[1],col=plot_colors[1],lend=1)

x<- seq(-10,5,by=0.01)

y<- rep(NA,1501*G1)
dim(y)<-c(1501,G1)
for(i in 1:G1){y[,i]=dnorm(x,mean=m[i],sd=sqrt(v[i]))*p[i]}
yt=rowSums(y)*0.1*n #calculates normalised total model from model components

#plot lines for components
for(i in 1:G1){lines(exp(x), y[,i]*0.1*n, type="l", lty=5, col=plot_colors[i+2])}

#plot line for total
lines(exp(x), yt, type="l", lty=plot_lty[1], col=plot_colors[2])

# Turn off device driver (to flush output to PDF)
dev.off()

# Restore default margins
par(mar=c(5, 4, 4, 2) + 0.1)
