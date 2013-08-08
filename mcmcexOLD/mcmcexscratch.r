 %%This contains a bunch of old code for creating an array of tables that is interesting

%% V,W tables

library(xtable)
load("scors.Rdata")
load("dcors.Rdata")
load("ecors.Rdata")
cnam <- colnames(scors[[2]][[1]])
rnam <- rownames(scors[[2]][[1]])
algn <- paste(paste("l",paste(rep("r",length(rnam)), collapse=""), sep=""),"", sep="")
sclbx <- 0.79

%V,W tables for T=10
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT10Vs, echo=FALSE, results='asis'>>=
T10Vs <- xtable(scors[[2]][[1]], digit=4, align=algn)
print(T10Vs,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10Ws, echo=FALSE, results='asis'>>=
T10Ws <- xtable(scors[[3]][[1]], digit=4, align=algn)
print(T10Ws,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10Vd, echo=FALSE, results='asis'>>=
T10Vd <- xtable(dcors[[2]][[1]], digit=4, align=algn)
print(T10Vd,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10Wd, echo=FALSE, results='asis'>>=
T10Wd <- xtable(dcors[[3]][[1]], digit=4, align=algn)
print(T10Wd,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10Ve, echo=FALSE, results='asis'>>=
T10Ve <- xtable(ecors[[2]][[1]], digit=4, align=algn)
print(T10Ve,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10We, echo=FALSE, results='asis'>>=
T10We <- xtable(ecors[[3]][[1]], digit=4, align=algn)
print(T10We,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=10$ for $V$ and $W$ and for the state, scaled
  disturbance, and scaled error samplers. Row and column names
  indicate the true values of $W$ and $V$ respectively for the
  simulated data. Note that the signal-to-noise ratio is constant
  across any diagonal.}
\label{tableVWT10}
\end{table}

%V,W tables for T=100
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT100Vs, echo=FALSE, results='asis'>>=
T100Vs <- xtable(scors[[2]][[2]], digit=4, align=algn)
print(T100Vs,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100Ws, echo=FALSE, results='asis'>>=
T100Ws <- xtable(scors[[3]][[2]], digit=4, align=algn)
print(T100Ws,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100Vd, echo=FALSE, results='asis'>>=
T100Vd <- xtable(dcors[[2]][[2]], digit=4, align=algn)
print(T100Vd,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100Wd, echo=FALSE, results='asis'>>=
T100Wd <- xtable(dcors[[3]][[2]], digit=4, align=algn)
print(T100Wd,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100Ve, echo=FALSE, results='asis'>>=
T100Ve <- xtable(ecors[[2]][[2]], digit=4, align=algn)
print(T100Ve,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100We, echo=FALSE, results='asis'>>=
T100We <- xtable(ecors[[3]][[2]], digit=4, align=algn)
print(T100We,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=100$ for $V$ and $W$ and for the state, scaled
  disturbance, and scaled error samplers. Row and column names
  indicate the true values of $W$ and $V$ respectively for the
  simulated data. Note that the signal-to-noise ratio is constant
  across any diagonal.}
\label{tableVWT100}
\end{table}

%V,W tables for T=1000
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000Vs, echo=FALSE, results='asis'>>=
T1000Vs <- xtable(scors[[2]][[3]], digit=4, align=algn)
print(T1000Vs,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000Ws, echo=FALSE, results='asis'>>=
T1000Ws <- xtable(scors[[3]][[3]], digit=4, align=algn)
print(T1000Ws,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000Vd, echo=FALSE, results='asis'>>=
T1000Vd <- xtable(dcors[[2]][[3]], digit=4, align=algn)
print(T1000Vd,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000Wd, echo=FALSE, results='asis'>>=
T1000Wd <- xtable(dcors[[3]][[3]], digit=4, align=algn)
print(T1000Wd,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000Ve, echo=FALSE, results='asis'>>=
T1000Ve <- xtable(ecors[[2]][[3]], digit=4, align=algn)
print(T1000Ve,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{V in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000We, echo=FALSE, results='asis'>>=
T1000We <- xtable(ecors[[3]][[3]], digit=4, align=algn)
print(T1000We,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{W in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=1000$ for $V$ and $W$ and for the state, scaled
  disturbance, and scaled error samplers. Row and column names
  indicate the true values of $W$ and $V$ respectively for the
  simulated data. Note that the signal-to-noise ratio is constant
  across any diagonal.}
\label{tableVWT1000}
\end{table}

<<corrplots, echo=FALSE>>=
library(ggplot2)
plotdat <- data.frame(corr=0, sampler=0, variable=0, T=0, V=0, W=0)
k <- 1
Ts <- c(10, 100, 1000)
Vs <- c(0.01, 0.1, 1, 10, 100, 1000)
Ws <- Vs
for(i in 1:6){
  for(j in 1:6){
    for(t in 1:3){
      plotdat[k,] <- c(scors[[2]][[t]][i,7-j], "state", "V", Ts[t], rnam[i], cnam[7-j])
      k <- k+1
      plotdat[k,] <- c(scors[[3]][[t]][i,7-j], "state", "W", Ts[t], rnam[i], cnam[7-j])
      k <- k+1
      plotdat[k,] <- c(dcors[[2]][[t]][i,7-j], "dist", "V", Ts[t], rnam[i], cnam[7-j])
      k <- k+1
      plotdat[k,] <- c(dcors[[3]][[t]][i,7-j], "dist", "W", Ts[t], rnam[i], cnam[7-j])
      k <- k+1
      plotdat[k,] <- c(ecors[[2]][[t]][i,7-j], "error", "V", Ts[t], rnam[i], cnam[7-j])
      k <- k+1
      plotdat[k,] <- c(ecors[[3]][[t]][i,7-j], "error", "W", Ts[t], rnam[i], cnam[7-j])
      k <- k+1
    }
  }
}
plotdat$corr <- as.numeric(plotdat$corr)
plotdat$T <- as.numeric(plotdat$T)
@

\begin{figure}[htb]
  \centering
<<stateplot, echo=FALSE>>=
m <- qplot(T, corr, data=plotdat[plotdat$sampler=="state",], facets=V~W, color=variable, geom="line", size=I(1), log="x", ylab="Autocorrelation", main="State Sampler")
m + scale_y_continuous( breaks=c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.45))
@
\caption{Plot of autocorrelations vs. length of time series ($T$) for
  the State Sampler. Across a diagonal, the signal-to-noise ratio
  ($W/V$) is constant. Note that when the signal-to-noise ratio is
  high, $W$ has high autocorrelation and $V$ has low autocorrelation,
  and when the signal-to-noise ratio is low, $V$ has high
  autocorrelation and $W$ has low autocorrelation. When the
  signal-to-noise ratio is near 1, both $V$ and $W$ have moderate to
  low autocorrelation. Also note that increasing $T$ seems to increase
  autocorrelation for both $W$ and $V$.}
\label{stateplot}
\end{figure}

\begin{figure}[htb]
  \centering
<<distplot, echo=FALSE>>=
m <- qplot(T, corr, data=plotdat[plotdat$sampler=="dist",], facets=V~W, color=variable, geom="line", size=I(1), log="x", ylab="Autocorrelation", main="Scaled Disturbance Sampler")
m + scale_y_continuous( breaks=c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.45))
@
\caption{Plot of autocorrelations vs. length of time series ($T$) for
  the Scaled Disturbance Sampler. Across a diagonal, the signal-to-noise ratio
  ($W/V$) is constant. Note that when the signal-to-noise ratio is
  high, both $V$ and $W$ have low autocorrelation and when the signal-to-noise
  ratio is low, both $V$ and $W$ have high autocorrelation. Also note
  that increasing $T$ seems to increase autocorrelation for both $W$
  and $V$.}
\label{distplot}
\end{figure}

\begin{figure}[htb]
  \centering
<<errorplot, echo=FALSE>>=
m <- qplot(T, corr, data=plotdat[plotdat$sampler=="error",], facets=V~W, color=variable, geom="line", size=I(1), log="x", ylab="Autocorrelation", main="Scaled Error Sampler")
m + scale_y_continuous( breaks=c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust=0.45))
@
\caption{Plot of autocorrelations vs. length of time series ($T$) for
  the Scaled Error Sampler. Across a diagonal, the signal-to-noise ratio
  ($W/V$) is constant. Note that when the signal-to-noise ratio is
  high, both $V$ and $W$ have high autocorrelation and when the signal-to-noise
  ratio is low, both $V$ and $W$ have low autocorrelation. Also note
  that increasing $T$ seems to increase autocorrelation for both $W$
  and $V$.}
\label{errorplot}
\end{figure}


\subsection{Results for the $\theta$'s}
<<maxcor, echo=FALSE>>=
mscor <- round(max(scors[[1]][[1]], scors[[1]][[2]], scors[[1]][[3]]),2)
mdcor <- round(max(dcors[[1]][[1]], dcors[[1]][[2]], dcors[[1]][[3]]),2)
mecor <- round(max(ecors[[1]][[1]], ecors[[1]][[2]], ecors[[1]][[3]]),2)
@


%Theta tables for T=10
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT10TH0s, echo=FALSE, results='asis'>>=
sTH0T10 <- scors[[1]][[1]][,,1]
rownames(sTH0T10) <- rnam
colnames(sTH0T10) <- cnam
TH0T10s <- xtable(sTH0T10, digit=4, align=algn)
print(TH0T10s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10THTs, echo=FALSE, results='asis'>>=
sTHTT10 <- scors[[1]][[1]][,,10+1]
rownames(sTHTT10) <- rnam
colnames(sTHTT10) <- cnam
THTT10s <- xtable(sTHTT10, digit=4, align=algn)
print(THTT10s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10TH0d, echo=FALSE, results='asis'>>=
dTH0T10 <- dcors[[1]][[1]][,,1]
rownames(dTH0T10) <- rnam
colnames(dTH0T10) <- cnam
TH0T10d <- xtable(dTH0T10, digit=4, align=algn)
print(TH0T10d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10THTd, echo=FALSE, results='asis'>>=
dTHTT10 <- dcors[[1]][[1]][,,10+1]
rownames(dTHTT10) <- rnam
colnames(dTHTT10) <- cnam
THTT10d <- xtable(dTHTT10, digit=4, align=algn)
print(THTT10d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10TH0e, echo=FALSE, results='asis'>>=
eTH0T10 <- ecors[[1]][[1]][,,1]
rownames(eTH0T10) <- rnam
colnames(dTH0T10) <- cnam
TH0T10e <- xtable(eTH0T10, digit=4, align=algn)
print(TH0T10e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10THTe, echo=FALSE, results='asis'>>=
eTHTT10 <- ecors[[1]][[1]][,,10+1]
rownames(eTHTT10) <- rnam
colnames(eTHTT10) <- cnam
THTT10e <- xtable(eTHTT10, digit=4, align=algn)
print(THTT10e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=10$ for $\theta_0$ and $\theta_T$ for the state, scaled
  disturbance, and scaled error samplers. Row and column names indicate
  the true values of $W$ and $V$ respectively for the simulated
  data. Note that the signal-to-noise ratio is constant across any diagonal.}
\label{tableTHT10}
\end{table}

%Theta tables for T=100
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT100TH0s, echo=FALSE, results='asis'>>=
sTH0T100 <- scors[[1]][[2]][,,1]
rownames(sTH0T100) <- rnam
colnames(sTH0T100) <- cnam
TH0T100s <- xtable(sTH0T100, digit=4, align=algn)
print(TH0T100s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100THTs, echo=FALSE, results='asis'>>=
sTHTT100 <- scors[[1]][[2]][,,100+1]
rownames(sTHTT100) <- rnam
colnames(sTHTT100) <- cnam
THTT100s <- xtable(sTHTT100, digit=4, align=algn)
print(THTT100s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100TH0d, echo=FALSE, results='asis'>>=
dTH0T100 <- dcors[[1]][[2]][,,1]
rownames(dTH0T100) <- rnam
colnames(dTH0T100) <- cnam
TH0T100d <- xtable(dTH0T100, digit=4, align=algn)
print(TH0T100d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100THTd, echo=FALSE, results='asis'>>=
dTHTT100 <- dcors[[1]][[2]][,,100+1]
rownames(dTHTT100) <- rnam
colnames(dTHTT100) <- cnam
THTT100d <- xtable(dTHTT100, digit=4, align=algn)
print(THTT100d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100TH0e, echo=FALSE, results='asis'>>=
eTH0T100 <- ecors[[1]][[2]][,,1]
rownames(eTH0T100) <- rnam
colnames(dTH0T100) <- cnam
TH0T100e <- xtable(eTH0T100, digit=4, align=algn)
print(TH0T100e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100THTe, echo=FALSE, results='asis'>>=
eTHTT100 <- ecors[[1]][[2]][,,100+1]
rownames(eTHTT100) <- rnam
colnames(eTHTT100) <- cnam
THTT100e <- xtable(eTHTT100, digit=4, align=algn)
print(THTT100e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=100$ for $\theta_0$ and $\theta_T$ for the state, scaled
  disturbance, and scaled error samplers. Row and column names indicate
  the true values of $W$ and $V$ respectively for the simulated
  data. Note that the signal-to-noise ratio is constant across any diagonal.}
\label{tableTHT100}
\end{table}


%Theta tables for T=1000
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000TH0s, echo=FALSE, results='asis'>>=
sTH0T1000 <- scors[[1]][[3]][,,1]
rownames(sTH0T1000) <- rnam
colnames(sTH0T1000) <- cnam
TH0T1000s <- xtable(sTH0T1000, digit=4, align=algn)
print(TH0T1000s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000THTs, echo=FALSE, results='asis'>>=
sTHTT1000 <- scors[[1]][[3]][,,1000+1]
rownames(sTHTT1000) <- rnam
colnames(sTHTT1000) <- cnam
THTT1000s <- xtable(sTHTT1000, digit=4, align=algn)
print(THTT1000s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000TH0d, echo=FALSE, results='asis'>>=
dTH0T1000 <- dcors[[1]][[3]][,,1]
rownames(dTH0T1000) <- rnam
colnames(dTH0T1000) <- cnam
TH0T1000d <- xtable(dTH0T1000, digit=4, align=algn)
print(TH0T1000d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000THTd, echo=FALSE, results='asis'>>=
dTHTT1000 <- dcors[[1]][[3]][,,1000+1]
rownames(dTHTT1000) <- rnam
colnames(dTHTT1000) <- cnam
THTT1000d <- xtable(dTHTT1000, digit=4, align=algn)
print(THTT1000d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000TH0e, echo=FALSE, results='asis'>>=
eTH0T1000 <- ecors[[1]][[3]][,,1]
rownames(eTH0T1000) <- rnam
colnames(dTH0T1000) <- cnam
TH0T1000e <- xtable(eTH0T1000, digit=4, align=algn)
print(TH0T1000e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000THTe, echo=FALSE, results='asis'>>=
eTHTT1000 <- ecors[[1]][[3]][,,1000+1]
rownames(eTHTT1000) <- rnam
colnames(eTHTT1000) <- cnam
THTT1000e <- xtable(eTHTT1000, digit=4, align=algn)
print(THTT1000e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=1000$ for $\theta_0$ and $\theta_T$ for the state, scaled
  disturbance, and scaled error samplers. Row and column names indicate
  the true values of $W$ and $V$ respectively for the simulated
  data. Note that the signal-to-noise ratio is constant across any diagonal.}
\label{tableTHT1000}
\end{table}


%% Theta Tables %%

%Theta tables for T=10
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT10TH0s, echo=FALSE, results='asis'>>=
sTH0T10 <- scors[[1]][[1]][,,1]
rownames(sTH0T10) <- rnam
colnames(sTH0T10) <- cnam
TH0T10s <- xtable(sTH0T10, digit=4, align=algn)
print(TH0T10s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10THTs, echo=FALSE, results='asis'>>=
sTHTT10 <- scors[[1]][[1]][,,10+1]
rownames(sTHTT10) <- rnam
colnames(sTHTT10) <- cnam
THTT10s <- xtable(sTHTT10, digit=4, align=algn)
print(THTT10s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10TH0d, echo=FALSE, results='asis'>>=
dTH0T10 <- dcors[[1]][[1]][,,1]
rownames(dTH0T10) <- rnam
colnames(dTH0T10) <- cnam
TH0T10d <- xtable(dTH0T10, digit=4, align=algn)
print(TH0T10d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10THTd, echo=FALSE, results='asis'>>=
dTHTT10 <- dcors[[1]][[1]][,,10+1]
rownames(dTHTT10) <- rnam
colnames(dTHTT10) <- cnam
THTT10d <- xtable(dTHTT10, digit=4, align=algn)
print(THTT10d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10TH0e, echo=FALSE, results='asis'>>=
eTH0T10 <- ecors[[1]][[1]][,,1]
rownames(eTH0T10) <- rnam
colnames(dTH0T10) <- cnam
TH0T10e <- xtable(eTH0T10, digit=4, align=algn)
print(TH0T10e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT10THTe, echo=FALSE, results='asis'>>=
eTHTT10 <- ecors[[1]][[1]][,,10+1]
rownames(eTHTT10) <- rnam
colnames(eTHTT10) <- cnam
THTT10e <- xtable(eTHTT10, digit=4, align=algn)
print(THTT10e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=10$ for $\theta_0$ and $\theta_T$ for the state, scaled
  disturbance, and scaled error samplers. Row and column names indicate
  the true values of $W$ and $V$ respectively for the simulated
  data. Note that the signal-to-noise ratio is constant across any diagonal.}
\label{tableTHT10}
\end{table}

%Theta tables for T=100
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT100TH0s, echo=FALSE, results='asis'>>=
sTH0T100 <- scors[[1]][[2]][,,1]
rownames(sTH0T100) <- rnam
colnames(sTH0T100) <- cnam
TH0T100s <- xtable(sTH0T100, digit=4, align=algn)
print(TH0T100s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100THTs, echo=FALSE, results='asis'>>=
sTHTT100 <- scors[[1]][[2]][,,100+1]
rownames(sTHTT100) <- rnam
colnames(sTHTT100) <- cnam
THTT100s <- xtable(sTHTT100, digit=4, align=algn)
print(THTT100s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100TH0d, echo=FALSE, results='asis'>>=
dTH0T100 <- dcors[[1]][[2]][,,1]
rownames(dTH0T100) <- rnam
colnames(dTH0T100) <- cnam
TH0T100d <- xtable(dTH0T100, digit=4, align=algn)
print(TH0T100d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100THTd, echo=FALSE, results='asis'>>=
dTHTT100 <- dcors[[1]][[2]][,,100+1]
rownames(dTHTT100) <- rnam
colnames(dTHTT100) <- cnam
THTT100d <- xtable(dTHTT100, digit=4, align=algn)
print(THTT100d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100TH0e, echo=FALSE, results='asis'>>=
eTH0T100 <- ecors[[1]][[2]][,,1]
rownames(eTH0T100) <- rnam
colnames(dTH0T100) <- cnam
TH0T100e <- xtable(eTH0T100, digit=4, align=algn)
print(TH0T100e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT100THTe, echo=FALSE, results='asis'>>=
eTHTT100 <- ecors[[1]][[2]][,,100+1]
rownames(eTHTT100) <- rnam
colnames(eTHTT100) <- cnam
THTT100e <- xtable(eTHTT100, digit=4, align=algn)
print(THTT100e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=100$ for $\theta_0$ and $\theta_T$ for the state, scaled
  disturbance, and scaled error samplers. Row and column names indicate
  the true values of $W$ and $V$ respectively for the simulated
  data. Note that the signal-to-noise ratio is constant across any diagonal.}
\label{tableTHT100}
\end{table}


%Theta tables for T=1000
\begin{table}[htb]
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000TH0s, echo=FALSE, results='asis'>>=
sTH0T1000 <- scors[[1]][[3]][,,1]
rownames(sTH0T1000) <- rnam
colnames(sTH0T1000) <- cnam
TH0T1000s <- xtable(sTH0T1000, digit=4, align=algn)
print(TH0T1000s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000THTs, echo=FALSE, results='asis'>>=
sTHTT1000 <- scors[[1]][[3]][,,1000+1]
rownames(sTHTT1000) <- rnam
colnames(sTHTT1000) <- cnam
THTT1000s <- xtable(sTHTT1000, digit=4, align=algn)
print(THTT1000s,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in State Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000TH0d, echo=FALSE, results='asis'>>=
dTH0T1000 <- dcors[[1]][[3]][,,1]
rownames(dTH0T1000) <- rnam
colnames(dTH0T1000) <- cnam
TH0T1000d <- xtable(dTH0T1000, digit=4, align=algn)
print(TH0T1000d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000THTd, echo=FALSE, results='asis'>>=
dTHTT1000 <- dcors[[1]][[3]][,,1000+1]
rownames(dTHTT1000) <- rnam
colnames(dTHTT1000) <- cnam
THTT1000d <- xtable(dTHTT1000, digit=4, align=algn)
print(THTT1000d,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Disturbance Sampler}
\vspace{.19 in}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000TH0e, echo=FALSE, results='asis'>>=
eTH0T1000 <- ecors[[1]][[3]][,,1]
rownames(eTH0T1000) <- rnam
colnames(dTH0T1000) <- cnam
TH0T1000e <- xtable(eTH0T1000, digit=4, align=algn)
print(TH0T1000e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_0$ in Scaled Error Sampler}
\end{minipage}
\begin{minipage}{.5\textwidth}
\centering
<<tableT1000THTe, echo=FALSE, results='asis'>>=
eTHTT1000 <- ecors[[1]][[3]][,,1000+1]
rownames(eTHTT1000) <- rnam
colnames(eTHTT1000) <- cnam
THTT1000e <- xtable(eTHTT1000, digit=4, align=algn)
print(THTT1000e,
  floating=FALSE,
  hline.after=NULL,
  scalebox=sclbx,
  add.to.row=list(pos=list(-1,0, length(rnam)),
  command=c('\\toprule\n','\\midrule\n','\\bottomrule\n')))
@
\captionof*{table}{$\theta_T$ in Scaled Error Sampler}
\end{minipage}
\caption{Autocorrelation in posterior sampler for a time series of
  length $T=1000$ for $\theta_0$ and $\theta_T$ for the state, scaled
  disturbance, and scaled error samplers. Row and column names indicate
  the true values of $W$ and $V$ respectively for the simulated
  data. Note that the signal-to-noise ratio is constant across any diagonal.}
\label{tableTHT1000}
\end{table}
