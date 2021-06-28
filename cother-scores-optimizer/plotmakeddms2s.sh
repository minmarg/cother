#!/bin/bash

basename=$( basename $0 )
dirname=$( dirname $0 )
[[ "${dirname:0:1}" != "/" ]] && dirname="$( pwd )/$dirname"

usage="
Plot ddms2s scores using R.
2021(C)Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$basename <Parameters>

Parameters:

-i <input>     Input file of ddms2s scores (produced by makeddms2s++.pl).
-o <output>    Output pdf filename.
-a <slopes>    Comma-separated list of slopes for additional lines.
-b <intercepts> Comma-separated list of intercepts for additional lines.
-h             This text.
"


while getopts "i:o:a:b:h" Option
do
    case $Option in
        i ) INPUT=${OPTARG} ;;
        o ) OUTPUT=${OPTARG} ;;
        a ) SLOPES=${OPTARG} ;;
        b ) INTERS=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo Error: Unrecognized argument.; exit 1 ;;
    esac
done
shift $(( $OPTIND - 1 ))

if [[ -z "${OUTPUT}" ]]; then echo "ERROR: No output filename given."; exit 1; fi
if [[ -z "${INPUT}" ]]; then echo "ERROR: No input filename given."; exit 1; fi
if [[ ! -f "${INPUT}" ]]; then echo "ERROR: Input file not found: ${INPUT}"; exit 1; fi


Rassign="$(perl -e '
@anms=("s","str","fp","fn");
while(<>){
  next unless /^(?:#+)?\s+([\-\d\.]+)\s+([\-\d\.]+)/;
  s/^#+(.+)$/$1/;
  last if $c++ > $#anms;
  push @{$a[$c-1]}, $1 while(/\s+(\S+)/g);}
exit 1 if $#a < $#anms;
for($i=0; $i<=$#{$a[0]}; $i++) {
  $b=0;
  for($j=0;$j<=$#a;$j++){if($a[$j][$i]=~/[^\-\d\.]/){$b=1;last}}
  next if $b;
  for($j=0;$j<=$#a;$j++){@out[$j].=",$a[$j][$i]"} }
for($j=0; $j<=$#out; $j++){
  substr($out[$j],0,1)=" ";
  $out[$j] = "$anms[$j]<-c(". $out[$j] . ")\n";
  print $out[$j]}' ${INPUT})"


Rtext="

# ib<-2;
# ie<-56;

${Rassign};

# s<-s[ib:ie];
# str<-str[ib:ie];
# fp<-fp[ib:ie];
# fn<-fn[ib:ie];

sf<-s[(fp>0.001) & (fn>0.001)];
strf<-str[(fp>0.001) & (fn>0.001)];

pcc<-cor(sf,strf,method='pearson');

print(s);
print(sum(fp));
print(sum(fn));
print(pcc);

fit<-lm(strf~sf);
a<-fit\$coefficients['(Intercept)'];
b<-fit\$coefficients['sf'];

print(fit);
print(anova(fit));

pdf(\"${OUTPUT}\", width=3.0, height=3.0, pointsize=1);

par(mar=c(4.9, 4.9, 1.2, 1))
par(cex=1,cex.main=1.5,cex.axis=1.8)

plot(s,str,type='l',log='',
    main='',
    axes=F,xlab='',ylab='',lwd=0.2,col='black');

lines(sf,strf,log='',lwd=2,col='black');

abline(a,b,lty=5,cex=1,col='grey50');

aaux<-c(${INTERS});
baux<-c(${SLOPES});

naux<-length(aaux)
nauxcols<-rainbow(naux);

if(length(aaux) && length(aaux)==length(baux)) {
    sapply(1:naux, 
      function(j){
        abline(aaux[j],baux[j],lty=5,cex=1,col=nauxcols[j])
    });
}

sat0<-s[which.min(abs(s+a/b))][1]
abline(h=0,lwd=0.2,col='grey80');
abline(v=sat0,lwd=0.2,col='grey80');
text(sat0,min(str),substitute(v,list(v=round(sat0,2))),adj=c(1,0),cex=1.5);

par(new=T)

axis(2, col='black');
mtext(expression(paste(italic('S')['tr'])),side=2,las=2,adj=1.2,line=2.5,cex=1.8);

box();
axis(1,col='black');
mtext(expression(paste(italic('S'))),side=1,padj=0.5,col='black',line=2.5,cex=1.8);

cds<-par('usr');
dx<-cds[2]-cds[1];
dy<-cds[4]-cds[3];
text(cds[1]+dx*.05,(cds[4]-dy*.05),substitute(italic(r)==c,list(c=round(pcc,2))),adj=c(0,1),cex=1.9);

if(naux<1) {
    legend('bottomright',lty=5,bty='n',cex=1.5,col=c('grey50'),
        legend=substitute(b*italic(S)~s~a,list(a=round(abs(a),3),b=round(b,3),s=ifelse(a<0,'-','+'))));
} else {
    legend('bottomright',lty=5,bty='n',cex=1.5,col=c('grey50',nauxcols),
        legend=c(
          as.expression(
            substitute(b*italic(S)~s~a,list(a=round(a,3),b=round(b,3),s=ifelse(a<0,'','+')))),
          sapply(1:naux,function(j){
            as.expression(
              substitute(b*italic(S)~s~a,
                list(a=round(abs(aaux[j]),3),b=round(baux[j],3),s=ifelse(aaux[j]<0,'-','+'))))
          })
    ));
}

dev.off();
"
( echo "${Rtext}" | R --vanilla )

