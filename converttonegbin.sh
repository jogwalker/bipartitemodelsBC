#!/bin/bash
#
#PBS -N jgwZINBconvert
#PBS -o output
#PBS -e error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash


filepath="$HOME/bipartitemodelsBC/finalmodels/jags/phase4"
outdir="$HOME/bipartitemodelsBC/finalmodels/jagsNB/phase4"

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})


for j in "${array[@]}"
do
	sed -e 's/count\[i\] ~ dpois(mu\[i\])/count\[i\] ~ dnegbin(Pp\[i\], r)/' -e 's_mu\[i\] <- lambda\[i\] \* use\[host\.sp\[i\], par\.sp\[i\]\]_Pp\[i\] <- r/(r+mu.eff\[i\]) \n mu.eff\[i\] <- lambda\[i\]* use\[host.sp\[i\], par.sp\[i\]\]_' -e 's_\(for (k in 1\:[0-9])\)_r~dgamma(0.1,0.1) \n \1_' < ${filepath}/${j} > ${outdir}/${j}
done

