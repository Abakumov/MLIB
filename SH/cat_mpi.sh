#!/bin/bash

DATANAME=${1:-"*"}
NUMDIP=${2:-19}
FIRSTINDX=${3:-0}

DATAFILES="${DATANAME}.angle \
           ${DATANAME}.RN \
           ${DATANAME}.Rnip \
           ${DATANAME}.singledipstack \
           ${DATANAME}.FuncEval \
           ${DATANAME}.iter \
           ${DATANAME}.coher \
           ${DATANAME}.Vnmo"

for DATAFILE in ${DATAFILES} ; do
#   for NR in {0..14} ; do 
    for ((NR=${FIRSTINDX};NR<$NUMDIP;NR++)); do
      rm -f ${DATAFILE}_${NR}.su
      FILES=$( ls -1 ${DATAFILE}_${NR}.* )
      FILESSORT=$( echo ${FILES[@]} | awk 'BEGIN{RS=" ";} {print $1}' | sort --version-sort -f )
      for FILE in ${FILESSORT} ; do
         cat $FILE >> ${DATAFILE}_${NR}_tmp.su
      done
	susort < ${DATAFILE}_${NR}_tmp.su > ${DATAFILE}_${NR}.su
	rm -f ${DATAFILE}_${NR}_tmp.su
   done
done

DATAFILES="${DATANAME}.angle \
           ${DATANAME}.RN \
           ${DATANAME}.Rnip \
           ${DATANAME}.stack \
           ${DATANAME}.Fstack \
           ${DATANAME}.FuncEval \
           ${DATANAME}.iter \
           ${DATANAME}.maxind \
           ${DATANAME}.coher \
           ${DATANAME}.Vnmo \
           ${DATANAME}.DiffStack \
           ${DATANAME}.DiffbinStack \
           ${DATANAME}.coh_weight_stack"

for DATAFILE in ${DATAFILES} ; do
   rm -f ${DATAFILE}.su
   FILES=$( ls ${DATAFILE}[0-9]* )
   FILESSORT=$( echo ${FILES[@]} | awk 'BEGIN{RS=" ";} {print $1}' | sort --version-sort -f )
   for FILE in ${FILESSORT} ; do
      cat $FILE >> ${DATAFILE}_tmp.su
   done
	susort < ${DATAFILE}_tmp.su > ${DATAFILE}.su
	rm -f ${DATAFILE}_tmp.su
done
