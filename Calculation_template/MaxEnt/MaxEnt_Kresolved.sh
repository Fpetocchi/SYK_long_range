#!/bin/bash
#

################################################################################
#                              DEFAULT ARGUMENTS                               #
################################################################################
#
FIELD="G"
#
err="0.01"
mesh="3000"
width="200"
#
Nbatch="70"
Nktot="350"
#
DEFAULT=-10000
STARTK=$DEFAULT
STOPK=$DEFAULT
#
MODE="path"
PAD="K"
#
Nspin=1
SPIN_DW="F"
SINGSUB="F"
MODEL="F"
#
TRACE="F"
#
BINPATH=/home/petocchif/1_GW_EDMFT/nobackup/1_production/Calculation_template





################################################################################
#                             PROVIDED ARGUMENTS                               #
################################################################################
while getopts ":e:w:W:F:N:B:i:f:m:d:t:p:M:" o; do
   case ${o} in
      e)
         err="${OPTARG}"
         ;;
      w)
         mesh="${OPTARG}"
         ;;
      W)
         width="${OPTARG}"
         ;;
      F)
         FIELD="${OPTARG}"
         if [ "$FIELD"  != "G" ] && [ "$FIELD"  != "S" ]; then echo "Option Error - F" ; exit 1 ; fi
         ;;
      N)
         Nktot="${OPTARG}"
         ;;
      B)
         Nbatch="${OPTARG}"
         ;;
      i)
         STARTK="${OPTARG}"
         SINGSUB="T"
         ;;
      f)
         STOPK="${OPTARG}"
         SINGSUB="T"
         ;;
      m)
         MODE="${OPTARG}"
         if [ "$MODE"  != "path" ] && [ "$MODE"  != "full" ] && [ "$MODE"  != "plane" ]; then echo "Option Error - m" ; exit 1 ; fi
         if [ "$MODE"  == "path" ] ; then PAD="K" ; fi
         if [ "$MODE"  == "full" ] ; then PAD="F" ; fi
         if [ "$MODE"  == "plane" ] ; then PAD="P" ; fi
         ;;
      d)
         SPIN_DW="${OPTARG}"
         #if [ "$SPIN_DW"  != "T" ] && [ "$SPIN_DW"  != "F" ]; then echo "Option Error - s" ; exit 1 ; fi
         ;;
      t)
         TRACE="${OPTARG}"
         ;;
      p)
         BINPATH="${OPTARG}"
         ;;
      M)
         MODEL="${OPTARG}"
         if [ "$MODEL"  != "T" ] && [ "$MODEL"  != "F" ]; then echo "Option Error - M" ; exit 1 ; fi
         ;;
      \? )
         echo "Invalid option: $OPTARG" 1>&2
         exit 1
         ;;
      : )
         echo "Invalid option: $OPTARG requires an argument" 1>&2
         exit 1
         ;;
    esac
done
#
NAME=${FIELD}${PAD}
#
if [ "$SPIN_DW"  == "T" ]; then Nspin=2 ; fi
RUNOPTIONS=" -s "$err" -w "$mesh" -W "$width
#
TRPAD=""
if [ "${TRACE}"  == "T" ]; then TRPAD="_Tr" ; fi





################################################################################
#                                  PRINT INFOS                                 #
################################################################################
BIN=${BINPATH}/MaxEnt
echo
echo "Binary from: " ${BIN}
echo "Run options: "${RUNOPTIONS}
echo "Name: "${NAME}
echo "Nspin: "${Nspin}
echo "orbital Trace: "${TRACE}
if [ "$SINGSUB"  == "T" ]; then
   echo "startK: "${STARTK}
   echo "stopK: "${STOPK}
else
   echo "NKtot: "${Nktot}
   echo "Nbatch: "${Nbatch}
fi
echo




################################################################################
#                                  SUBMIT LOOP                                 #
################################################################################
for ispin in `seq 1 1 ${Nspin}` ; do

   #
   #Check if the I/O folder are present
   Gsource=MaxEnt_${FIELD}k_${MODE}_t_s${ispin}
   reports=${Gsource}_report
   if [ ! -d ${Gsource} ]; then
       echo ${Gsource} " (source) does not exists."
       exit 1
   fi
   if [ ! -d ${reports} ]; then mkdir ${reports} ; fi
   cd ${reports}

   #
   if [ ${STARTK} == ${DEFAULT} ] && [ ${STOPK} == ${DEFAULT} ]; then
      #
      #
      #Loop on batches
      kperbatch=`LC_ALL=en awk "BEGIN{print ${Nktot}/${Nbatch} }"`
      for btndx in `seq 1 1 ${Nbatch}` ; do
         #
         #Boundaries
         startk=`LC_ALL=en awk "BEGIN{print ${kperbatch}*(${btndx}-1) + 1 }"`
         stopk=` LC_ALL=en awk "BEGIN{print ${kperbatch}*${btndx}         }"`
         #
         #Info
         echo "batch: "${btndx}" startk: "${startk}" stopk: "${stopk}
         #
         #Edit one submit script for all the batches
         cat << EOF > submit_MaxEnt_${NAME}_B_${btndx}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N ${NAME}_s${ispin}_B${btndx}${TRPAD}
#$ -e error${NAME}_bt${btndx}.out
#$ -o log${NAME}_bt${btndx}.out
#$ -pe  smp 1
#$ -q   big.q,new.q

export PYTHONPATH=\${PYTHONPATH}:${BIN}/docopt/
export OMP_NUM_THREADS=1

for i in \`seq  ${startk} ${stopk}\`; do
   #
   echo K_\${i} > job_Kb_\${i}.out
   #
   if [ "$MODEL"  == "T" ] && [ -f ../${Gsource}/${FIELD}k_t_k\${i}.DAT_dos.dat ]; then
      #
      mpiexec -np 1 python3.6  ${BIN}/bryan.py ${RUNOPTIONS} -m  ../${Gsource}/${FIELD}k_t_k\${i}.DAT_dos.dat ../${Gsource}/${FIELD}k_t_k\${i}${TRPAD}.DAT >> job_Kb_\${i}.out
      #
   else
      #
      mpiexec -np 1 python3.6  ${BIN}/bryan.py ${RUNOPTIONS} ../${Gsource}/${FIELD}k_t_k\${i}${TRPAD}.DAT >> job_Kb_\${i}.out
      #
   fi
   #
   echo "MaxEnt on K_\${i} done" >> job_Kb_\${i}.out
   #
done
EOF
      #
      done
      #
      #submit all batches
      for btndx in `seq 1 1 ${Nbatch}` ; do
         qsub submit_MaxEnt_${NAME}_B_${btndx}
      done
      #
      #
   else
      #
      #
      #Loop on user-provided K-points
      for kp in `seq ${STARTK} ${STOPK}`; do
         #
         #Info
         echo "k point: "${kp}
         #
         #Edit one submit script for all the user-provided K-points
         cat << EOF > submit_MaxEnt_${NAME}_K_${kp}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N ${NAME}_s${ispin}_K${kp}${TRPAD}
#$ -e error${NAME}_K${kp}.out
#$ -o log${NAME}_K${kp}.out
#$ -pe  smp 1
#$ -q   big.q,new.q

export PYTHONPATH=\${PYTHONPATH}:${BIN}/docopt/
export OMP_NUM_THREADS=1

#
echo K_${kp} > job_K_${kp}.out
#
if [ "$MODEL"  == "T" ] && [ -f ../${Gsource}/${FIELD}k_t_k${kp}.DAT_dos.dat ]; then
   #
   mpiexec -np 1 python3.6  ${BIN}/bryan.py ${RUNOPTIONS} -m ../${Gsource}/${FIELD}k_t_k${kp}.DAT_dos.dat ../${Gsource}/${FIELD}k_t_k${kp}${TRPAD}.DAT >> job_K_${kp}.out
   #
else
   #
   mpiexec -np 1 python3.6  ${BIN}/bryan.py ${RUNOPTIONS} ../${Gsource}/${FIELD}k_t_k${kp}${TRPAD}.DAT >> job_K_${kp}.out
   #
fi
#
echo "MaxEnt on K_${kp} done" >> job_K_${kp}.out
#
EOF
      done
      #
      #submit all the user-provided K-points
      for kp in `seq ${STARTK} ${STOPK}`; do # -w
         qsub submit_MaxEnt_${NAME}_K_${kp}
      done
      #
      #
   fi

   #
   #Exiting from reports
   cd ../


done

exit 1
