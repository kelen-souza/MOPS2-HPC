#!/bin/bash
#SBATCH --nodes=1            			#Número de Nós
#SBATCH --ntasks-per-node=1  			#Número de tarefas por Nó
#SBATCH --ntasks=1           			#Número de tarefas
#SBATCH -p sequana_cpu_dev   			#Fila (partition) a ser utilizada
#SBATCH --time=00:20:00			    	#Tempo máximo de execução
#SBATCH -J masa_test   		#Nome do job (job name) (EDITAR)
#SBATCH --exclusive        
#SBATCH -o slurm_%j.out                   #Output padrão do SLURM
#SBATCH -e slurm_%j.err                   #Output padrão do SLURM

#export OMP_NUM_THREADS=1
module load gcc/9.3_sequana
module load boost/1.87.0_gnu+openmpi-4.1.4_sequana
module load masa-openmp/1.0.1.1024_sequana

# SEQ_PATH="/scratch/cenapadrjsd/rafael.terra2/Workflow-MASAOpenMP-PAStar/src/saida/sequences"
SEQ_PATH="/scratch/cenapadrjsd/kelen.souza/workflow_sscad_2025/Workflow-MASAOpenMP-PAStar/resultados_genomas_completos_18seq_20h_1/sequences/"
SEQS=($(ls $SEQ_PATH)) # Utiliza as sequências separadas pelo workflow (X sequências, cada uma em um .fasta)
OUTDIR="/scratch/cenapadrjsd/rafael.terra2/Workflow-MASAOpenMP-PAStar/teste" # Salva em uma pasta qualquer
# Realiza a combinação par-a-par das sequencias
for ((i = 0; i < ${#SEQS[@]}; i++)); do 
    for ((j = i + 1; j < ${#SEQS[@]}; j++)); do
        # Seleciona dois arquivos do vetor de sequências 
        seq1="${SEQ_PATH}/${SEQS[i]}" 
        seq2="${SEQ_PATH}/${SEQS[j]}"
        # Verifica se os arquivos existem
        [ ! -f "${seq1}" ] && echo "Arquivo seq1 não encontrado: ${seq1}" && exit 1
        [ ! -f "${seq2}" ] && echo "Arquivo seq2 não encontrado: ${seq2}" && exit 1
        # Limpa a execução antiga e executa o MASA
        rm -rf ${OUTDIR}
        mkdir -p ${OUTDIR}
        time OMP_NUM_THREADS=1 masa-openmp --work-dir "${OUTDIR}" "${seq1}" "${seq2}"
        # time masa-openmp --work-dir "${OUTDIR}" "${seq1}" "${seq2}"
    done
done