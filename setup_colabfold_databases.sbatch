#!/bin/bash
#SBATCH --job-name=colabfold_setup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=491G
#SBATCH --time=24:00:00
#SBATCH --output=colabfold_setup_%j.out
#SBATCH --error=colabfold_setup_%j.err

# Modified version of setup_databases.sh from ColabFold GitHub that uses a colabfold container
# This script downloads and decompresses the sequence databases, then uses mmseqs to create indexes for faster searching
# These indexes are compatible with both CPU and GPU sequence searches

#################
# CONFIGURATION #
#################

# Load any required modules (adjust as needed for your HPC)
module load apptainer

# Set working directory where databases will be stored (EDIT ME, or override: WORK_DIR=/path sbatch ...)
WORK_DIR="${WORK_DIR:-$PWD/colabfolddbs_gpu}"

# Path to the Apptainer container with ColabFold/MMseqs2, built from apptainer/colabfold.def
# (EDIT ME, or override: CONTAINER=/path/colabfold.sif sbatch ...)
CONTAINER="${CONTAINER:-$PWD/colabfold.sif}"

# Set number of threads to use (default: SBATCH --cpus-per-task)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

###############
# BEGIN SETUP #
###############
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# Copy the setup script into the working directory
cat > setup_databases.sh << 'SETUP_SCRIPT'
#!/bin/bash -e
# Setup everything for using mmseqs locally
mkdir -p tmp
cd tmp

if [ ! -f DOWNLOADS_READY ]; then
    aria2c --max-connection-per-server=8 --allow-overwrite=true "https://opendata.mmseqs.org/colabfold/uniref30_2302.db.tar.gz"
    aria2c --max-connection-per-server=8 --allow-overwrite=true "https://opendata.mmseqs.org/colabfold/colabfold_envdb_202108.db.tar.gz"
    aria2c --max-connection-per-server=8 --allow-overwrite=true "https://opendata.mmseqs.org/colabfold/uniref30_2302_newtaxonomy.tar.gz"
    touch DOWNLOADS_READY
fi

if [ ! -f ALL_FILES_EXTRACTED ]; then
  tar -xzvf "uniref30_2302.db.tar.gz"
  tar -xzvf "uniref30_2302_newtaxonomy.tar.gz"
  # contains replacement uniref30_2302_db_mapping, uniref30_2302_db_taxonomy, see https://github.com/sokrypton/ColabFold/wiki/MSA-Server-Database-History#2025-08-04-updated-uniref100_2302-taxonomypairing-files
  tar -xzvf "colabfold_envdb_202108.db.tar.gz"
  touch ALL_FILES_EXTRACTED
fi

# Make MMseqs2 merge the databases to avoid spamming the folder with files
export MMSEQS_FORCE_MERGE=1

if [ ! -f UNIREF30_READY ]; then
  mmseqs createindex "uniref30_2302_db" tmp1 --remove-tmp-files 1 ${INDEX_PAR} --split 1 --index-subset 2

  # Convert mapping to binary
  mmseqs createbintaxmapping "uniref30_2302_db_mapping" "uniref30_2302_db_mapping.bin"
  mv -f "uniref30_2302_db_mapping.bin" "uniref30_2302_db_mapping"
  # Create hardlinks (with new taxonomy)
  ln -f uniref30_2302_db_mapping uniref30_2302_db.idx_mapping
  ln -f uniref30_2302_db_taxonomy uniref30_2302_db.idx_taxonomy
  touch UNIREF30_READY
fi

if [ ! -f COLABDB_READY ]; then
    mmseqs createindex "colabfold_envdb_202108_db" tmp2 --remove-tmp-files 1 --split 1 --index-subset 2
    touch COLABDB_READY
fi

# Move index files and other necessary database files to main directory
cd ..
mv tmp/*.idx .
# These files are needed for ColabFold to work
cp tmp/uniref30_2302_db_mapping .
cp tmp/uniref30_2302_db_taxonomy .
cp tmp/colabfold_envdb_202108_db.dbtype .
cp tmp/uniref30_2302_db.dbtype .
# Other files in tmp are not needed after indexing
#rm -rf tmp

echo "Database setup complete! You can remove temporary files in ${WORK_DIR}/tmp if needed."
SETUP_SCRIPT

chmod +x setup_databases.sh

# Run the setup script inside the container
echo "Starting MMseqs2 database download and indexing..."
echo "Job started at: $(date)"
echo "Using ${SLURM_CPUS_PER_TASK} CPU threads"
echo "Allocated memory: ${SLURM_MEM_PER_NODE}M"

apptainer exec --bind ${WORK_DIR}:${WORK_DIR} \
    ${CONTAINER} \
    bash -c "cd ${WORK_DIR} && ./setup_databases.sh"

echo "Job completed at: $(date)"

# Print database statistics
echo ""
echo "=== Database files created ==="
ls -lh ${WORK_DIR}
