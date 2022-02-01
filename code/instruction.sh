# This is an instruction of how to set up my l501.F code on Klone

# In bash shell, load necessary modules
module load chem/pgi/18.10
module load chem/gdv/j14p

# Add the following lines in ~/.tcshrc
source /opt/ohpc/admin/lmod/lmod/init/csh
setenv LD_LIBRARY_PATH64 $gdvroot/gdv/bsd:$gdvroot/gdv/local:$gdvroot/gdv/extras:$gdvroot/gdv
source $gdvroot/gdv/bsd/gdv.login
setenv GAUSS_SCRDIR "/scr"

# Switch to a compute node with tcsh

# Make tree build
mktree -n gau_diss l502 l301
cd gau_diss/
cd nutil/
gau-get evlind utilam

# Then put the files into respective directories
scp l502_files.tar.gz al777@klone.hyak.uw.edu:/gscratch/chem/al777/gau_diss/l502
scp ginput.F al777@klone.hyak.uw.edu:/gscratch/chem/al777/gau_diss/l301

# Now the codes are in place, need to compile on a compute node
# Make sure pgi and gdv are loaded, then ask for a compute node
srun -p compute --job-name "InteractiveJob" --cpus-per-task 20 --mem=50G --time 4:00:00 --account chem  --pty /bin/tcsh

# Go to the top of the tree and compile
mk

# To run jobs, we can
# Option 1. At the top of the tree, run
gt FILENAME.gjf
# Option 2. From any directory, run
/sw/contrib/chem-src/gdv/j14p/gdv/bsd/gdvtest -w /gscratch/chem/al777/gau_diss -scrdir='/tmp' coh2_sto3g_protsp_new3.gjf





