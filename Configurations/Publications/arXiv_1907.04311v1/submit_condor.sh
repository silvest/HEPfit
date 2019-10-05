######################################################
# HTCondor Submit Description File.                  #
# COMMON TEMPLATE FOR FILE TRANSFER                  #
# Next commands that can be added to submit files    #
# -- Ayan Paul -- July 2018 -- DESY, Hamburg --      #
######################################################
# Picking up the NAME in a variable
FAME = NAME

# Picking up the ID
ID                      = $(Cluster).$(Process)

# Define universe and transfer protocols
# NOTE: Be careful with getenv
universe                = vanilla
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
getenv                  = True
# home_dir                = $ENV(HOME)   # If you want to point to your home

# Define directories and IO
# initialdir              = $(RUNDIR)
output                  = $(FAME).$(ID).out
error                   = $(FAME).$(ID).err
log                     = $(FAME).$(Cluster).log

# Set up notifications
notify_user             = ayan.paul@desy.de
notification            = Complete

# Hardware requests
request_cpus            = 16
request_memory          = 6 GB
requirements            = OpSysAndVer == "CentOS7" && Machine != "bird859.desy.de" && Machine != "batch1054.desy.de" && Machine != "batch1066.desy.de" && Machine != "batch1064.desy.de"
+RequestRuntime         = 324000

# Set up the run with exe and arguments. "config" is a directory.
# This is a MPI run within a node. Look into mpi-wrap.sh. 
transfer_input_files    = config
executable              = mpi-wrap.sh
arguments               = config/model.conf config/MonteCarlo.conf

# generate the queue from a list in another file
queue
# FOR MORE DETIALS:
# https://indico.cern.ch/event/611296/contributions/2604376/attachments/1471164/2276521/TannenbaumT_UserTutorial.pdf
