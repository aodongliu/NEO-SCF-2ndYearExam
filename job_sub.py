#!/usr/bin/python

from __future__ import print_function

import sys
import os
import re

import textwrap
try:
    from subprocess import check_output
    from subprocess import Popen, PIPE
except ImportError:
    print('\nERROR: Must load at least anaconda_2.3 environment on Ikt\n')
    sys.exit()

'''
  Patrick J. Lestrange 2017

  gaussian-sub.py: Builds a PBS submission script to run
                   Gaussian on the Hyak Ikt/Mox clusters.

                   Also used as a teaching aid in UW's
                   CHEM 465/565 Computations in Chemistry.

  Maintained by Andrew Wildman
'''

#----------------------------------------------------------------------------
def get_user_input():
    """Grab input from the user about what type of job to run."""

    global f_input, gdv, queue, allocation, version, n_nodes
    global linda, n_cores, time, f_output, gen, memory
    global lclScr, email, memory_write
    gdv = False

    #--------------------------------------
    # Determine which generation machine we're on
    gen  = 'ikt'
    host = Popen('hostname',stdout=PIPE).stdout.read()
    if 'mox' in host: gen = 'mox'
    #--------------------------------------

    #--------------------------------------
    # Check arguments
    if len(sys.argv) == 1:
        print('No files specified')
        exit()
    if len(sys.argv) != 2:
        print("Using the same parameters for all files")
        print(sys.argv[1:])
        print()
    if sys.argv[1] == '-h' or 'help' in sys.argv[1]:
        print_help()
        sys.exit()
    #--------------------------------------

    #--------------------------------------
    # Check that the argument is a .com or .gjf and that it exists.
    f_input = []
    for fil in sys.argv[1:]:
        f_temp= fil.split('.')
        if not os.path.isfile(fil):
            print('ERROR: The input file '+fil+' does not exist')
            sys.exit()
        if len(f_temp) != 2:
            print('ERROR: The filename must include the extension')
            sys.exit()
        if f_temp[-1] != 'com' and f_temp[-1] != 'gjf':
            print('ERROR: The file extension must be .com or .gjf')
            sys.exit()
        f_input.append(f_temp)
    #--------------------------------------

    #--------------------------------------
    # Check that the user has the right permissions to use Gaussian.
    username = Popen('whoami',stdout=PIPE).stdout.read().strip()
    command = 'groups '+username
    groups = Popen(command,stdout=PIPE,shell=True).stdout.read().split(' ')
    groups[-1] = groups[-1].strip()
    if 'ligroup-gaussian' not in groups:
        print(textwrap.fill(textwrap.dedent("""\
            ERROR: You must be part of the ligroup-gaussian Unix group
            to use Gaussian. Contact Prof. Xiaosong Li (xsli@uw.edu) to
            be added to the group."""), 60))
        sys.exit()
    if 'ligroup-gdv' in groups: gdv = True
    #--------------------------------------

    #--------------------------------------
    # Determine which queue to submit to.
    queue = raw_input(textwrap.fill(textwrap.dedent("""\
            Which queue would you like to submit to?
            [batch] - (default) or [ckpt] : """),100))
    if queue == '': queue = 'batch'
    if queue == 'bf': queue = 'ckpt'
    if queue != 'batch' and queue != 'ckpt':
        print('ERROR: Invalid option for a queue. Must be batch or ckpt')
        sys.exit()
    print('Using the '+queue+' queue\n')
    #--------------------------------------

    #--------------------------------------
    # Determine which group's nodes to use.
    allocation = ''
    if queue == 'batch' or queue == 'ckpt':
        allocs= []
        command = 'sinfo -o "%P"'
        partitions = Popen(command,stdout=PIPE,shell=True).stdout.read().split()[1:]
        for group in groups:
            if 'hyak-' in group:
                if 'test' not in group and 'highmem' not in group:
                    if group.split('-')[1] in partitions or group.split('-')[1] == 'genpool':
                        allocs.append(group)
        print('Whose allocation would you like to use?')
        for allocation in allocs:
            print('['+allocation+']',end=' ')
        print(': ',end='')
        allocation = raw_input('')
        if allocation == '':
            print(textwrap.fill(textwrap.dedent("""\
                ERROR: You must specify an allocation"""),100))
            sys.exit()
            #allocation = 'hyak-stf'
        elif allocation not in allocs:
            print(textwrap.fill(textwrap.dedent("""\
                ERROR: You must choose an allocation that you
                are a part of"""),100))
            sys.exit()
        print('Submitting to the '+allocation+' allocation\n')
    #--------------------------------------

    #--------------------------------------
    # Ask how many nodes to use.
    if queue == 'batch':
        print('Checking how many nodes are in this allocation...')
        allocation_name = allocation.split('-')
        if allocation_name[1] == 'genpool':
            allocation_name[1] = 'hpc'
        command = 'hyakalloc | grep '+allocation_name[1]
        print(command)
        specs = Popen(command,stdout=PIPE,shell=True).stdout.read().split()
        max_nodes = int(specs[1])
        smallest_mem = int(specs[2][:-1])
    else:
        max_nodes = 1000

    n_nodes = raw_input('How many nodes do you want to use? (default=1) : ')
    if n_nodes == '':
        n_nodes = 1
    else:
        n_nodes = int(n_nodes)

    if max_nodes == 0:
        print(textwrap.fill(textwrap.dedent("""\
            ERROR: There are no nodes available for this allocation"""),100))
        sys.exit()
    elif n_nodes < 1 or n_nodes >= max_nodes:
        print(textwrap.fill(textwrap.dedent("""\
            ERROR: You must select at least one node and
            less than %d""" % max_nodes),100))
        sys.exit()

    if n_nodes > 1:
        linda = True
    else:
        linda = False

    if linda:
        print(textwrap.fill(textwrap.dedent("""\
            NOTE: You must include "%lindaworker" in your input file when
            using more than one node. Some version also require including
            %UseSSH in your Gaussian input file."""),60))
    #--------------------------------------
    # Ask to set the local scratch during printing
    lclScr = raw_input('Set scratch locally? (y or n): ')
    writeScr = None
    if lclScr == '' or lclScr == 'n': lclScr = 0
    else: lclScr = 1

    #--------------------------------------
    #Ask about where to send email notifications
    whoami = re.sub('\n', '', os.popen('whoami').readlines()[0])
    email = raw_input("Is %s@uw.edu the correct email, else what is?: " %whoami)
    if email == '': email = whoami

    #--------------------------------------
    # Ask how many cores/memory on each node to use.
    print('Checking what types of nodes are in this allocation...')
    if queue == 'batch':
        command = 'sinfo -p ' + allocation_name[1] + ' -e -O "nodes,cpus"'
        specs = Popen(command, stdout=PIPE, shell=True).stdout.read().split()
    else:
        command = 'sinfo -p ckpt -e -O "nodes,cpus"'
        specs = Popen(command, stdout=PIPE, shell=True).stdout.read().split()
    nodes_per_cpu = {}
    for i in range(1, len(specs)//2):
        nodes_per_cpu[int(specs[2*i+1])] =  int(specs[2*i])
    smallest_node = min(nodes_per_cpu)
    n_cores = raw_input(textwrap.fill(textwrap.dedent("""\
              How many cores do you want to use
              on each node? (default=%d) : """ % smallest_node).strip()))

    if n_cores == '':
        n_cores = 0
    else:
        n_cores = int(n_cores)

    if n_cores > max(nodes_per_cpu):
        print(textwrap.fill(textwrap.dedent("""\
            You cannot request more than the maximum
            number of cores on this allocation. (%d)""" % (max(nodes_per_cpu))),100))
        exit()

    try:
        max_nodes = nodes_per_cpu[n_cores]
    except KeyError:
        for cores in sorted(nodes_per_cpu):
            if n_cores < cores:
                max_nodes = nodes_per_cpu[cores]
                print(textwrap.fill(textwrap.dedent("""\
                    Requesting all the cores on the smallest node (%d)
                    that fits the requested number (%d)""" % (cores, n_cores)),100))
                n_cores = cores
                break

    if n_nodes > max_nodes:
        n_nodes = max_nodes
        print('You requested too many nodes with '+n_cores+' cores')
        print('Resetting to maximum number of nodes with that many cores (%d)' %(max_nodes))


    print('Checking available memory for the requested nodes...')
    if queue == 'batch':
        command = 'sinfo -p ' + allocation_name[1] + ' -e -O "cpus,memory" | grep ^' + str(n_cores)
        specs = Popen(command, stdout=PIPE, shell=True).stdout.read().split()
    else:
        command = 'sinfo -p ckpt -e -O "cpus,memory" | grep ^' + str(n_cores)
        specs = Popen(command, stdout=PIPE, shell=True).stdout.read().split()
    # - 10 for os overhead
    mem_types = [int(specs[2*i+1])//1000 - 10 for i in range(len(specs)//2)]
    #smallest_mem = min(mem_types)
    smallest_mem = 0 #Default all mem
    max_mem = max(mem_types)

    memory = raw_input(textwrap.fill(textwrap.dedent("""\
              How much memory, in Gb, do you want to use
              on each node? (default=All Available) : """).strip()))
    if memory == '':
        memory = smallest_mem
        memory_write='0'
    else:
        memory = int(memory)
        memory_write=str(memory)+'G'
        if memory > max_mem:
            print(textwrap.fill(textwrap.dedent("""\
                You cannot request more than the maximum
                memory on this type of node. (%d)""" % (max_mem),100)))
            exit()

    print('Using %d node(s) with %d cores and %d Gb\n' %
        (n_nodes, n_cores, memory))
    #--------------------------------------

    #--------------------------------------
    # Determine which version of Gaussian to use.
    if gen == 'ikt':
        g09_versions = ['d01','e01']
        g16_versions = ['a03']
        gdv_versions = ['i03','i03p','i04p','i06','i06p','i09']
    elif gen == 'mox':
        g09_versions = ['a02','e01']
        g16_versions = ['a03', 'b01']
        gdv_versions = ['i03','i03p','i04p','i06p','i10pp','i11p','i13','i13p']
    print('Which version of Gaussian would you like to use?')
    for version in g09_versions: print('[g09.'+version+']',end=' ')
    for version in g16_versions: print('[g16.'+version+']',end=' ')
    if gdv:
        for version in gdv_versions: print('[gdv.'+version+']',end=' ')
    print('- (default) : ',end='')
    version = raw_input('')
    if version == '':
        if gdv: version = 'gdv.'+gdv_versions[-1]
        else: version = 'g16.'+g16_versions[-1]
    version_name = version.split('.')
    if len(version_name) != 2:
        print("Can't recognize the specificed version")
        sys.exit()
    bad_version = False
    if version_name[0] == 'g09':
        if version_name[1] not in g09_versions:
            bad_version = True
    elif version_name[0] == 'g16':
        if version_name[1] not in g16_versions:
            bad_version = True
    elif version_name[0] == 'gdv':
        if version_name[1] not in gdv_versions:
            bad_version = True
    else:
        bad_version = True
    if bad_version:
        print('Choose a valid version of Gaussian. Not %s.%s'
            % (version_name[0], version_name[1]))
        sys.exit()
    print('Using the '+version+' version of Gaussian\n')
    #--------------------------------------

    #--------------------------------------
    # Check what the max walltime should be.
    max_time = 750000
    lim = 240
    if queue != 'bf' and queue != 'ckpt':
        default = 1
        unit    = 'hr'
        time = raw_input(textwrap.fill(textwrap.dedent("""\
               For how many hours do you want to run your
               calculation? (default=%d hr) : """ % default).strip(),100))
    else:
        default = 6
        unit    = 'hr'
        time = raw_input(textwrap.fill(textwrap.dedent("""\
               For how many hours do you want to run your
               calculation? (default=%d hr) : """ % default).strip(),100))
    if time == '': time = default
    else: time = int(time)
    if queue == 'bf' and time < 6:
        print(textwrap.fill(textwrap.dedent("""\
            If you want your job to be resubmitted automatically on
            the ckpt partition, you need to specify greater than
            5 hours of runtime. This is just a warning."""),100))
    print('Running the calculation for %d %s(s)\n' % (time, unit))

    #Check on MAXTIME
    time_lim = max_time / 60 / n_nodes / n_cores
    if time > time_lim  or time > lim:#Checking doesn't exceed max time or 10 d
      print(textwrap.fill(textwrap.dedent("""\
            WARNING:
            You have asked for %i hours on %i nodes with %i processing cores.
            If you are submitting to the stf queue this will exceed your
            timelimit (%d hours).  To correct this either request less time,
            fewer nodes, or fewer processing cores. If this causes an issue,
            or you have questions please contact hpcc@uw.edu.\n""" % \
            (time, n_nodes, n_cores, min(time_lim, lim))),100))


    #--------------------------------------

    #--------------------------------------
    # Get a name for the .pbs script.
    extension = 'sh'
    if len(f_input) == 1:
        f_output = raw_input(textwrap.fill(textwrap.dedent("""\
                   What should the .%s script be named?
                   (default=%s.%s) : """ % ( extension, f_input[0][0].strip(),
                                             extension) )))
        if f_output == '': f_output = [f_input[0][0]+'.'+extension]
        else: f_output = [f_output+'.'+extension]
    # AW - assume default if multiple files submitted
    else:
        f_output = [f_i[0]+'.'+extension for f_i in f_input]
    #--------------------------------------

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
def check_Gaussian_input():
    """Read the Gaussian input file and check for problems."""

    for fil in f_input:
        gauss_input = str(fil[0])+'.'+str(fil[1])
        f = open(gauss_input,'r')
        contents = f.readlines()
        found_linda, found_ssh, exit, gb = False, False, False, False
        memory, nproc = 0, 1
        warnings = []

        for line in contents:
            if 'lindaworker' in line.lower(): found_linda = True
            if 'usessh' in line.lower(): found_ssh = True
            if 'mem' in line.lower():
                mem_line = line.split('=')
                mem_line[-1] = mem_line[-1].strip()
                if 'gb' in mem_line[1].lower():
                    gb = True
                    mem = mem_line[1].lower().split('gb')
                    memory = int(mem[0])
                else:
                    gb = False
                    warning = textwrap.dedent("""\
                        This script only checks the memory specfication if it
                        is in Gb. Your calculation may still be fine, but
                        this script won't check. This is just a warning.""")
                    warnings.append(warning)
            if 'nproc' in line.lower():
                nproc_line = line.split('=')
                nproc_line[-1] = nproc_line[-1].strip()
                nproc = int(nproc_line[1])
            if 'chk' in line.lower():
                pwd = Popen('pwd',stdout=PIPE,shell=True).stdout.read().strip()
                if 'c:' in line.lower():
                    warning = textwrap.dedent("""\
                        Your checkpoint file includes the C: drive and there is
                        no C: drive on this machine. Please update the path for
                        your checkpoint file.""")
                    warnings.append(warning)
                    exit = True
                if '/' in line.lower() or '\\' in line.lower() and \
                    pwd not in line.lower():
                    warning = textwrap.dedent("""\
                        The path specified for your checkpoint file is not the
                        current directory. You may want to change this. This is
                        just a warning.""")
                    warnings.append(warning)

        if gb and allocation == 'hyak-stf':
            if gen == 'ikt' and memory > 32:
                warning = textwrap.dedent("""\
                    Generally you don't want to specify more than half the
                    memory on a node. You've asked for %dGb and most STF
                    nodes only have 64Gb.
                    This is just a warning.""" % memory)
                warnings.append(warning)
            elif gen == 'mox' and memory > 64:
                warning = textwrap.dedent("""\
                    Generally you don't want to specify more than half the
                    memory on a node. You've asked for %dGb and most STF
                    nodes only have 128Gb.
                    This is just a warning.""" % memory)
                warnings.append(warning)

        if allocation == 'hyak-stf':
            if gen == 'ikt' and nproc > 16:
                warning = textwrap.dedent("""\
                    You should not specify to use more cores than the number
                    available on your node. The STF nodes have 16 cores
                    and you've asked for %d cores. Please lower the number
                    of cores you've requested in your input file.
                    Not forming PBS script.""" % nproc)
                warnings.append(warning)
                exit = True
            elif gen == 'mox' and nproc > 28:
                warning = textwrap.dedent("""\
                    You should not specify to use more cores than the number
                    available on your node. The STF nodes have 28 cores
                    and you've asked for %d cores. Please lower the number
                    of cores you've requested in your input file.
                    Not forming SBATCH script.""" % nproc)
                warnings.append(warning)
                exit = True

        if allocation == 'hyak-stf':
            if gen == 'ikt' and nproc < 16:
                warning = textwrap.dedent("""\
                    You usually want to use all the cores on a node. The
                    STF nodes have 16 cores and you've asked for %d core(s).
                    This is just a warning.""" % nproc)
                warnings.append(warning)
            elif gen == 'mox' and nproc < 28:
                warning = textwrap.dedent("""\
                    You usually want to use all the cores on a node. The
                    STF nodes have 28 cores and you've asked for %d cores.
                    This is just a warning.""" % nproc)
                warnings.append(warning)

        if linda and not found_linda:
            warning = textwrap.dedent("""\
                Your input file does not contain %lindaworker, but
                you have asked to use more than one node. Please add this
                line or request only one node. Not forming PBS script.""")
            warnings.append(warning)
            exit = True

        if linda and not found_ssh:
            warning = textwrap.dedent("""\
                Your input file does not contain %UseSSH, but
                you have asked to use more than one node. This is can be
                a problem for some versions of Gaussian. Please add this
                line or request only one node. This is just a warning.""")
            warnings.append(warning)

        if found_linda and n_nodes == 1:
            warning = textwrap.dedent("""\
                Your input file contains %lindaworker, but you have
                only asked to use one node. Please remove this line or
                request more than one node. Not forming PBS script.""")
            warnings.append(warning)
            exit = True

        # Print warnings and exit if there are too many errors.
        if len(warnings) > 0:
          print('\n'+'#'*40+'\n'
               +' '*15+'WARNINGS\n'
               +'#'*40)
        for warning in warnings:
            print('\n'+textwrap.fill(warning, 60))
        if exit:
          print("\nExiting without writing PBS file.\n")
          sys.exit()
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
def write_torque_script():
    """Make a .pbs script based on user specifications."""

    for i in range(len(f_input)):

        gauss_input = str(f_input[i][0])+'.'+str(f_input[i][1])
        print('Writing to '+f_output[i]+'\n')
        pwd = Popen('pwd',stdout=PIPE,shell=True).stdout.read().strip()
        f = open(f_output[i],'w')

        f.write(textwrap.dedent("""\
            #!/bin/bash
            #PBS -N %s
            #PBS -l nodes=%d:ppn=%d,feature=%dcore"""
            % (f_input[i][0], n_nodes, n_cores, n_cores)))
        if queue != 'bf':
            f.write('\n#PBS -l walltime=%d:00:00\n' % time)
        else:
            f.write('\n#PBS -l walltime=0:%d:00\n' % time)
        f.write(textwrap.dedent("""\
            #PBS -j oe
            #PBS -o %s
            #PBS -d %s\n""" % (pwd, pwd)))
        if queue == 'batch':
            f.write('#PBS -W group_list=%s\n' % allocation)
        f.write(textwrap.dedent("""\
            #PBS -q %s

            # load Gaussian environment
            module load contrib/%s

            # debugging information
            HYAK_NPE=$(wc -l < $PBS_NODEFILE)
            HYAK_NNODES=$(uniq $PBS_NODEFILE | wc -l )
            echo "**** Job Debugging Information ****"
            echo "This job will run on $HYAK_NPE CPUs on $HYAK_NNODES nodes"
            echo ""
            echo Node:CPUs Used
            uniq -c $PBS_NODEFILE | awk '{print $2 ":" $1}'
            echo "SHARED LIBRARY CHECK"
            ldd ./test
            echo "ENVIRONMENT VARIABLES"
            set
            echo "**********************************************" """
            % (queue, version)))

        if linda:
            f.write(textwrap.dedent("""\
                \n
                # add linda nodes
                HYAK_NNODES=$(uniq $PBS_NODEFILE | wc -l )
                nodes=()
                nodes+=(`uniq -c $PBS_NODEFILE | awk '{print $2}'`)
                for ((i=0; i<${#nodes[*]}-1; i++));
                do
                \tstring+=${nodes[$i]}
                \tstring+=","
                done
                string+=${nodes[$HYAK_NNODES-1]}
                sed -i -e "s/%%LindaWorker.*/%%LindaWorker=$string/Ig" %s

                # check that the Linda nodes are correct
                lindaline=(`grep -i 'lindaworker' %s`)
                if [[ $lindaline == *$string ]]
                then
                \techo "Using the correct nodes for Linda"
                else
                \techo "Using the wrong nodes for Linda"
                \techo "Nodes assigned by scheduler = $string"
                \techo "Line in Gaussian input file = $lindaline"
                \texit 1
                fi """ % (gauss_input, gauss_input)))

        if queue == 'bf':
            f.write(textwrap.dedent("""\
              \n
              # copy last log file to another name
              num=`ls -l %s*.log | wc -l`
              cp %s.log %s$num.log""" % (f_input[i][0], f_input[i][0],
              f_input[i][0])))

        if 'gdv' in version:
            command = 'gdv'
        elif 'g16' in version:
            command = 'g16'
        else:
            command = 'g09'
        f.write(textwrap.dedent("""\
            \n
            # run Gaussian
            %s %s

            exit 0 """ % (command, gauss_input)))

        print("""Please run 'qsub %s' to submit to the scheduler\n"""
              % f_output[i])
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
def write_slurm_script():
    """Make a .sh script based on user specifications."""

    # Loop over the input files
    for i in range(len(f_input)):

        gauss_input = str(f_input[i][0])+'.'+str(f_input[i][1])
        print('Writing to '+f_output[i]+'\n')
        pwd = Popen('pwd',stdout=PIPE,shell=True).stdout.read().strip()
        f = open(f_output[i],'w')
        short_name = re.split('-',allocation)[1]
        partition, account = short_name, short_name
        if account == 'genpool':
            partition = 'hpc'
        if queue == 'bf' or queue == 'ckpt':
            partition = queue
            account = short_name+'-ckpt'

        f.write(textwrap.dedent("""\
            #!/bin/bash
            #SBATCH --job-name=%s
            #SBATCH --nodes=%d
            #SBATCH --cpus-per-task=%d""" % (f_input[i][0], n_nodes, n_cores)))
        f.write('\n#SBATCH --time=%d:00:00\n' % time)
        f.write(textwrap.dedent("""\
            #SBATCH --mem=%s
            #SBATCH --chdir=%s
            #SBATCH --mail-type=FAIL,END
            #SBATCH --mail-user=%s@uw.edu

            #SBATCH --partition=%s
            #SBATCH --account=%s\n\n"""
            % (memory_write, pwd, email, partition, account)))
        f.write(textwrap.dedent("""\
            # load Gaussian environment
            module load contrib/%s
            export inputfile='%s'


            # debugging information
            echo "**** Job Debugging Information ****"
            echo "This job will run on $SLURM_JOB_NODELIST"
            echo ""
            echo "ENVIRONMENT VARIABLES"
            set
            echo "**********************************************" """
            % (version,gauss_input)))

        if lclScr != 0:
            f.write(textwrap.dedent("""\
                \n
                # local scratch
                export GAUSS_SCRDIR='%s'
                """ % (pwd)))

        else:
            f.write(textwrap.dedent("""\
                \n
                # scrubbed scratch
                export scrDir='/gscratch/scrubbed/%s/'
                mkdir -p $scrDir
                export GAUSS_SCRDIR=$scrDir
                """ % (email)))

        f.write(textwrap.dedent("""\
            \n
            ## Memory
            gbmem=`expr $SLURM_MEM_PER_NODE / 1000`
            gbmem=`expr $gbmem - 10`
            echo "Parsed memory: $gbmem"
            sed -i "/mem/s/.*/%mem=${gbmem}GB/" $inputfile
            """))

        f.write(textwrap.dedent("""\
            \n
            ## Set number of threads
            export num_threads=$(echo $SLURM_JOB_CPUS_PER_NODE| cut -f1 -d"(" )
            sed -i "/nproc/s/.*/%nprocshared=${num_threads}/" $inputfile
            """))


        if linda:
            f.write(textwrap.dedent("""\
                \n
                # add linda nodes
                nodes=()
                nodes+=(`scontrol show hostnames $SLURM_JOB_NODELIST `)
                for ((i=0; i<${#nodes[*]}-1; i++));
                do
                \tstring+=${nodes[$i]}
                \tstring+=","
                done
                string+=${nodes[$SLURM_NNODES-1]}
                sed -i -e "s/\%LindaWorker.*/\%LindaWorker=$string/gI" "$inputfile"

                # check that the Linda nodes are correct
                lindaline=(`grep -i 'lindaworker' $inputfile`)
                if [[ $lindaline == *$string ]]
                then
                \techo "Using the correct nodes for Linda"
                else
                \techo "Using the wrong nodes for Linda"
                \techo "Nodes assigned by scheduler = $string"
                \techo "Line in Gaussian input file = $lindaline"
                \texit 1
                fi """ ))

        if queue == 'bf' or queue == 'ckpt':
            f.write(textwrap.dedent("""\
              \n
              # copy last log file to another name
              num=`ls -l %s*.log | wc -l`
              let "num += 1"
              cp %s.log %s$num.log""" % (f_input[i][0], f_input[i][0],
              f_input[i][0])))

        if 'gdv' in version:
            command = 'gdv'
        elif 'g16' in version:
            command = 'g16'
        else:
            command = 'g09'
        f.write(textwrap.dedent("""\
            \n
            # run Gaussian
            %s $inputfile

            exit 0 """ % (command)))

        print("""Please run 'sbatch %s' to submit to the scheduler\n"""
              % f_output[i])

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
def print_help():
    """Print a description of the script for the user."""

    print(textwrap.dedent("""\
        NAME
           \tGaussian-Submit
        DESCRIPTION
          \tThis program will help submit Gaussian calculations to
          \tthe Hyak Ikt and Mox clusters. The script will ask questions
          \tabout the calculation to help set up the .pbs/.sh script.

          \tIt will check for the types of nodes available to
          \tset appropriate defaults. It will also read your input
          \tfile to check for potential issues if using the STF
          \tallocation.
        EXAMPLES
          \tpython gaussian-sub.py input{.gjf,.com}
        AUTHOR
          \tPatrick J. Lestrange <patricklestrange@gmail.com>"""))
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
if __name__ == '__main__':
    get_user_input()
    check_Gaussian_input()
    write_slurm_script()
#----------------------------------------------------------------------------



