import os,sys
import multiprocessing as mp
import psutil
import shutil

from app.lib import getfiles,checkexists
import app.callDocker as cd

# trimming function
def q_trim(reads,jobs,cpu,outdir,tracker):
    minlength = 100
    windowsize = 4
    qscore = 30
    logfile = os.path.join(outdir,'qtrim.log')

    cmds = []
    read_path = ''
    for read_pair in reads:
        # main command
        main_cmd = 'java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {0}'.format(cpu)

        sid = os.path.basename(read_pair[0]).split('_')[0]

        args = '{read1} {read2} -baseout /output/{sid}.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'.format(minlength=minlength,windowsize=windowsize,qscore=qscore,read1=os.path.basename(read_pair[0]),read2=os.path.basename(read_pair[1]),sid=sid)

        cmds.append(main_cmd + args)

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read_pair[1]))
        elif read_path != os.path.dirname(os.path.abspath(read_pair[1])):
            print('Reads cannot be in multiple locations. Exiting.')
            sys.exit()
        else:
            pass
    checkexists(os.path.join(outdir,'trimmed'))

    # start multiprocessing
    pool = mp.Pool(processes=jobs)
    print('Begining quality trimming of reads:\n Number of Jobs: {0}\n CPUs/Job: {1}'.format(jobs,cpu))
    # denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Trimmomatic\n')
        # begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/trimmomatic:0.39',cmd,'/data',{read_path:'/data',os.path.join(outdir,'trimmed'):'/output'}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        # denote end of logs
        outlog.write('***********\n')

    # remove unpaired reads
    for root,dirs,files in os.walk(os.path.join(outdir,'trimmed')):
        for file in files:
            if 'U.fastq.gz' in file:
                os.remove(os.path.join(root,file))

    print('Finished quality trimming reads')

# assembly function
def assemble_reads(jobs,cpu_job,outdir):
    logfile = os.path.join(outdir,'assembly.log')
    input_path = os.path.join(outdir,'trimmed')

    # determine free ram
    free_ram = int(psutil.virtual_memory()[1]/1000000000)
    ram_job = int(free_ram / jobs)

    # assemble
    assemblies_path = os.path.join(outdir,'assemblies')
    checkexists(assemblies_path)

    # get trimmed reads
    fastqs = list(getfiles(input_path))[0]
    cmds = []
    read_path = ''
    for read_pair in fastqs:
        # main command
        read1 = os.path.basename(read_pair[0])
        read2 = os.path.basename(read_pair[1])
        sid = os.path.basename(read_pair[0]).split('_')[0]

        cmds.append('bash -c \"shovill --R1 /data/{0} --R2 /data/{1} --outdir /output/{2} --cpus {3} --ram {4} && mv /output/{2}/contigs.fa /output/{2}/{2}.fa\"'.format(read1,read2,sid,cpu_job,ram_job))

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read_pair[1]))
        elif read_path != os.path.dirname(os.path.abspath(read_pair[1])):
            print('Reads cannot be in multiple locations. Exiting.')
            sys.exit()
        else:
            pass

    # start multiprocessing
    pool = mp.Pool(processes=jobs)
    print('Begining assembly of reads:\n Number of Jobs: {0}\n CPUs/Job: {1}'.format(jobs,cpu_job))

    # denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Assembly\n')
        # begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/shovill:1.0.4',cmd,'/data',{read_path:'/data',os.path.join(outdir,'assemblies'):'/output'}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        # denote end of logs
        outlog.write('***********\n')
    print('Finished assembling reads')


# abricate function
def abricate(jobs,cpu_job,outdir):
    logfile = os.path.join(outdir,'abricate.log')
    input_path = os.path.join(outdir,'assemblies')
    # abricate
    abricated_path = os.path.join(outdir,'abricated')
    checkexists(abricated_path)

    # get assemblies
    fastas = list(getfiles(input_path))[0]
    # setup command list for annotating
    cmds = []
    for path in fastas:
        # main command
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        sid = os.path.basename(path).split('.')[0]
        cmds.append('abricate --db ncbi {0}/{1} > /output/{2}.tab'.format(assembly_dir,assembly_file,sid))

    # start multiprocessing
    pool = mp.Pool(processes=jobs)
    print('Begining abricate of assemblies:\n Number of Jobs: {0}\n CPUs/Job: {1}'.format(jobs,cpu_job))

    # denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Abricate\n')
        #begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/abricate',cmd,'/data',{input_path:'/data',os.path.join(outdir,'abricated'):'/output'}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')

    print('Finished running Abricate')

# ------------------------------------------------------

def pipeline(paired_reads,jobs,cpu_job,outdir,tracker):
    if not tracker.check_status('trimmed'):
        print('Starting trimming')
        print('Trimming reads with Trimmomatic')
        q_trim(paired_reads,jobs,cpu_job,outdir,tracker)
        tracker.update_status_done('trimmed')

    if not tracker.check_status('assemble'):
        print('Starting assembly')
        print('Assembling reads using Shovill')
        assemble_reads(jobs,cpu_job,outdir)
        tracker.update_status_done('assemble')

    if not tracker.check_status('abricate'):
        print('Starting Abricate')
        print('IDing AR genes with Abricate')
        abricate(jobs,cpu_job,outdir)
        tracker.update_status_done('abricate')
