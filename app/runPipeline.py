import os,sys
import multiprocessing as mp
import psutil
import shutil

from app.lib import getfiles,checkexists
import app.callDocker as cd

# trimming function
def q_trim(paired_reads,jobs,cpu,outdir,tracker):
    minlength = 100
    windowsize = 4
    qscore = 30
    logfile = os.path.join(outdir,'qtrim.log')

    # trimming output path
    trimmed_path = os.path.join(outdir,'trimmed')
    checkexists(trimmed_path)

    # setup command list for Trimmomatic
    cmds = []
    read_path = ''
    for read_pair in paired_reads:
        main_cmd = f'java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {cpu}'
        sid = os.path.basename(read_pair[0]).split('_')[0]
        read1 = os.path.basename(read_pair[0])
        read2 = os.path.basename(read_pair[1])

        # main command
        args = f' {read1} {read2} -baseout /output/{sid}.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
        cmds.append(main_cmd + args)

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read2))
        elif read_path != os.path.dirname(os.path.abspath(read2)):
            print('Reads cannot be in multiple locations. Exiting.')
            sys.exit()
        else:
            pass

    # start multiprocessing
    pool = mp.Pool(processes=jobs)
    print(f'Begining quality trimming of reads:\n Number of Jobs: {jobs}\n CPUs/Job: {cpu}')

    # denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Trimmomatic\n')
        # run commands and get results
        results = pool.starmap_async(cd.call,[['staphb/trimmomatic:0.39',cmd,'/data',{read_path:'/data',trimmed_path:'/output'}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
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

    # assembly output path
    assemblies_path = os.path.join(outdir,'assemblies')
    checkexists(assemblies_path)

    # get trimmed reads
    fastqs = list(getfiles(input_path))[0]

    # setup command list for shovill
    cmds = []
    read_path = ''
    for read_pair in fastqs:
        # main command
        read1 = os.path.basename(read_pair[0])
        read2 = os.path.basename(read_pair[1])
        sid = os.path.basename(read_pair[0]).split('_')[0]

        cmds.append(f'bash -c \"shovill --R1 /data/{read1} --R2 /data/{read2} --outdir /output/{sid} --cpus {cpu_job} --ram {ram_job} && mv /output/{sid}/contigs.fa /output/{sid}/{sid}.fa\"')

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read_pair[1]))
        elif read_path != os.path.dirname(os.path.abspath(read_pair[1])):
            print('Reads cannot be in multiple locations. Exiting.')
            sys.exit()
        else:
            pass

    # start multiprocessing
    pool = mp.Pool(processes=jobs)
    print(f'Begining assembly of reads:\n Number of Jobs: {jobs}\n CPUs/Job: {cpu_job}')

    # denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Assembly\n')
        # run commands and get results
        results = pool.starmap_async(cd.call,[['staphb/shovill:1.0.4',cmd,'/data',{read_path:'/data',os.path.join(outdir,'assemblies'):'/output'}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        outlog.write('***********\n')
    print('Finished assembling reads')


# abricate function
def abricate(jobs,cpu_job,outdir):
    input_path = os.path.join(outdir,'assemblies')
    abricate_path = os.path.join(outdir,'abricate')
    checkexists(abricate_path)

    # get assemblies
    fastas = list(getfiles(input_path))[0]

    # setup command list for abricate
    cmds = []
    for path in fastas:
        # main command
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        sid = os.path.basename(path).split('.')[0]
        cmds.append(f'abricate --db ncbi {assembly_dir}/{assembly_file}')

        # start multiprocessing
        pool = mp.Pool(processes=jobs)
        print('Beginning abricate')

        # output
        handle = sid + '.tab'
        logfile = os.path.join(abricate_path,handle)
        with open(logfile,'a') as outlog:
            # run commands and get results
            results = pool.starmap_async(cd.call,[['staphb/abricate',cmd,'/data',{input_path:'/data',abricate_path:'/output'}] for cmd in cmds])
            stdouts = results.get()
            for stdout in stdouts:
                outlog.write(stdout)
        cmds = []
        print('Finished running Abricate')

# quast function
def quast(ref,jobs,cpu_job,outdir):
    logfile = os.path.join(outdir,'quast.log')
    input_path = os.path.join(outdir,'assemblies')
    ref_path = os.path.dirname(ref)
    ref_genome = os.path.basename(ref)
    print(ref_path)
    print(ref_genome)
    # determine free ram
    free_ram = int(psutil.virtual_memory()[1]/1000000000)
    ram_job = int(free_ram / jobs)

    # quast output path
    quast_path = os.path.join(outdir,'quast')
    checkexists(quast_path)

    # get assemblies
    fastas = list(getfiles(input_path))[0]

    # setup command list and assembly list string
    cmds = []
    fasta_files = ''

    for path in fastas:
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        sid = os.path.basename(path).split('.')[0]
        fasta_files = fasta_files + f'/data/{assembly_dir}/{assembly_file} '

    # main command
    cmds.append(f'quast.py -t {cpu_job} -o /output/ -r /reference/{ref_genome} '+ fasta_files)
    print(cmds)
    # start multiprocessing
    pool = mp.Pool(processes=jobs)
    print('Begining Quast:\n Number of Jobs: {0}\n CPUs/Job: {1}'.format(jobs,cpu_job))

    # denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Assembly\n')
        # run commands and get results
        results = pool.starmap_async(cd.call,[['staphb/quast',cmd,'/data',{input_path:'/data',ref_path:'/reference',quast_path:'/output'}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        outlog.write('***********\n')
    print('Finished running quast')


# ------------------------------------------------------

def spriggan_pipeline(paired_reads,ref,jobs,cpu_job,outdir,tracker):
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

    if not tracker.check_status('quast'):
        print('Starting Quast')
        print('Evaluating assembly quality with Quast')
        quast(ref,jobs,cpu_job,outdir)
        tracker.update_status_done('quast')
