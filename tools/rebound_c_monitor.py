#!/usr/bin/env python3

import os
import sys
import time
import signal
import psutil
import zipfile
from glob import glob
from tqdm.auto import tqdm
from getopt import gnu_getopt as getopt
from subprocess import call, Popen, check_output
from watch_memory import start_watch_memory
# from lps_tools import *
from lps_tools import userhome, get_starids_from_rebInfo, get_starid_from_n6_out_filename, tee, restore_host_ptb_pkl_to_txt, wait_io_and_mem, confirm_ctrl_c, run_from_ipython, get_output, get_nplanet_nkuiper_from_crebound_input
my_path = os.path.abspath(__file__)
my_dir = os.path.dirname(my_path)
RootDir = os.path.dirname(my_dir)

############ USER OPTIONS #############
# user options can also input from the command line argument
## how many logical cores for the simulation? os.cpu_count() means all 
MAX_NPROCESS = os.cpu_count() - 2 
## whether to plot the results when simulation finishes
DOPLOT = False 
## how many perturbers in each system
nuseptb = 10
## path to the n6 result directory (with Host*.txt, Perturber*.txt, LPSdiag.txt, and the standard output named NB.out.txt)
n6_result_dir = ''
## path to the rebound problem file directory
PROBLEM_DIR = RootDir + '/Source/problemfile_rebound'
## name of the executable file. It can be found in the Makefile under PROBLEM_DIR. Default: "problem"
REB_EXEC = "problem"
## path to the input file
input_file_path = PROBLEM_DIR + '/input_P1.txt'
# compile with -g and leave executables without running. In formal run: False
CDEBUG = False           
############ USER OPTIONS END #########


def getnrunning():
    cmdreslst = get_output(r"ps -ax | grep problem | awk '{print $5}'")
    nrunning = 0
    for cmdres in cmdreslst:
        if cmdres.startswith("./problem"):
            nrunning += 1
    return nrunning

def get_MAX_NPROCESS():
    '''可以用来刷新，也可以取返回值'''
    global MAX_NPROCESS
    try:
        with open(userhome + "/creconf.txt", "r") as f:
            max_process_str = f.readlines()
        max_process_num = int(max_process_str[0].strip())
        MAX_NPROCESS = max_process_num
    except Exception:
        pass
    return MAX_NPROCESS

def generate_example_input_file(filepath=my_dir+'/example_input_file.txt'):
    with open(filepath, 'w') as f:
        # Write file format description
        f.write("# Example Input File for crebound\n")
        f.write("# Format: key=value, separated by 1 or more spaces\n")
        f.write("# Lines with leading '#' or spaces are comments\n")
        f.write("# Available keys: m_me, a_au, e, i_deg, peri_deg, node_deg, M_deg\n")
        f.write("# m_me: mass of the particle in earth mass\n")
        f.write("# a_au: semi-major axis of the particle in au\n")
        f.write("# e: eccentricity of the particle\n")
        f.write("# i_deg: inclination of the particle in degrees\n")
        f.write("# peri_deg: argument of pericenter of the particle in degrees\n")
        f.write("# node_deg: longitude of the ascending node of the particle in degrees\n")
        f.write("# M_deg: mean anomaly of the particle in degrees\n")
        f.write("# peri, node, and M can also be the word: random, which randomize between 0 to 360 degrees.\n\n\n")
        f.write("# Some examples particles below:\n")
        f.write("# 50-au-Jupiter\n")
        f.write("m_me=317.83  a_au=50.0  e=0.04839266  i_deg=0.0  peri_deg=14.75385  node_deg=100.55615  M_deg=34.40438\n")
        f.write("# Planet earth\n")
        f.write("m_me=1.0  a_au=1.0  e=0.0167  i_deg=0.0  peri_deg=102.93735  node_deg=-11.26064  M_deg=100.46435\n")
        f.write("# 3 comets (set mass=0 to greatly boost the simulation speed)\n")
        f.write("m_me=0.0  a_au=20.0  e=0.01  i_deg=0.01  peri_deg=random  node_deg=random  M_deg=random\n")
        f.write("m_me=0.0  a_au=25.0  e=0.01  i_deg=0.01  peri_deg=random  node_deg=random  M_deg=random\n")
        f.write("m_me=0.0  a_au=100.0  e=0.01  i_deg=0.01  peri_deg=random  node_deg=random  M_deg=random\n")
    print(f"Example input file generated at: {filepath}")

USAGE = f'''
Usage: python3 {sys.argv[0]} [options]

Prerequisites:
    1. REBOUND (https://github.com/hannorein/rebound) has been git cloned, and uses version == 3.x . Version >= 4 is not supported due to the change of the API.
        e.g., git clone https://github.com/hannorein/rebound && cd rebound && git checkout 3.28.4
        Alternatively, download https://github.com/hannorein/rebound/archive/refs/tags/3.28.4.tar.gz and extract

    2. The 'problem file', i.e., the C code that calls the REBOUND C API and starts a single planetary system simulation, has been prepared, and put beside the cloned REBOUND directory. 
        e.g., directory structure:
            ├── problemfile_rebound
            │   ├── Makefile
            │   ├── problem
            │   ├── problem.c
            │   ├── problem.h
            ├── rebound
            │   ├── examples
            │   ├── ipython_examples
            │   ├── rebound
            │   ├── LICENSE
            │   ├── README.md
            │   ├── ...
        You may consider using --problem-dir to specify the path to the problem file directory.
        (also possible not to keep this directory structure, but need to edit the Makefile in the problem file directory)

    3. The final executable should be named "problem" in Makefile.

Options:
    -u, --nuseptb
        How many perturbers in each system. Default: 10
    -i, --input-file
        Path to the input file. Default: {input_file_path}
    -d, --dir
        Path to the n6 result directory (with Host*.txt, Perturber*.txt, LPSdiag.txt, and the standard output named NB.out.txt)
    -y, --yes
        Delete the existing reb output directory without asking. Default: ask
    -n, --no
        Continue without deleting the existing reb output directory. Default: ask
    -r, --reverse
        Reverse the order of starids to run. Default: False
    -p, --max-process
        How many logical cores for the simulation? os.cpu_count() means all. Default: os.cpu_count() - 2
    -m, --problem-dir
        Path to the rebound problem file directory. Default: {PROBLEM_DIR}
    --exec-name
        Name of the executable file. It can be found in the Makefile under PROBLEM_DIR. Default: {REB_EXEC}
    -g, --cdebug
        Compile with -g and leave executables without running. Default: False
    --example-input
        Generate an example input file and exit
    -h, --help
        Print this help message and exit
Example:
    python3 {sys.argv[0]} --nuseptb=10 --dir=/path/to/n6_result_dir --input-file=/path/to/input_file.txt --problem-dir=/path/to/problem_file_dir
'''

if __name__ == '__main__':
    if not run_from_ipython():
        delete_yes_no = None
        reverse_starid_order = False
        input_filename = os.path.basename(input_file_path)
        opts, args = getopt(sys.argv[1:], 'u:i:d:ynrp:m:gh', ['nuseptb=', 'input-file=', 'dir=', 'yes', 'no', 'reverse', 'max-process=', 'problem-dir=', 'exec-name=', 'cdebug', 'example-input', 'help'])
        for opt_name, opt_value in opts:
            if opt_name in ('-u', '--nuseptb'):
                nuseptb = int(opt_value)
            if opt_name in ('-i', '--input-file'):
                input_file_path = opt_value
                input_filename = os.path.basename(input_file_path)
            if opt_name in ('-d', '--dir'):
                n6_result_dir = opt_value
            if opt_name in ('-y', '--yes'):
                delete_yes_no = 'Y'
            if opt_name in ('-n', '--no'):
                delete_yes_no = 'n'
            if opt_name in ('-r', '--reverse'):
                reverse_starid_order = True
            if opt_name in ('-p', '--max-process'):
                MAX_NPROCESS = int(opt_value)
            if opt_name in ('-m', '--problem-dir'):
                PROBLEM_DIR = opt_value
            if opt_name in ('--exec-name',):
                REB_EXEC = opt_value
            if opt_name in ('-g', '--cdebug'):
                CDEBUG = True
            if opt_name in ('--example-input',):
                generate_example_input_file()
                sys.exit(0)
            if opt_name in ('-h', '--help'):
                print(USAGE)
                sys.exit(0)
        if not os.path.isdir(PROBLEM_DIR):
            print(f"ERROR: Problem file directory {PROBLEM_DIR} does not exist. Use --problem-dir to specify the path to the problem file directory.")
            print(USAGE)
            sys.exit(1)
        if not os.path.isdir(n6_result_dir):
            print(f"ERROR: n6 result directory {n6_result_dir} does not exist. Use --dir to specify the path to the n6 result directory.")
            print(USAGE)
            sys.exit(1)
        start_watch_memory()

    # print all readed arguments
    print("nuseptb:", nuseptb)
    print("input_file_path:", input_file_path)
    print("input_filename:", input_filename)
    print("n6_result_dir:", n6_result_dir)
    print("delete_yes_no:", delete_yes_no)
    print("reverse_starid_order:", reverse_starid_order)
    print("MAX_NPROCESS:", MAX_NPROCESS)
    print("PROBLEM_DIR:", PROBLEM_DIR)
    print("REB_EXEC:", REB_EXEC)
    print("CDEBUG:", CDEBUG)
    
    nplanet, nkuiper = get_nplanet_nkuiper_from_crebound_input(input_file_path)
    print("nplanet:", nplanet)
    print("nkuiper:", nkuiper)
    reb_dirname = f'reb{nplanet}{nuseptb:02d}{nkuiper:03d}'
    reb_exec_dest = REB_EXEC + "_" + reb_dirname # 可执行文件位于n6目录下，不在reb结果目录下

    # Compile rebound
    if CDEBUG:
        print("......DEBUGGING..................")
        call(f"cd {PROBLEM_DIR} && rm -f {REB_EXEC} ; export CFLAGS='-g -O0 -fsanitize=address -Werror=array-bounds' && make", shell=True)
    else:
        call(f"cd {PROBLEM_DIR} && rm -f {REB_EXEC} ; make", shell=True)

    def Popen_exit(signum, frame):
        for process in processlist:
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass
        print("Interrupted. All rebound subprocesses killed.")
        sys.exit(0)
    # kill all simulation processes when this script is killed with Ctrl+C or SIGTERM
    signal.signal(signal.SIGTERM, Popen_exit) 
    confirm_ctrl_c(ctrl_c_handler=Popen_exit)

    processlist = []
    print(f"{n6_result_dir} start ...")
    check_output(f"cp -a {PROBLEM_DIR}/{REB_EXEC} {n6_result_dir}/{reb_exec_dest}", shell=True)
    check_output(f"cp -a {input_file_path} {n6_result_dir}/", shell=True)
    if not os.path.isfile(n6_result_dir + '/librebound.so'):
        check_output("cp " + PROBLEM_DIR + "/lib* " + n6_result_dir + '/', shell=True)
    def _delete():
        print("Deleting...")
        call("rm -rf "+ n6_result_dir + '/' + reb_dirname, shell=True)
        os.makedirs(n6_result_dir + '/' + reb_dirname)
    def _continue():
        print("Continue without delete.")
    if os.path.isdir(n6_result_dir + '/' + reb_dirname):
        if delete_yes_no is None:
            while 1:
                inp = input('-'*30 + f"\n{n6_result_dir + '/' + reb_dirname} exists. Delete? \n[Y] to delete \n[N/n] continue without delete \n[Ctrl+C] to abort: ")
                if inp in 'Yy':
                    delete_yes_no = 'Y'
                    break
                elif inp in 'Nn':
                    delete_yes_no = 'n'
                    break
                else:
                    print("Wrong input.")
        else: # specified by command line argument
            pass 
        if delete_yes_no in 'Yy':
            _delete()
        elif delete_yes_no in 'Nn':
            _continue()
    else:
        os.makedirs(n6_result_dir + '/' + reb_dirname)

    with zipfile.ZipFile(n6_result_dir + "/" + reb_dirname + "/zippedReb.zip", 'w') as z:
        z.write(PROBLEM_DIR + "/problem.c")
        z.write(PROBLEM_DIR + "/problem.h")
        z.write(PROBLEM_DIR + "/Makefile")
        z.write(PROBLEM_DIR + "/librebound.so")
        z.write(my_path)

    hostpaths = glob(n6_result_dir + "/*Host*")
    if reverse_starid_order:
        hostpaths = list(reversed(hostpaths))

    starids_finished = get_starids_from_rebInfo(n6_result_dir + '/' + reb_dirname, DISCARD_LARGE_DE=False) # once a simulation finished, it will write info to rebinfo.txt
    for hostpath in tqdm(hostpaths):
        if '.pkl' in hostpath: 
            if (hostpath.split('.pkl')[0] + '.txt') in hostpaths: # optimized, txt expanded: wait to use txt
                continue
            else: # optimized, no txt: expand to txt and run
                restore_host_ptb_pkl_to_txt(hostpath)
                restore_host_ptb_pkl_to_txt(hostpath.replace("Host", "Perturber"))

        starid = get_starid_from_n6_out_filename(hostpath)
        if int(starid) in starids_finished:
            print("Already finished...")
            continue

        wait_io_and_mem()
        while getnrunning() >= get_MAX_NPROCESS():
            time.sleep(10)
        
        # ./problem --n6-result-dir=/s16/planet_in_GC_data/8k-fixed/100.0_8000_0.5_3.0_0.0_0_0_1e-06_0/1 --starname=89 --nuseptb=10 --input-file-path=input_P1.txt --reb-output-subdir=ttt
        cmd = f"cd {n6_result_dir} && nice -n 2 ./{reb_exec_dest} --n6-result-dir={n6_result_dir} --starname={starid} --nuseptb={nuseptb} --reb-output-subdir={reb_dirname} --input-file-path={input_filename} > {n6_result_dir}/{reb_dirname}/reb_out_{starid}_consoleStdout.txt 2>&1"
        if CDEBUG:
            print(f"cd {n6_result_dir} && gdb -ex run --args ./{reb_exec_dest} --n6-result-dir={n6_result_dir} --starname={starid} --nuseptb={nuseptb} --reb-output-subdir={reb_dirname} --input-file-path={input_filename}")
            print('Please start debug manually.')
            sys.exit(0)
        call(f"python3 {my_dir}/save_disk_space.py -d {n6_result_dir}/{reb_dirname}", shell=True)
        processlist.append(Popen(cmd, shell=True, preexec_fn=os.setsid))
        time.sleep(1)


    for process in processlist:
        process.wait()
    call(f"python3 {my_dir}/save_disk_space.py -d {n6_result_dir}", shell=True)

    # statistics
    OLDPWD = os.getcwd()
    os.chdir(n6_result_dir + '/' + reb_dirname)
    tee('-'*70 + f"\nin {os.getcwd()}: ", filename='simulogReb.txt')
    tee("Total number of planetary systems: ", get_output("ls *Energy* -l | wc -l")[0], filename='simulogReb.txt')
    tee("Aborted because of nan/inf: ", get_output("grep nan *console* | wc -l")[0], filename='simulogReb.txt')
    tee("Reported 'At least 10 predictor corrector loops in IAS15 did not converge': ",
        get_output("grep least *console* | wc -l")[0], filename='simulogReb.txt')
    n_total_de_too_large = int(get_output("grep Too *consoleStdout* | wc -l")[0])
    tee("Total energy error too large: ", n_total_de_too_large, filename='simulogReb.txt')
    # tee("Unfinished and energy error too large: ", int(get_output(
    #     "grep -P '\s+.*\+[0-9][1-9]' *EnergyError* | wc -l")) - n_total_de_too_large, filename='simulogReb.txt')
    os.chdir(OLDPWD)

    foldersize_GB = int(check_output("du "+n6_result_dir+"/"+reb_dirname+" --max-depth=0", shell=True).split()[0]) / 1024 / 1024
    if foldersize_GB > 1.0 :
        tee(f"\nSize of reb output in {n6_result_dir} is {foldersize_GB:.2f} GB. Too large.", filename='simulogReb.txt')
    if DOPLOT:
        print("-"*30+"   Plotting   "+"-"*30)
        call(f"cd {RootDir}; trash 0scripts/__pycache__; python3 0scripts/newplotreb.py --reanalyze --dir={n6_result_dir}/{reb_dirname}", shell=True)
    
    print("\n".join(get_output('df -lh | grep -P -v "/dev/loop|tmpfs|udev"')))
