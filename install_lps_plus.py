#!/usr/bin/env python3
USAGE='''
Usage:    python3 install_lps_plus.py [options]

Possible options:
    -p, --dir_of_lps_files   path to fLPS file dir. Default: [path_to_this_file]/flps
    -c, --dir_of_creb_files    path to creb file fir. Default: [path_to_this_file]/creb
    -n, --n6_dir    path to NBODY6++GPU dir. Default: ./installed_lps_plus/Nbody6PPGPU-beijing
    -r, --reb_dir    path to REBOUND dir. Default: ./rebound
    -d, --creb_dest    path to creb destination dir. Default: put beside n6_dir, named creb_production
    -h, --help

in case of any problem, see https://github.com/kaiwu-astro/lps_plus, or contact Kai Wu (kaiwu.astro@gmail.com)
'''

import os
import sys
import time
import shutil
from subprocess import run, check_call
from getopt import gnu_getopt as getopt

def get_output(cmd:str, raise_error=False) -> list:
    '''
    get output of a shell cmd. ignore error by default
    return: list of strings (one line per item, no \\n in the end)
    '''
    return str(run(cmd, shell=True, check=raise_error, capture_output=True).stdout, encoding='utf-8').split("\n")[:-1]

def printsep(num=1):
    for _ in range(num):
        print("-"*30)

def get_line_number(cmd):
    return int(str(run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').split(":")[0]) - 1

def clone_nbody6ppgpu(branch='dev'):
    check_call(f'git clone -b {branch} https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing {n6_dir}', shell=True)

def clone_rebound():
    # download and switch to tag v3.28.4
    check_call(f'git clone https://github.com/hannorein/rebound {reb_dir}', shell=True)

if __name__ == "__main__":
    my_dir = os.path.abspath(os.path.dirname(__file__))
    dir_of_lps_files = my_dir + os.sep + 'flps'
    dir_of_creb_files = my_dir + os.sep + 'creb'
    n6_dir = os.getcwd() + os.sep + 'installed_lps_plus' + os.sep + 'Nbody6PPGPU-beijing' # n6 is short for NBODY6++GPU 
    reb_dir = os.getcwd() + os.sep + 'rebound'
    creb_dest = ''

    if sys.version_info < (3, 7):
        raise RuntimeError("Python>=3.7 is required")

    opts, args = getopt(sys.argv[1:], 'n:r:p:c:d:h', ['n6_dir=', 'reb_dir=', 'dir_of_lps_files=', 'dir_of_creb_files=', 'creb_dest=', 'help'])
    for opt_name, opt_value in opts:
        if opt_name in ('-n', '--n6_dir'):
            n6_dir = opt_value
        if opt_name in ('-r', '--reb_dir'):
            reb_dir = opt_value
        if opt_name in ('-p', '--dir_of_lps_files'):
            dir_of_lps_files = opt_value
        if opt_name in ('-c', '--dir_of_creb_files'):
            dir_of_creb_files = opt_value
        if opt_name in ('-d', '--creb_dest'):
            creb_dest = opt_value
        if opt_name in ('-h', '--help'):
            print(USAGE)
            sys.exit(0)
    if creb_dest == '':
        creb_dest = os.path.dirname(n6_dir) + os.sep + 'creb_production'

    # print these args; wait for 3s and continue
    print('Your input: \n  NBODY6++GPU path: ', n6_dir, '\n  REBOUND path: ', reb_dir, '\n  fLPS files path: ', dir_of_lps_files, '\n  CREB files path: ', dir_of_creb_files)
    time.sleep(1)

    # get NBODY6++GPU and REBOUND
    if not os.path.isdir(n6_dir):
        os.makedirs(n6_dir)
        clone_nbody6ppgpu()
    else:
        print('NBODY6++GPU path already exists. Skip downloading.')
    ## check
    if not os.path.isfile(n6_dir + '/build/Makefile.in'):
        raise IOError("Wrong NBODY6++GPU path because did not find build/Makefile.in inside. Use -n to specify NBODY6++GPU path, or remove the NBODY6++GPU directory and try again.")
    
    if not os.path.isdir(reb_dir):
        os.makedirs(reb_dir)
        clone_rebound()
    else:
        print('REBOUND path already exists. Skip downloading.')
    print('Switching REBOUND to version 3.28.4 ...')
    check_call(f'cd {reb_dir} && git checkout 3.28.4', shell=True)

    # 使用sed命令，查找dir_of_creb_files + 'Makefile' 文件内开头为REB_DIR=的行，修改为REB_DIR=reb_dir
    printsep(3)
    print(f'Copying creb files to {creb_dest} ...')
    if os.path.isdir(creb_dest):
        input(f"Directory {creb_dest} already exists. Press Enter to remove it and continue, or Ctrl+C to exit.")
        shutil.rmtree(creb_dest)
    shutil.copytree(dir_of_creb_files, creb_dest)
    print(f'Editing Makefile of creb to set REB_DIR= to {reb_dir} ...')
    check_call(f"sed -i 's#REB_DIR=.*#REB_DIR={reb_dir}#' {creb_dest}/Makefile", shell=True)

    printsep(3)
    print(f'Patching NBODY6++GPU {n6_dir} with fLPS...')
    if os.path.isfile(n6_dir + '/src/Main/fLPS.f') and os.path.isfile(n6_dir + '/include/fLPS.h'):
        input(r"fLPS files already exist. Press any key to reset to original NBODY6++GPU, or Ctrl+C to abort.")
        check_call(f'cd {n6_dir} && git reset --hard HEAD && git clean -fd', shell=True)
    print(f"Source dir: {dir_of_lps_files}")
    if not os.path.isfile(dir_of_lps_files + '/fLPS_in.f') or not os.path.isfile(dir_of_lps_files + '/fLPS_in.h'):
        raise IOError(f"fLPS files not found in {dir_of_lps_files}. Use -p to specify the directory of fLPS files")
    # 1. add files
    printsep(3)
    print("1. add fLPS files")
    shutil.copy2(dir_of_lps_files + '/fLPS_in.f', n6_dir + '/src/Main/fLPS.f')
    shutil.copy2(dir_of_lps_files + '/fLPS_in.h', n6_dir + '/include/fLPS.h')
    if os.path.isfile(n6_dir + '/src/Main/fLPS.f') and os.path.isfile(n6_dir + '/include/fLPS.h'):
        print("fLPS.h and fLPS.f copied")
    else:
        raise OSError("Copy file failed. Check directory permissions of NBODY6++GPU")

    # 2. modify files for flps
    # 修改文件的方法模拟人工修改：找到关建行，然后添加代码
    #   a. build/Makefile.in: SOURCE = 一行，最后加入fLPS.f
    printsep(3)
    print("2. insert code to original files for flps")
    print("Modified files are backed up to FILENAME.flpsbackup")
    print("2. a. build/Makefile.in")
    subdir = "/build/Makefile.in"
    with open(n6_dir + subdir, 'r') as f:
        fcontent = f.readlines()
    shell_cmd = 'grep -in "\\bSOURCE =" ' + n6_dir + subdir
    source_file_line_number = get_line_number(shell_cmd)
    for line_number in range(source_file_line_number, len(fcontent)):
        if fcontent[line_number].isspace() == False and fcontent[line_number].split()[-1][-1] != '\\':
            # 找到第一个不是空行而且结尾不为反斜杠（续行符号）的行，就是SOURCE的最后一行
            break
    if line_number == len(fcontent):
        raise OSError("Did not find SOURCE in " + subdir)
    print('Inserting "fLPS.f" to ' + subdir)
    fcontent[line_number] = fcontent[line_number][:-1] + ' fLPS.f \n'

    os.rename(n6_dir + subdir, n6_dir + subdir + '.flpsbackup')
    with open(n6_dir + subdir, 'w') as f:
        f.write(''.join(fcontent))

    print("Result given by diff:")
    print(str(run("diff -u " + n6_dir + subdir + '.flpsbackup' + " " +
                        n6_dir + subdir, shell=True, capture_output=True).stdout, 
            encoding='utf-8'))

    #   b. include/common6.h: 开头加入 include 'fLPS.h'
    printsep(3)
    print("2. b. include/common6.h")
    subdir = "/include/common6.h"
    with open(n6_dir + subdir, 'r') as f:
        fcontent = f.readlines()
    shell_cmd = 'grep -in "mpi_base" ' + n6_dir + subdir
    include_mpibase_linenumber = get_line_number(shell_cmd)
    print("Inserting INCLUDE 'fLPS.h' to " + subdir)
    line_to_insert = " "*6 + "INCLUDE 'fLPS.h'\n"

    os.rename(n6_dir + subdir, n6_dir + subdir + '.flpsbackup')
    with open(n6_dir + subdir, 'w') as f:
        f.write(''.join(fcontent[:include_mpibase_linenumber]) +
                line_to_insert + ''.join(fcontent[include_mpibase_linenumber:]))

    print("Result given by diff:")
    print(str(run("diff -u " + n6_dir + subdir + '.flpsbackup' + " " +
                        n6_dir + subdir, shell=True, capture_output=True).stdout, 
            encoding='utf-8'))
    
    #   c. src/Main/intgrt.f: line 139附近，CALL custom_output之前：加入代码段
    #   d. src/Main/intgrt.f: line 1276附近，CALL custom_output之前：加入 IF(KZ(43).GE.1) call LPSOutput
    #      实际实现方法：找custom_output。第一个结果作为init，后面的结果都拿来做输出
    printsep(3)
    print("2. c.&d. src/Main/intgrt.f")
    subdir = "/src/Main/intgrt.f"
    if not os.path.isfile(n6_dir + subdir):
        subdir = "/src/Main/intgrt.F"
    with open(n6_dir + subdir, 'r') as f:
        fcontent = f.readlines()
    shell_cmd = 'grep -in "custom_output" ' + n6_dir + subdir
    custom_output_found = str(run(shell_cmd, shell=True, capture_output=True).stdout, encoding='utf-8').split('\n')[:-1] #舍弃最后一个空字符串
    print(f"Found {len(custom_output_found)} custom_output in " + subdir)
    if len(custom_output_found) != 2:
        print("WARNING: more custom_output found in file")
    line_numbers = []
    line_nspace = []
    for ifound in range(len(custom_output_found)):
        # 对于call custom_output上一行/两行是if内的情况，把行号移到if前
        custom_line = int(custom_output_found[ifound].split(":")[0]) - 1
        if fcontent[custom_line-1].split() != []:
            if fcontent[custom_line-1].split()[0][0:2].lower() == 'if':
                custom_line -= 1
        elif fcontent[custom_line-2].split() != []:
            if fcontent[custom_line-2].split()[0][0:2].lower() == 'if':
                custom_line -= 2
        line_numbers.append(custom_line)
        nspace = 0
        while fcontent[line_numbers[-1]][nspace].isspace():
            nspace += 1
        line_nspace.append(nspace)

    code_block_to_insert = [
        "*" + " "*(line_nspace[0]-1) + "fLPS Output init start\n",
            " "* line_nspace[0]    + " IF(KZ(43).GE.1) THEN\n",
            " "* line_nspace[0]    + "     CALL LPSinit\n",
            " "* line_nspace[0]    + "     CALL LPSOutput\n",
            " "* line_nspace[0]    + " END IF\n",
        "*" + " "*(line_nspace[0]-1) + "fLPS Output init end\n"]
        
    os.rename(n6_dir + subdir, n6_dir + subdir + '.flpsbackup')
    with open(n6_dir + subdir, 'w') as f:
        f.write(
            ''.join(fcontent[:line_numbers[0]]) + 
            ''.join(code_block_to_insert)
        )
        for ifound in range(1, len(custom_output_found)):
            f.write(
                ''.join(fcontent[line_numbers[ifound-1] : line_numbers[ifound]]) + 
                " " * line_nspace[ifound] + "IF(KZ(43).GE.1) call LPSOutput\n"
            )
        f.write(''.join(fcontent[line_numbers[-1] : ]))
        
    print("Result given by diff:")
    print(str(run("diff -u " + n6_dir + subdir + '.flpsbackup' + " " +
                        n6_dir + subdir, shell=True, capture_output=True).stdout, 
            encoding='utf-8'))
    
    #    2. e. src/Main/escape.F: line 323附近，call flush(11)之前：加入 IF(KZ(43).GE.1) call RecordEscaper(NAMEI)
    printsep(3)
    print("2. c.&d. src/Main/escape.F")
    subdir = "/src/Main/escape.F"
    if not os.path.isfile(n6_dir + subdir):
        subdir = "/src/Main/escape.f"
    with open(n6_dir + subdir, 'r') as f:
        fcontent = f.readlines()

    shell_cmd = 'grep -in "flush(11)" ' + n6_dir + subdir
    flush11_found = str(run(shell_cmd, shell=True, capture_output=True).stdout, encoding='utf-8').split('\n')[:-1] #舍弃最后一个空字符串
    print(f"Found {len(flush11_found)} flush(11) in " + subdir)
    if len(flush11_found) > 1:
        print("WARNING: more than one flush(11) found in escape.F . This file is potentially modified.")
    line_numbers = []
    for ifound in range(len(flush11_found)):
        # 在flush(11)后面立即刷新escaper List
        flush11_line = int(flush11_found[ifound].split(":")[0]) - 1
        line_numbers.append(flush11_line + 1)

    os.rename(n6_dir + subdir, n6_dir + subdir + '.flpsbackup')
    with open(n6_dir + subdir, 'w') as f:
        f.write(
            ''.join(fcontent[:line_numbers[0]]) + 
            "      IF(KZ(43).GE.1) call RecordEscaper(NAMEI)\n"
        )
        for ifound in range(1, len(flush11_found)):
            f.write(
                ''.join(fcontent[line_numbers[ifound-1] : line_numbers[ifound]]) + 
                "      IF(KZ(43).GE.1) call RecordEscaper(NAMEI)\n"
            )
        f.write(''.join(fcontent[line_numbers[-1] : ]))

    print("Result given by diff:")
    print(str(run("diff -u " + n6_dir + subdir + '.flpsbackup' + " " +
                        n6_dir + subdir, shell=True, capture_output=True).stdout, 
            encoding='utf-8'))

    printsep()
    printsep()
    print("Succeed! fLPS has been patched to NBODY6++GPU, and creb is ready for planetary simulations.\nWhat's next? See https://github.com/kaiwu-astro/lps_plus")
