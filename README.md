
```
 __       ____    ____       __      
/\ \     /\  _`\ /\  _`\    /\ \     
\ \ \    \ \ \L\ \ \,\L\_\  \_\ \___ 
 \ \ \  __\ \ ,__/\/_\__ \ /\___  __\
  \ \ \L\ \\ \ \/   /\ \L\ \/__/\ \_/
   \ \____/ \ \_\   \ `\____\  \ \_\ 
    \/___/   \/_/    \/_____/   \/_/ 
```
# Welcome to *LPS+* !
**LPS+** is a code to bridge NBODY6++GPU, a star cluster *N*-body simulator, and REBOUND, a planetary system *N*-body simulator. The code is an upgrade to the LPS (LonelyPlanetS) code ([Cai et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.4337C); [Flammini Dotti et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.2280F)).

**LPS+** contains mainly five files:
- `fLPS.f`: an additional routine for NBODY6++GPU, to 
- 1. select planetary system host stars, 
- 2. identify their perturbing stars (hereafter, perturbers) from the neighbour list (variable `LIST` of NBODY6++GPU) during the star cluster simulation, and
- 3. output the motion of host stars and perturbers with high time-resolution.
- `fLPS.h`: FORTRAN "header" file, to put variables into the NBODY6++GPU common block
- `creb.c`: a "problem file" of REBOUND, to call REBOUND via its `C` API [link to REBOUND API documentation](https://rebound.readthedocs.io/en/3.28.4/api/), to 
- 1. read initial conditions (hereafter, ic) of planetary systems, 
- 2. run the planetary system simulation with the perturbing stars, and
- 3. *reset* the motion of the perturbers at intervals.
- `Makefile` of the "problem file"
- `creb.h` a header file for `creb.c`

# Dependency

- NBODY6++GPU, version Dec2019 or higher: https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing/tree/dev 
- REBOUND, 3.0 <= version <= 3.28.4: https://github.com/hannorein/rebound/tree/3.28.4
- A Linux (*nix) operating system (All previous tests are performed on Linux; MacOS may or may not work)
- C compiler and FORTRAN compiler (which are probably already installed in most Linux and MacOS)
- Python >= 3.7 (for the one-key installation script)

# Installation

Easily install **LPS+** with the `install_lps_plus.py` script! 

1.  Clone this repo: `git clone https://github.com/kaiwu-astro/lps_plus`
2.  ```python3  install_lps_plus.py```

The script will download NBODY6++GPU & REBOUND and patch them with **LPS+**. It is also possible to use your downloaded NBODY6++GPU & REBOUND version to install **LPS+**. Check the options with `python3  install_lps_plus.py --help`.


# Usage 

## Prerequisite

In the following instructions, I assume that:
- You have some experience with NBODY6++GPU. Perhaps you need to successfully (i) download the code, (ii) compile the code, and (iii) run a star cluster simulation at least once. (See https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing)
- You have some experience with REBOUND. Maybe at a similar level of proficiency as above. (See https://github.com/hannorein/rebound)

I am sorry for the prerequisite, but simulating planetary systems in star clusters is a bit of cutting-edge science 🎩 . 

## Compile NBODY6++GPU 

Follow instructions in https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing .

TL;DR: `./configure --enable-mcmodel=large --with-par=b1m --disable-gpu --disable-mpi && make clean && make -j`

## Compile creb.c

Go to the directory of installed creb (by default: `./installed_lps_plus/creb_production`), and simply type `make`

## Run a simulation 

### Run the patched NBODY6++GPU
1. Set KZ(43)=1 in the initial condition file to enable fLPS output, and run as normal.
2. By default, fLPS searches 100 stars with mass closest to 1 solar mass as planetary system host star, and for each of them, search 10 nearest neighbor stars as perturbers. The output files are LPS_Host_[star_name].txt x100, LPS_Perturber_[star_name].txt x100, and LPSdiag.txt x1. If you want to change any parameter related to host stars and perturbers, modify {n6_dir}/src/Main/fLPS.f and {n6_dir}/include/fLPS.h, and recompile NBODY6++GPU.
3. How to use the output of fLPS?
    - For a preview, host and perturber files can be read in Python with `pandas.read_csv(path, sep=r'\s+', header=None, index_col=False, names=['time[Myr]', 'star_index', 'star_name', 'mass[solar_mass]', 'x[pc]', 'y[pc]', 'z[pc]', 'vx[kmps]', 'vy[kmps]', 'vz[kmps]'])`  (Watch out! reading the entire file may take a lot of memory!)
    - Note that in NBODY6++GPU, star_name is the only identifier of a star, and star_index should hardly be used.
    - Each line in the host file corresponds to 10 lines in the perturber file, which are the 10 nearest neighbor stars of the host star. For example, the first line (line 0) of host file corresponds to the first 10 lines (line 0 - 9) of perturber file; the second line (line 1) of host file corresponds to the second 10 lines (line 10 - 19) of perturber file; the line X of host file corresponds to the 10 lines (line 10*X - 10*X+9) of perturber file. 
    - Besides host and perturber files, LPSdiag.txt save some diagnostics especially when the host star escapes from the star cluster.

### Run creb 

For one planetary system simulation perturbed by stellar flybys, use

`./creb --n6-result-dir=/path/to/n6/result --starname=starname --nuseptb=nuseptb --input-file-path=/path/to/input/file --reb-output-subdir=subdirname`

The creb file cannot be paralleled so it uses one CPU thread (because the OPENMP in REBOUND has some restrictions). If you would like to fully utilize all CPU cores, you are recommended to run multiple instances at the same time, using the monitor script `tools/rebound_c_monitor.py`. See usage using `python3 rebound_c_monitor.py --help` .

# FAQ

Oh I have not gotten any FAQ yet :P . You are welcomed to ask any question in our discussion, which will go to my email automatically ~

# List of publications that ultilised *LPS+*

- Wu et al. 2023  https://ui.adsabs.harvard.edu/abs/2023MNRAS.523.4801W/abstract (Currently I do not a "code paper" of LPS+, and I will consider submitting one soon. If you use this code, you may cite this <-- paper, where this code is for the first time. )
- Wu et al. 2024, submitted to MNRAS

The following papers ultilised the predecessor code LonelyPlanetS 
- Cai et al. 2017  https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.4337C
- Cai et al. 2019  https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4311C
- Flammini Dotti et al. 2019  https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.2280F
- Veras et al. 2020  https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.5062V
- Stock et al. 2020  https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.1807S
- Stock et al. 2022  https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.2460S
- Benkendorff et al. 2024 https://ui.adsabs.harvard.edu/abs/2024MNRAS.528.2834B


# TODO

- [ ] Add support to the latest REBOUND API of >=4.0


# License
[CC-BY-SA-4.0 license](https://github.com/kaiwu-astro/lps_plus/blob/main/LICENSE).
