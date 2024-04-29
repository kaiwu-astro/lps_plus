#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <omp.h>
#include <inttypes.h>
#include "rebound.h"
#include "creb.h"

int use_thread = 1; // for OpenMP
int nplanet=0, nuseptb=-1, nkuiper=0; // from input
int NPLANET0, NUSEPTB0, NKUIPER0; // initial value
int particle_i_escape[N_PARTICLE_MAX] = {0}; // 0: not escaped; 1: escaped
double ttot;               // ttot is in input, or read from fLPS.f output
double resetTime = 1e-4;   // reset-time: small integration time after reset. in Myr. The smaller, the more accurate
double output_interval = 1e-3; // in Myr. Need to be >= resetTime. if <, will be forced to =resetTime
int output_N = 1000;           // specify output_interval or output_N. If both, use output_interval
#define BUFFER_SIZE 1048576 // 1048576 means 1MB

const double MAX_DE_STEP = 1.0e-5, MAX_DE_TOTAL = 1.0e-4;
const bool DO_INTERPOL_PTB = True; // if True, will interpolate ptb position and velocity. If False, will use the previous one
const bool ALLOW_large_de_cannot_solve = True; // if False, will exit when energe error is too large, and cannot be solved by decreasing integration time step
const bool USE_MERCURY_INTEGRATOR = False; // will force to use mercuries integrator
const bool USE_SYMPLECTIC_INTEGRATOR = False; // Usually: False for IAS15. Switch to True for WHFAST
const double symp_dt_division_factor = 10.0;
const bool OUTPUT_TAILLINE_ONLY = False; // Should: False. Turn to True to debug if disk IO takes too much time. 
const bool PLANETS_as_MLP = False; // Should: False: Include planet-planet force.   True: Planets feels force only.
const bool KUIPER_as_MLP = True;   // Should: True. DEBUG only for False.
const bool CHECK_PTB_OFFSET = False; // if True, will check ptb offset. This is for debug
const bool DO_OUTPUT = True; // if False, will not output. This is for debug
const bool OUTPUT_PTBS = True; // if False, will not output perturbers. This is for debug

// data read from NBODY6 output
double ptb_masses[N_PTB_MAX], ptb_xs[3][N_PTB_MAX];
int nhost, n_ptb, ptb_indexes[100], planet_index_offset;
double Mscale_Mdot, Rscale_pc, Tscale_Myr, Vscale_kmps, Mscale_me, Rscale_au, Tscale_yr; // scale factors are read from LPSdiag.txt. Approximate value: Rscale_pc = 1, Tscale_Myr = 0.2
char filename[MAX_PATH_LENGTH], starname[128], logprefix[128], n6_result_dir[MAX_PATH_LENGTH], reb_output_subdir[MAX_PATH_LENGTH], i_str[5];
FILE *timefile, *infofile, *semiMajorAxisfile, *hostfile, *ptbfile;

double hosttime[NROW_LPS_HOST_MAX];
struct nbparticle host1, host2, ptb1[N_PTB_MAX], ptb2[N_PTB_MAX];
int linen = 0;

// store planetary system input
int numInputParticles = 0;
inputparticle inputParticles[N_PARTICLE_MAX];


int is_particle_escape(const struct reb_particle *particle) {
    const int ESCAPE_DISTANCE = 10000; // AU
    if (sqrt(particle->x * particle->x + particle->y * particle->y + particle->z * particle->z) > ESCAPE_DISTANCE) {
        printf("[ESCAPE] Particle hash=%" PRIu32 " escaped due to r > %d au. Coordinates: %.10e, %.10e, %.10e\n", particle->hash, ESCAPE_DISTANCE, particle->x, particle->y, particle->z);
        return 1;
    }
    return 0;
}

int output(struct reb_simulation *const r, const double simutime_earth, FILE *f) {
    if(!DO_OUTPUT) return 0;
    struct reb_orbit orbit;
    struct reb_particle *_p;
    double total = 0.0;

    if (OUTPUT_PTBS) {
        for (int iptb = 1; iptb < nuseptb + 1; iptb++) {
            orbit = reb_tools_particle_to_orbit(r->G, r->particles[iptb], r->particles[0]);
            fprintf(f, "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
                    simutime_earth / myr_to_earth_system_time, orbit.a, orbit.e, orbit.inc, orbit.omega, orbit.Omega, orbit.f,
                    r->particles[iptb].m,
                    (r->particles[iptb].x - r->particles[0].x), (r->particles[iptb].y - r->particles[0].y), (r->particles[iptb].z - r->particles[0].z),
                    (r->particles[iptb].vx - r->particles[0].vx) / kmps_to_earth_system_speed, (r->particles[iptb].vy - r->particles[0].vy) / kmps_to_earth_system_speed, (r->particles[iptb].vz - r->particles[0].vz) / kmps_to_earth_system_speed,
                    0.0);
            total += simutime_earth / myr_to_earth_system_time + orbit.a + orbit.e + orbit.inc + orbit.omega + orbit.Omega + orbit.f + r->particles[iptb].m + (r->particles[iptb].x - r->particles[0].x) + (r->particles[iptb].y - r->particles[0].y) + (r->particles[iptb].z - r->particles[0].z) + (r->particles[iptb].vx - r->particles[0].vx) / kmps_to_earth_system_speed + (r->particles[iptb].vy - r->particles[0].vy) / kmps_to_earth_system_speed + (r->particles[iptb].vz - r->particles[0].vz) / kmps_to_earth_system_speed;
        }
    }

    if (CHECK_PTB_OFFSET) return 0; 

    for (int name_pl = 0; name_pl < NPLANET0; name_pl++) {
        sprintf(i_str, "%d", name_pl);
        // printf("name_pl: %d, i_str: %s\n", name_pl, i_str);
        _p = reb_get_particle_by_hash(r, reb_hash(i_str));
        if (_p) {
            orbit = reb_tools_particle_to_orbit(r->G, *_p, r->particles[0]);
            fprintf(f, "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
                    simutime_earth / myr_to_earth_system_time, orbit.a, orbit.e, orbit.inc, orbit.omega, orbit.Omega, orbit.f,
                    _p->m,
                    (_p->x - r->particles[0].x), (_p->y - r->particles[0].y), (_p->z - r->particles[0].z),
                    (_p->vx - r->particles[0].vx) / kmps_to_earth_system_speed, (_p->vy - r->particles[0].vy) / kmps_to_earth_system_speed, (_p->vz - r->particles[0].vz) / kmps_to_earth_system_speed,
                    orbit.P / yr_to_earth_system_time);
            total -= simutime_earth / myr_to_earth_system_time + orbit.a + orbit.e + orbit.inc + orbit.omega + orbit.Omega + orbit.f + _p->m + (_p->x - r->particles[0].x) + (_p->y - r->particles[0].y) + (_p->z - r->particles[0].z) + (_p->vx - r->particles[0].vx) / kmps_to_earth_system_speed + (_p->vy - r->particles[0].vy) / kmps_to_earth_system_speed + (_p->vz - r->particles[0].vz) / kmps_to_earth_system_speed + orbit.P / yr_to_earth_system_time;
        }
        else {
            // escaped. make every output number to be 6, but keep the format
            fprintf(f, "%.10e\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\n",
                    simutime_earth / myr_to_earth_system_time);
        }
    }

    /*
    column names:
    ['t', 'a', 'e', 'i', 'peri', 'node', 'true_anomaly', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    */
    for (int name_mlp = NPLANET0; name_mlp < NKUIPER0 + NPLANET0; name_mlp++) {
        sprintf(i_str, "%d", name_mlp);
        // printf("name_mlp: %d, i_str: %s\n", name_mlp, i_str);
        _p = reb_get_particle_by_hash(r, reb_hash(i_str));
        if (_p) {
            orbit = reb_tools_particle_to_orbit(r->G, *_p, r->particles[0]);
            fprintf(f, "%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
                    simutime_earth / myr_to_earth_system_time, orbit.a, orbit.e, orbit.inc, orbit.omega, orbit.Omega, orbit.f,
                    _p->m,
                    (_p->x - r->particles[0].x), (_p->y - r->particles[0].y), (_p->z - r->particles[0].z),
                    (_p->vx - r->particles[0].vx) / kmps_to_earth_system_speed, (_p->vy - r->particles[0].vy) / kmps_to_earth_system_speed, (_p->vz - r->particles[0].vz) / kmps_to_earth_system_speed,
                    orbit.P / yr_to_earth_system_time);
            total -= simutime_earth / myr_to_earth_system_time + orbit.a + orbit.e + orbit.inc + orbit.omega + orbit.Omega + orbit.f + _p->m + (_p->x - r->particles[0].x) + (_p->y - r->particles[0].y) + (_p->z - r->particles[0].z) + (_p->vx - r->particles[0].vx) / kmps_to_earth_system_speed + (_p->vy - r->particles[0].vy) / kmps_to_earth_system_speed + (_p->vz - r->particles[0].vz) / kmps_to_earth_system_speed + orbit.P / yr_to_earth_system_time;
        }
        else {
            // escaped. make every output number to be 6, but keep the format
            fprintf(f, "%.10e\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\t6\n",
                    simutime_earth / myr_to_earth_system_time);
        }
    }
    if (isnan(total) || isinf(total)) return 1;
    return 0;
}

void read_2_lines() {
    char buf[1024] = "";

    if (fgets(buf, sizeof(buf), hostfile) != NULL) {
        sscanf(buf, "%lf %ld %ld %lf %lf %lf %lf %lf %lf %lf",
                &host1.time, &host1.index, &host1.name, &host1.mass,
                &host1.x1, &host1.x2, &host1.x3,
                &host1.v1, &host1.v2, &host1.v3);
        host1.time *= myr_to_earth_system_time;
        host1.x1 *= parsec_to_au;
        host1.x2 *= parsec_to_au;
        host1.x3 *= parsec_to_au;
        host1.v1 *= kmps_to_earth_system_speed;
        host1.v2 *= kmps_to_earth_system_speed;
        host1.v3 *= kmps_to_earth_system_speed;
        host1.linen = linen;
        linen++;
    }
    if (fgets(buf, sizeof(buf), hostfile) != NULL) {
        sscanf(buf, "%lf %ld %ld %lf %lf %lf %lf %lf %lf %lf",
                &host2.time, &host2.index, &host2.name, &host2.mass,
                &host2.x1, &host2.x2, &host2.x3,
                &host2.v1, &host2.v2, &host2.v3);
        host2.time *= myr_to_earth_system_time;
        host2.x1 *= parsec_to_au;
        host2.x2 *= parsec_to_au;
        host2.x3 *= parsec_to_au;
        host2.v1 *= kmps_to_earth_system_speed;
        host2.v2 *= kmps_to_earth_system_speed;
        host2.v3 *= kmps_to_earth_system_speed;
        host2.linen = linen;
        linen++;
    }

    for (int i = 0; i < nuseptb; i++) {
        if (fgets(buf, sizeof(buf), ptbfile) != NULL) {
            sscanf(buf, "%lf %ld %ld %lf %lf %lf %lf %lf %lf %lf",
                    &ptb1[i].time, &ptb1[i].index, &ptb1[i].name, &ptb1[i].mass,
                    &ptb1[i].x1, &ptb1[i].x2, &ptb1[i].x3,
                    &ptb1[i].v1, &ptb1[i].v2, &ptb1[i].v3);
            ptb1[i].time *= myr_to_earth_system_time;
            ptb1[i].x1 *= parsec_to_au;
            ptb1[i].x2 *= parsec_to_au;
            ptb1[i].x3 *= parsec_to_au;
            ptb1[i].v1 *= kmps_to_earth_system_speed;
            ptb1[i].v2 *= kmps_to_earth_system_speed;
            ptb1[i].v3 *= kmps_to_earth_system_speed;
        }
    }
    for (int i = 0; i < nuseptb; i++) {
        if (fgets(buf, sizeof(buf), ptbfile) != NULL) {
            sscanf(buf, "%lf %ld %ld %lf %lf %lf %lf %lf %lf %lf",
                    &ptb2[i].time, &ptb2[i].index, &ptb2[i].name, &ptb2[i].mass,
                    &ptb2[i].x1, &ptb2[i].x2, &ptb2[i].x3,
                    &ptb2[i].v1, &ptb2[i].v2, &ptb2[i].v3);
            ptb2[i].time *= myr_to_earth_system_time;
            ptb2[i].x1 *= parsec_to_au;
            ptb2[i].x2 *= parsec_to_au;
            ptb2[i].x3 *= parsec_to_au;
            ptb2[i].v1 *= kmps_to_earth_system_speed;
            ptb2[i].v2 *= kmps_to_earth_system_speed;
            ptb2[i].v3 *= kmps_to_earth_system_speed;
        }
    }
}

void read_one_line() {
    // 移动一行，然后读入一行
    char buf[1024] = "";

    host1 = host2;
    if (fgets(buf, sizeof(buf), hostfile) != NULL) {
        sscanf(buf, "%lf %ld %ld %lf %lf %lf %lf %lf %lf %lf",
                &host2.time, &host2.index, &host2.name, &host2.mass,
                &host2.x1, &host2.x2, &host2.x3,
                &host2.v1, &host2.v2, &host2.v3);
        host2.time *= myr_to_earth_system_time;
        host2.x1 *= parsec_to_au;
        host2.x2 *= parsec_to_au;
        host2.x3 *= parsec_to_au;
        host2.v1 *= kmps_to_earth_system_speed;
        host2.v2 *= kmps_to_earth_system_speed;
        host2.v3 *= kmps_to_earth_system_speed;
        host2.linen = linen;
        linen++;
    }

    for (int i = 0; i < nuseptb; i++) {
        ptb1[i] = ptb2[i];
        if (fgets(buf, sizeof(buf), ptbfile) != NULL) {
            sscanf(buf, "%lf %ld %ld %lf %lf %lf %lf %lf %lf %lf",
                    &ptb2[i].time, &ptb2[i].index, &ptb2[i].name, &ptb2[i].mass,
                    &ptb2[i].x1, &ptb2[i].x2, &ptb2[i].x3,
                    &ptb2[i].v1, &ptb2[i].v2, &ptb2[i].v3);
            ptb2[i].time *= myr_to_earth_system_time;
            ptb2[i].x1 *= parsec_to_au;
            ptb2[i].x2 *= parsec_to_au;
            ptb2[i].x3 *= parsec_to_au;
            ptb2[i].v1 *= kmps_to_earth_system_speed;
            ptb2[i].v2 *= kmps_to_earth_system_speed;
            ptb2[i].v3 *= kmps_to_earth_system_speed;
        }
    }
}

struct reb_simulation *construct_tmpsim(FILE *outputfile) {
    // tmpsim init: save initial info of simulation
    // double mean_anomaly, True_anomaly, e, mlp_a_list[N_PARTICLE_MAX];
    struct reb_simulation *tmpsim = reb_create_simulation();
    struct reb_particle host = {0}, planetary;
    // tmpsim HOST
    read_2_lines();
    host.m = host1.mass;
    reb_add(tmpsim, host);
    // tmpsim ptb.
    // this is unnecessary, but to make the output function consistent (for the first output), and the whole code structure consistent, I add some fake ptbs here
    for (int ptbind = 0; ptbind < nuseptb; ptbind++)
        reb_add(tmpsim, reb_tools_orbit_to_particle(tmpsim->G, host, 0.2333, 100 + ptbind * 1, 0.0, 0.0, 0.0, 0.0, 0.0));

    // add particle from inputParticles
    // printf("numInputParticles: %d\n", numInputParticles);
    for (int i = 0; i < numInputParticles; i++)
    {
        int mean_anomaly_rad = inputParticles[i].M_deg * PI / 180.0;
        int true_anomaly_rad = reb_tools_M_to_f(inputParticles[i].e, mean_anomaly_rad);
        planetary = reb_tools_orbit_to_particle(tmpsim->G, host, 
            inputParticles[i].m_me / Mdot_to_Me, 
            inputParticles[i].a_au, 
            inputParticles[i].e, 
            inputParticles[i].i_deg * PI / 180.0, 
            inputParticles[i].node_deg * PI / 180.0, 
            inputParticles[i].peri_deg * PI / 180.0, 
            true_anomaly_rad);
        sprintf(i_str, "%d", i);
        planetary.hash = reb_hash(i_str);
        reb_add(tmpsim, planetary);
        // print planetary m, x, y, z, vx, vy, vz
        // printf("m: %.10e, x: %.10e, y: %.10e, z: %.10e, vx: %.10e, vy: %.10e, vz: %.10e, hash: %" PRIu32 "\n",
        //         planetary.m,
        //         planetary.x, planetary.y, planetary.z,
        //         planetary.vx, planetary.vy, planetary.vz, planetary.hash);
    }
    return tmpsim;
}

int last_indbefore = 0;
int *get_interval(const double *time) {
    /* return a[0] = indBefore, and a[1] = indAfter.
    */
    static int a[2];
    int last_indafter;
    last_indafter = last_indbefore + 1;
    if (hosttime[last_indafter] > *time) { // still in the last interval
        a[0] = last_indbefore; 
        a[1] = last_indafter;
    }
    else {
        a[1] = last_indafter; // find next indAfter starting from last_indafter
        while (hosttime[a[1]] < *time) a[1]++;
        a[0] = a[1] - 1;
    }

    last_indbefore = a[0];
    return a;
}

double get_fraction(const double *t1, const double *t2, const double *tin) {
    return 1.0 - (*tin - *t1) / (*t2 - *t1);
}

struct nbparticle *nb_interpol_host(const int *indBefore, const int *indAfter, const double *time, const double *frac1){
    static struct nbparticle hostp;
    double frac2 = 1.0 - *frac1;

    hostp.time = *time;
    hostp.index = host1.index;
    hostp.name = host1.name;
    hostp.mass = host1.mass * *frac1 + host2.mass * frac2;
    hostp.x1 = host1.x1 * *frac1 + host2.x1 * frac2;
    hostp.x2 = host1.x2 * *frac1 + host2.x2 * frac2;
    hostp.x3 = host1.x3 * *frac1 + host2.x3 * frac2;
    hostp.v1 = host1.v1 * *frac1 + host2.v1 * frac2;
    hostp.v2 = host1.v2 * *frac1 + host2.v2 * frac2;
    hostp.v3 = host1.v3 * *frac1 + host2.v3 * frac2;
    return &hostp;
}

struct nbparticle *nb_interpol_ptb(const int *indBefore_host, const int *indAfter_host, const double *time, const double *frac1, const int *ptbid){
    static struct nbparticle ptbp;
    double frac2 = 1.0 - *frac1;
    // int indBefore_ptb = *indBefore_host * n_ptb + *ptbid;
    // int indAfter_ptb = *indAfter_host * n_ptb + *ptbid;

    if (DO_INTERPOL_PTB == True && ptb1[*ptbid].index == ptb2[*ptbid].index && ptb1[*ptbid].name == ptb2[*ptbid].name)
    {
        ptbp.time = *time;
        ptbp.index = ptb1[*ptbid].index;
        ptbp.name = ptb1[*ptbid].name;
        ptbp.mass = ptb1[*ptbid].mass * *frac1 + ptb2[*ptbid].mass * frac2;
        ptbp.x1 = ptb1[*ptbid].x1 * *frac1 + ptb2[*ptbid].x1 * frac2;
        ptbp.x2 = ptb1[*ptbid].x2 * *frac1 + ptb2[*ptbid].x2 * frac2;
        ptbp.x3 = ptb1[*ptbid].x3 * *frac1 + ptb2[*ptbid].x3 * frac2;
        ptbp.v1 = ptb1[*ptbid].v1 * *frac1 + ptb2[*ptbid].v1 * frac2;
        ptbp.v2 = ptb1[*ptbid].v2 * *frac1 + ptb2[*ptbid].v2 * frac2;
        ptbp.v3 = ptb1[*ptbid].v3 * *frac1 + ptb2[*ptbid].v3 * frac2;
    }
    else { // [ptb changed - use the previous one] OR [ptb interpol disabled]
        ptbp.time = *time;
        ptbp.index = ptb1[*ptbid].index;
        ptbp.name = ptb1[*ptbid].name;
        ptbp.mass = ptb1[*ptbid].mass;
        ptbp.x1 = ptb1[*ptbid].x1;
        ptbp.x2 = ptb1[*ptbid].x2;
        ptbp.x3 = ptb1[*ptbid].x3;
        ptbp.v1 = ptb1[*ptbid].v1;
        ptbp.v2 = ptb1[*ptbid].v2;
        ptbp.v3 = ptb1[*ptbid].v3;
    }

    return &ptbp;
}

void refresh_host_ptb(int line1, int line2) {
    if (line1 == host1.linen && line2 == host2.linen)
        return;
    else if (line1 == host2.linen) {
        read_one_line();
        return;
    }
    else if (line1 > host2.linen) {
        read_2_lines();
        return;
    }
    else {
        printf("ERROR when reading new line.\nThey are: infb4: %d, indafter: %d, host1.linen: %ld, host2.linen: %ld\n", line1, line2, host1.linen, host2.linen);
        exit(1);
    }
}

int add_particles_from_NB6data(struct reb_simulation *sim, const struct reb_particle *prevhost, const struct reb_particle prevps[], const struct reb_particle prevmlps[], double ctp) { // current time (pointer)
    int *b4aft, indBefore, indAfter;
    b4aft = get_interval(&ctp); indBefore = b4aft[0]; indAfter = b4aft[1];
    //决定是否要重新读文件
    refresh_host_ptb(indBefore, indAfter);
    struct nbparticle *hostp, *ptbp;
    double frac1 = get_fraction(&hosttime[indBefore], &hosttime[indAfter], &ctp);
    if (frac1 > 1 || frac1 < 0) {
        printf("%s\tERROR: frac1 = %lf ", logprefix, frac1);
        exit(1);
    }
    hostp = nb_interpol_host(&indBefore, &indAfter, &ctp, &frac1);
    int iptb, iplanet, imlp;
    // unsigned int hash_int;
    struct reb_particle planet, ptb, mlp;
    double total_value = 0.0;
    for (iptb = 0; iptb < nuseptb; iptb++) {
        ptbp = nb_interpol_ptb(&indBefore, &indAfter, &ctp, &frac1, &iptb);
        ptb.m = ptbp->mass;
        ptb.x = ptbp->x1 - hostp->x1;
        ptb.y = ptbp->x2 - hostp->x2;
        ptb.z = ptbp->x3 - hostp->x3;
        ptb.vx = ptbp->v1 - hostp->v1;
        ptb.vy = ptbp->v2 - hostp->v2;
        ptb.vz = ptbp->v3 - hostp->v3;
        int safe = True;
        for (int addPtbInd = 1; addPtbInd < iptb + 1; addPtbInd++) {// start from 1: skip central particle
            if (abs(sim->particles[addPtbInd].x - ptb.x) < 0.01) {
                safe = False; 
                break;
            }
            // distance < 0.01 AU for ptb is too close. This must be the same particle.
            // This may happen when NBList of current host in Nbody6 is too short (< 10), 
            //    so the rest position of NBList is filled with neighbors of previous host
            //    which can be regarded as random to current host
            // To solve this problem at smallest cost of modifying code, the solution here is to
            //    put this particle 100 times further, which actually remove the influence of this particle
        }
        if (!safe) {
            // + iptb * 100 is to solve the case: same starid appear in NBList many times
            // in the extreme case: 
            //    NBLIST: [1, 702, 702, 702, 702, 702, 702 ...]
            ptb.x  *= (1000.0 + iptb * 100);
            ptb.y  *= (1000.0 + iptb * 100);
            ptb.z  *= (1000.0 + iptb * 100);
        }
        reb_add(sim, ptb);
        total_value += ptb.m + ptb.x + ptb.y + ptb.z + ptb.vx + ptb.vy + ptb.vz;
    }

    int _nescaped = 0;

    for (iplanet = 0; iplanet < nplanet; iplanet++) {
        planet.m = prevps[iplanet].m;
        planet.x = prevps[iplanet].x - prevhost->x;
        planet.y = prevps[iplanet].y - prevhost->y;
        planet.z = prevps[iplanet].z - prevhost->z;
        planet.vx = prevps[iplanet].vx - prevhost->vx;
        planet.vy = prevps[iplanet].vy - prevhost->vy;
        planet.vz = prevps[iplanet].vz - prevhost->vz;
        planet.hash = prevps[iplanet].hash;
        // printf("hash3: %" PRIu32 "\n", planet.hash);
        if (is_particle_escape(&planet) == 0) {
            reb_add(sim, planet);
            total_value += planet.m + planet.x + planet.y + planet.z + planet.vx + planet.vy + planet.vz;
        }
        else {
            _nescaped += 1;
        }
    }
    nplanet -= _nescaped;

    _nescaped = 0;
    for (imlp = 0; imlp < nkuiper; imlp++) {
        mlp.m = prevmlps[imlp].m;
        mlp.x = prevmlps[imlp].x - prevhost->x;
        mlp.y = prevmlps[imlp].y - prevhost->y;
        mlp.z = prevmlps[imlp].z - prevhost->z;
        mlp.vx = prevmlps[imlp].vx - prevhost->vx;
        mlp.vy = prevmlps[imlp].vy - prevhost->vy;
        mlp.vz = prevmlps[imlp].vz - prevhost->vz;
        mlp.hash = prevmlps[imlp].hash;
        // printf("hash4: %" PRIu32 "\n", mlp.hash);
        if (is_particle_escape(&mlp) == 0) {
            reb_add(sim, mlp);
            total_value += mlp.m + mlp.x + mlp.y + mlp.z + mlp.vx + mlp.vy + mlp.vz;
        }
        else {
            // hash_int = strtoimax(mlp.hash, NULL, 10);
            // particle_i_escape[hash_int] = True;
            _nescaped += 1;
        }
    }
    nkuiper -= _nescaped;

    if (isnan(total_value) || isinf(total_value)){
        printf("problem.c abort because of nan/inf value detected when adding particles\n");
        printf("Check NBODY6++GPU output around %.10e\n", ctp/myr_to_earth_system_time);
        return 1;
    }
    return 0;
}

void start_simulation() {
    clock_t starttime, endtime, simustart, simuclock = 0;
    starttime = clock();
    last_indbefore = 0;
    int iplanet, imlp;
    char printbuf[1024];
    struct reb_particle centralhost, prevhost, prevps[8], prevmlps[N_PARTICLE_MAX];

    sprintf(filename, "%s/reb_out_%s.txt", reb_output_subdir, starname);
    // 如果后面要修改输出文件名，请注意不要在文件名中加入dot(.)。也就是，不要是用reb.out.123.txt这样的文件名。这是因为后面数据分析时，为了压缩输出文件，会采用.pkl.zst的后缀，而wukai现有的数据分析脚本，会把文件名用python:os.path.basename().split('.')[0]来获取，因此多于的点会导致文件名被错误识别。
    // If the output filename is to be changed later, be careful not to add dot(.) to the filename. That is, do not use reb.out.123.txt. This is because when you analyze the data later, in order to compress the output file, you will use the .pkl.zst suffix, and wukai's existing data analysis script, you will use python:os.path.basename().split('.') for the file name. [0] to get it, so more dots than that will cause the filename to be misrecognized.
    FILE *f = fopen(filename, "w+");
    if (f == NULL) {
        printf("Error: Could not open file %s\n", filename);
        return 1; 
    }
    setvbuf(f, NULL, _IOFBF, BUFFER_SIZE);

    sprintf(filename, "%s/reb_out_%s_EnergyError.txt", reb_output_subdir, starname);
    FILE *de_file = fopen(filename, "w+");
    setvbuf(de_file, NULL, _IOFBF, BUFFER_SIZE);

    struct reb_simulation *tmpsim = construct_tmpsim(f);

    // prev particles array init
    centralhost = tmpsim->particles[0]; 
    prevhost = tmpsim->particles[0];
    for (iplanet = 0; iplanet < nplanet; iplanet++)
        prevps[iplanet] = tmpsim->particles[iplanet + planet_index_offset];
    for (imlp = 0; imlp < nkuiper; imlp++)
        prevmlps[imlp] = tmpsim->particles[imlp + nplanet + planet_index_offset];

    // simulation: loop over dataframe
    struct reb_simulation *sim;
    double energy_init, energy_end, energy_error, energy_error_total = 0.0, energy_error_abs_total = 0.0;
    double simutime_earth = 0.0, next_output_time;
    int outputbad = False, nanfound = False, large_de_cannot_solve = False, n_de_too_large = 0, retry = False; 
    sim = reb_create_simulation(); // to be freed by the first line
    next_output_time = 0.0;
    while(simutime_earth < ttot) {
        reb_free_simulation(sim);
        sim = reb_create_simulation();
        reb_add(sim, centralhost);
        // print nuseptb, nplanet, nkuiper
        // printf("simutime_earth: %lf, nuseptb: %d, nplanet: %d, nkuiper: %d\n", simutime_earth/myr_to_earth_system_time, nuseptb, nplanet, nkuiper);
        // // if nplanet is 0 and nkuiper is 0, then exit the program
        // if (nplanet == 0 && nkuiper == 0) {
        //     printf("nplanet and nkuiper are both 0. Exit.\n");
        //     exit(1);
        // }
        if (add_particles_from_NB6data(sim, &prevhost, prevps, prevmlps, simutime_earth) == 1){
            nanfound = 1;
            break;
        }

        sim->testparticle_type = 0; // only feel force. no give force 
        sim->N_active = 1 + nuseptb;

        // boundary and collisions. See https://rebound.readthedocs.io/en/latest/collisions/
        // radius of each particle is needed for collision.
        // sim->collision = REB_COLLISION_LINE;
        // sim->collision = REB_COLLISION_TREE;
        // reb_configure_box(sim, 12000, 1, 1, 1);
        // sim->collision_resolve = reb_collision_resolve_merge;

        if (!PLANETS_as_MLP) sim->N_active += nplanet;
        if (!KUIPER_as_MLP) sim->N_active += nkuiper;

        if (USE_SYMPLECTIC_INTEGRATOR == False) {
            sim->integrator = REB_INTEGRATOR_IAS15;
            // sim->dt = 10 / 1.0e6 * myr_to_earth_system_time;
            if (!retry) sim->ri_ias15.epsilon = 1.0e-7;
            else sim->ri_ias15.epsilon = 1.0e-7 / pow(2.0, (double)n_de_too_large);
            sim->exact_finish_time = 0; // use 0 to accelerate
            // struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
            // sim->ri_ias15.min_dt = min(orbit.P / 30 , resetTime);
        }
        else {
            sim->integrator = REB_INTEGRATOR_WHFAST;
            // sim->integrator = REB_INTEGRATOR_LEAPFROG;
            // sim->collision = REB_COLLISION_TREE;
            // sim->collision_resolve = reb_collision_resolve_hardsphere;
            // sim->boundary = REB_BOUNDARY_OPEN;
            // sim->softening = 0.01;
            sim->dt = 90.0 / symp_dt_division_factor;
            sim->exact_finish_time = 0; // 0 is necessary for symplectic integrator
        }
        if (USE_MERCURY_INTEGRATOR == True) {
            sim->integrator = REB_INTEGRATOR_MERCURIUS;
        }

        energy_init = reb_tools_energy(sim);

        // output info at t=0
        if (simutime_earth == 0.0 || CHECK_PTB_OFFSET)
            output(sim, simutime_earth, f);

        simustart = clock();
        if (!retry)
            reb_integrate(sim, resetTime);
        else
            reb_integrate(sim, resetTime / (double)(n_de_too_large+1));
        simuclock += clock() - simustart;

        simutime_earth += sim->t;

        energy_end = reb_tools_energy(sim);
        energy_error = fabs((energy_init - energy_end) / energy_init);
        if (isnan(energy_error)) energy_error = 0.0;
        if (abs(energy_error) > MAX_DE_STEP)
            n_de_too_large++;
        else {
            n_de_too_large = 0;
            retry = False;
        }

        if (n_de_too_large > 0) {
            if (n_de_too_large <= 4) {
                // retry this reset-step will higher integration precision, and smaller resetTime
                retry = True;
                printf("de too large at %.10e is %.10e. Retry... %d times\n", simutime_earth / myr_to_earth_system_time, energy_error, n_de_too_large);
                simutime_earth -= sim->t;
                continue;
            }
            else
                large_de_cannot_solve = True;
        }


        fprintf(de_file, "%.10e\t %.10e\n", simutime_earth / myr_to_earth_system_time, energy_error);
        energy_error_total += energy_error;
        energy_error_abs_total += abs(energy_error);

        outputbad = False;

        if (simutime_earth >= next_output_time && !OUTPUT_TAILLINE_ONLY){
            outputbad = output(sim, simutime_earth, f);
            next_output_time += output_interval;
        }
        else if (OUTPUT_TAILLINE_ONLY && simutime_earth > ttot)
            outputbad = output(sim, simutime_earth, f);

        if (large_de_cannot_solve && ALLOW_large_de_cannot_solve) {
            printf("energy error too large during simulation and cannot solve, but is allowed, at %.10e Myr\n", simutime_earth/myr_to_earth_system_time);
            n_de_too_large = 0;
            retry = False;
            large_de_cannot_solve = False;
        }

        if (outputbad || (large_de_cannot_solve && !ALLOW_large_de_cannot_solve)) {
            if (outputbad) {
                printf("nan/inf occured during integration. \n");
                nanfound = 1;
            }
            if (large_de_cannot_solve)
                printf("energy error too large during simulation. \n");
            printf("Data for re-produce this error:\n");
            struct reb_simulation *tsim = reb_create_simulation();
            reb_add(tsim, centralhost);
            add_particles_from_NB6data(tsim, &prevhost, prevps, prevmlps, simutime_earth - sim->t);
            output(tsim, simutime_earth - sim->t, stdout);
            reb_free_simulation(tsim);
            break;
        }

        prevhost = sim->particles[0];
        for (iplanet = 0; iplanet < nplanet; iplanet++)
            prevps[iplanet] = sim->particles[iplanet + planet_index_offset];
        for (imlp = 0; imlp < nkuiper; imlp++)
            prevmlps[imlp] = sim->particles[imlp + nplanet + planet_index_offset];
        // printf("\rdt = %.10e",sim->dt);

        // reb_free_simulation(sim); // moved to the beginning of the program, in order to keep the last sim for log
        // a sep-simulation finishes here. Then ptbs will be reset to their position
        // for nuseptb = 0, resetTime = ttot, there is only 1 sep-simulation
    } 

    reb_free_simulation(tmpsim);
    fclose(f);

    fprintf(de_file, "Arithmetic sum of relative energy error: %.10e\n", energy_error_total);
    fprintf(de_file, "Arithmetic sum of absolute value of relative energy error: %.10e\n", energy_error_abs_total);
    fclose(de_file);
    // record info (to be read by plot)
    sprintf(filename, "%s/reb_time.txt", reb_output_subdir);
    timefile = fopen(filename, "a");
    sprintf(filename, "%s/reb_info.txt", reb_output_subdir);
    infofile = fopen(filename, "a");

    sprintf(printbuf, "%s %d, %d, %d\n", logprefix, nuseptb, nplanet, nkuiper);
    fprintf(infofile, "%s", printbuf);
    sprintf(printbuf, "%s Integrator: %d, last dt = %.10e year\n", logprefix, sim->integrator, sim->dt / yr_to_earth_system_time);
    fprintf(infofile, "%s", printbuf);

    // record problem
    if (energy_error_abs_total > MAX_DE_TOTAL){
        sprintf(printbuf, "%s Arithmetic sum of energy error: %.2lf%%. Sum of abs: %.2lf%%. Too large!!!", logprefix, energy_error_total / 100.0, energy_error_abs_total / 100.0);
        if (large_de_cannot_solve) strcat(printbuf, " Abort.\n");
        else strcat(printbuf, "\n");
        printf("%s", printbuf);
        fprintf(infofile, "%s", printbuf);
    }
    if (nanfound == 1){
        sprintf(printbuf, "%s abort because of nan/inf\n", logprefix);
        printf("%s", printbuf);
        fprintf(infofile, "%s", printbuf);
    }

    // record Host star mass
    sprintf(printbuf, "%s Host mass: %.15e\n", logprefix, host1.mass);
    printf("%s", printbuf);
    // printf("----------------\n");
    fprintf(infofile, "%s", printbuf);

    // record simulation time
    sprintf(printbuf, "%s Simu used: %lf seconds\n", logprefix, ((double)simuclock) / CLOCKS_PER_SEC / (double)use_thread);
    printf("%s", printbuf);
    fprintf(infofile, "%s", printbuf);
    fprintf(timefile, "%s", printbuf);

    endtime = clock();
    sprintf(printbuf, "%s Total: %lf seconds\n", logprefix, ((double)(endtime - starttime)) / CLOCKS_PER_SEC / (double)use_thread);
    printf("%s", printbuf);
    fprintf(infofile, "%s", printbuf);
    fprintf(timefile, "%s", printbuf);

    fprintf(infofile, "----------------\n");

    reb_free_simulation(sim); // free the last sim in loop
    fclose(timefile);
    fclose(infofile);
    // fclose(console_stdout);
    // stdout = fdopen(STDOUT_FILENO, "w");
    // setvbuf(stdout, NULL, _IONBF, 0);
    // fclose(console_stderr);
    // stderr = fdopen(STDERR_FILENO, "w");
    // setvbuf(stderr, NULL, _IONBF, 0);
    }

void read_scale() {
    // read scaling factors from nb6 output

    char buf[1024] = "", fname[MAX_PATH_LENGTH] = "", junk[1024] = "";
    FILE *scalefile;

    sprintf(fname, "%s%s", n6_result_dir, "/LPSdiag.txt");
    scalefile = fopen(fname, "r");
    char* _result = fgets(buf, sizeof(buf), scalefile);
    sscanf(buf, "%s %d", junk, &nhost);
    _result = fgets(buf, sizeof(buf), scalefile);
    sscanf(buf, "%s %d", junk, &n_ptb);
    _result = fgets(buf, sizeof(buf), scalefile);
    sscanf(buf, "%s %lf", junk, &Mscale_Mdot);
    _result = fgets(buf, sizeof(buf), scalefile);
    sscanf(buf, "%s %lf", junk, &Rscale_pc);
    _result = fgets(buf, sizeof(buf), scalefile);
    sscanf(buf, "%s %lf", junk, &Vscale_kmps);
    _result = fgets(buf, sizeof(buf), scalefile);
    sscanf(buf, "%s %lf", junk, &Tscale_Myr);
    fclose(scalefile);

    printf("%s:\tnhost: %d\tnptb: %d\tMscale: %lf\tRscale: %lf\tVscale: %lf\tTscale: %lf\t\n",
           starname, nhost, n_ptb, Mscale_Mdot, Rscale_pc, Vscale_kmps, Tscale_Myr);

    Mscale_me = Mscale_Mdot * Mdot_to_Me;
    Rscale_au = Rscale_pc * parsec_to_au;
    Tscale_yr = Tscale_Myr * 1.e6;
    return;
}

void read_time_params(int argc, char *argv[]) {
    // read host time
    FILE *hf;
    double maxtime = -1.0; // in NBODY
    char fname[MAX_PATH_LENGTH] = "", buf[1024] = "";
    sprintf(fname, "%s%s%s%s", n6_result_dir, "/LPS_Host_", starname, ".txt");
    hf = fopen(fname, "r");
    for (int i = 0; i < NROW_LPS_HOST_MAX; i++) {
        if (fgets(buf, sizeof(buf), hf) != NULL) {
            sscanf(buf, "%lf", &hosttime[i]);
            hosttime[i] *= myr_to_earth_system_time;
            maxtime = hosttime[i];
        }
        else break;
    }
    fclose(hf);

    switch(argc) {
        case 8:  // time given by input: ./rebound, dir, num, nuseptb, npl, nkui, time
            ttot = atof(argv[7]) * myr_to_earth_system_time;  // time input in Myr. Convert to earth_unit_system.
            printf("%s TTOT overridden by input!!!: %lf Myr\n", logprefix, ttot / myr_to_earth_system_time);
            break;
        case 4:  // time given by input: ./rebound, dir, num, time
            ttot = atof(argv[3]) * myr_to_earth_system_time;  // time input in Myr. Convert to earth_unit_system.
            printf("%s TTOT overridden by input!!!: %lf Myr\n", logprefix, ttot / myr_to_earth_system_time);
            break;
        default: // no inputed time. use maxtime in file
            ttot = maxtime; // time read is in earth-system time.
            printf("%s TTOT read from host file: %lf Myr\n", logprefix, ttot / myr_to_earth_system_time);
            break;
    }

    resetTime *= myr_to_earth_system_time;
    output_interval *= myr_to_earth_system_time;
    if (output_interval == 0)
        output_interval = ttot / (double)output_N;
    else
        output_N = 1 + (int)ceil(ttot / output_interval);
    if (output_interval < resetTime)
        output_interval = resetTime;

    return;
}

void open_host_ptb_file() {
    char hostFilePath[MAX_PATH_LENGTH] = "";
    char perturberFilePath[MAX_PATH_LENGTH] = "";
    const char *HOST_FILE_PREFIX = "/LPS_Host_";
    const char *PERTURBER_FILE_PREFIX = "/LPS_Perturber_";

    sprintf(hostFilePath, "%s%s%s%s", n6_result_dir, HOST_FILE_PREFIX, starname, ".txt");
    sprintf(perturberFilePath, "%s%s%s%s", n6_result_dir, PERTURBER_FILE_PREFIX, starname, ".txt");

    if ((hostfile = fopen(hostFilePath, "r")) == NULL)
    {
        printf("Failed to open host file: %s\n", hostFilePath);
        exit(1);
    }
    // Enable full buffering
    if (setvbuf(hostfile, NULL, _IOFBF, BUFFER_SIZE) != 0) {
        perror("Failed to set buffer");
        return 1;
    }

    if ((ptbfile = fopen(perturberFilePath, "r")) == NULL)
    {
        printf("Failed to open perturber file: %s\n", perturberFilePath);
        fclose(hostfile);
        exit(1);
    }
    if (setvbuf(ptbfile, NULL, _IOFBF, BUFFER_SIZE) != 0) {
        perror("Failed to set buffer");
        return 1;
    }
    return;
}

void close_host_ptb_file(){
    fclose(hostfile);
    fclose(ptbfile);
}

int read_input_file(const char* filepath, inputparticle particles[N_PARTICLE_MAX]) {
    FILE* file = fopen(filepath, "r");
    if (file == NULL) {
        printf("Failed to open input file: %s\n", filepath);
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int nparticle = 0;

    while (fgets(line, sizeof(line), file) != NULL) {
        if (nparticle >= N_PARTICLE_MAX) {
            printf("Maximum number of particles exceeded.\n");
            break;
        }
        // skip comment line. Comment line starts with #. Skip empty line or line contains only spaces.
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '\t' || line[0] == ' ' || line[0] == '#') {
            continue;
        }

        char* key;
        char* value;
        inputparticle particle;
    
        // key = strtok(line, "=");
        // value = strtok(NULL, " ");
        char *save_ptr1, *save_ptr2; 
        char *kvpair = strtok_r(line, " ", &save_ptr1);

        // while (key != NULL && value != NULL) {
        while (kvpair != NULL) {
            key = strtok_r(kvpair, "=", &save_ptr2);
            value = strtok_r(NULL, "=", &save_ptr2);
            // printf("kvpair: %s ", kvpair);
            // printf("key: %s, value: %s\n", key, value);
            if (strcmp(key, "m_me") == 0) {
                particle.m_me = atof(value);
            } else if (strcmp(key, "a_au") == 0) {
                particle.a_au = atof(value);
            } else if (strcmp(key, "e") == 0) {
                particle.e = atof(value);
            } else if (strcmp(key, "i_deg") == 0) {
                particle.i_deg = atof(value);
            } else if (strcmp(key, "peri_deg") == 0) {
                // if value is string "random", generate random value from 0 to 360
                if (strcmp(value, "random") == 0) {
                    particle.peri_deg = getRandom(0, 360);
                } else {
                    particle.peri_deg = atof(value);
                }
            } else if (strcmp(key, "node_deg") == 0) {
                if (strcmp(value, "random") == 0) {
                    particle.node_deg = getRandom(0, 360);
                } else {
                    particle.node_deg = atof(value);
                }
            } else if (strcmp(key, "M_deg") == 0) {
                if (strcmp(value, "random") == 0) {
                    particle.M_deg = getRandom(0, 360);
                } else {
                    particle.M_deg = atof(value);
                }
            }

            kvpair = strtok_r(NULL, " ", &save_ptr1);
        }
        particles[nparticle] = particle;
        nparticle++;
    }
    print_input_params(particles, nparticle);
    return nparticle;
}

void print_input_params(inputparticle *particles, int nparticle) {
    printf("Input particle information:\n");
    for (int i = 0; i < nparticle; i++) {
        printf("Particle %d: \t", i + 1);
        printf("m_me: %.2f \t", particles[i].m_me);
        printf("a_au: %.2f \t", particles[i].a_au);
        printf("e: %.2f \t", particles[i].e);
        printf("i_deg: %.2f \t", particles[i].i_deg);
        printf("peri_deg: %.2f \t", particles[i].peri_deg);
        printf("node_deg: %.2f \t", particles[i].node_deg);
        printf("M_deg: %.2f \t", particles[i].M_deg);
        printf("\n");
    }
}

void generate_example_input_file(const char* filepath) {
    FILE* file = fopen(filepath, "w");
    if (file == NULL) {
        printf("Failed to create example input file: %s\n", filepath);
        return;
    }

    // Write file format description
    fprintf(file, "# Example Input File\n");
    fprintf(file, "# Format: key=value, seperated by spaces\n");
    fprintf(file, "# Lines with leading '#' or spaces are comments\n");
    fprintf(file, "# Available keys: m_me, a_au, e, i_deg, peri_deg, node_deg, M_deg\n");
    fprintf(file, "# m_me: mass of the particle in earth mass\n");
    fprintf(file, "# a_au: semi-major axis of the particle in au\n");
    fprintf(file, "# e: eccentricity of the particle\n");
    fprintf(file, "# i_deg: inclination of the particle in degrees\n");
    fprintf(file, "# peri_deg: argument of pericenter of the particle in degrees\n");
    fprintf(file, "# node_deg: longitude of the ascending node of the particle in degrees\n");
    fprintf(file, "# M_deg: mean anomaly of the particle in degrees\n");
    fprintf(file, "# peri, node, and M can also be the word: random, which randomize between 0 to 360 degrees.\n\n\n");

    // Write example input lines
    fprintf(file, "# 50-au-Jupiter\n");
    fprintf(file, "m_me=317.83  a_au=50.0  e=0.04839266  i_deg=0.0  peri_deg=14.75385  node_deg=100.55615  M_deg=34.40438\n");
    fprintf(file, "# Planet earth\n");
    fprintf(file, "m_me=1.0  a_au=1.0  e=0.0167  i_deg=0.0  peri_deg=102.93735  node_deg=-11.26064  M_deg=100.46435\n");
    fprintf(file, "# 3 comets (set mass=0 to greatly boost the simulation speed)\n");
    fprintf(file, "m_me=0.0  a_au=20.0  e=0.01  i_deg=0.01  peri_deg=random  node_deg=random  M_deg=random\n");
    fprintf(file, "m_me=0.0  a_au=25.0  e=0.01  i_deg=0.01  peri_deg=random  node_deg=random  M_deg=random\n");
    fprintf(file, "m_me=0.0  a_au=100.0  e=0.01  i_deg=0.01  peri_deg=random  node_deg=random  M_deg=random\n");

    fclose(file);
}

int main(int argc, char* argv[]){

    char inputFilePath[MAX_PATH_LENGTH] = "";
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--n6-result-dir=", 16) == 0) {
            sprintf(n6_result_dir, "%s", argv[i] + 16);
        } else if (strncmp(argv[i], "--starname=", 11) == 0) {
            sprintf(starname, "%s", argv[i] + 11);
        } else if (strncmp(argv[i], "--nuseptb=", 10) == 0) {
            nuseptb = atoi(argv[i] + 10);
        } else if (strncmp(argv[i], "--input-file-path=", 18) == 0) {
            sprintf(inputFilePath, "%s", argv[i] + 18);
        } else if (strncmp(argv[i], "--reb-output-subdir=", 20) == 0) {
            sprintf(reb_output_subdir, "%s", argv[i] + 20);
        }
        // on --example-input: call generate_example_input_file() and exit
        else if (strncmp(argv[i], "--example-input", 15) == 0) {
            generate_example_input_file("example_input.txt");
            // print full path
            printf("Example input file generated: %s/example_input.txt\n", getcwd(NULL, 0));
            exit(0);
        }
    }
    
    // 禁用写入缓存not working. Works the same as &>
    // setvbuf(stdout, NULL, _IONBF, 0);
    // setvbuf(stderr, NULL, _IONBF, 0);
    // sprintf(filename, "%s/reb_out_%s_consoleStdout.txt", reb_output_subdir, starname);
    // FILE *console_stdout = freopen(filename, "w+", stdout);
    // setvbuf(console_stdout, NULL, _IOFBF, 0); 
    // sprintf(filename, "%s/reb_out_%s_consoleStderr.txt", reb_output_subdir, starname);
    // FILE *console_stderr = freopen(filename, "w+", stderr);
    // setvbuf(console_stderr, NULL, _IOFBF, 0); 

    // print all input params
    printf("n6_result_dir: %s\n", n6_result_dir);
    printf("starname: %s\n", starname);
    printf("nuseptb: %d\n", nuseptb);
    printf("inputFilePath: %s\n", inputFilePath);
    printf("reb_output_subdir: %s\n", reb_output_subdir);

    // check input. If not enough, print usage and exit
    if (argc < 6) {
        printf("Usage: %s --n6-result-dir=/path/to/n6/result --starname=starname --nuseptb=nuseptb --input-file-path=/path/to/input/file --reb-output-subdir=subdirname\n", argv[0]);
        printf("Example: %s --n6-result-dir=/path/to/n6/result --starname=89 --nuseptb=10 --input-file-path=/path/to/input.txt --reb-output-subdir=reb010200\n", argv[0]);
        // indicate that --example-input is available
        printf("If input file is not prepared, use --example-input to generate an example input file.\n");
        exit(1);
    }

    // read input file
    numInputParticles = read_input_file(inputFilePath, inputParticles);
    if (numInputParticles == 0) {
        printf("Failed to read input file.\n");
        exit(1);
    }
    // 统计输入文件中的行星和Kuiper数量。依据：行星质量不为0，Kuiper质量为0
    for (int i = 0; i < numInputParticles; i++) {
        if (fabs(inputParticles[i].m_me) > 1e-9) {
            nplanet++;
        } else {
            nkuiper++;
        }
    }

    NUSEPTB0 = nuseptb;
    NPLANET0 = nplanet;
    NKUIPER0 = nkuiper;
    printf("NUSEPTB0: %d, NPLANET0: %d, NKUIPER0: %d\n", NUSEPTB0, NPLANET0, NKUIPER0);

    sprintf(logprefix, "%6s:", starname);
    planet_index_offset = 1 + nuseptb; // 1 for central star

    read_scale();
    // process_planet_info();
    read_time_params(argc, argv);

    #ifdef OPENMP
    int total_thread = omp_get_num_procs();
    use_thread = total_thread / nhost;
    omp_set_num_threads(use_thread);
    printf("Using OpenMP with %d thread", use_thread);
    #endif

    open_host_ptb_file();
    start_simulation();
    close_host_ptb_file();

    // fclose(console_stdout);
    // fclose(console_stderr);

    return 0;
}