#include <stdio.h>

#define bool int
#define True 1
#define False 0
#define MAX_PATH_LENGTH 2048
#define N_PTB_MAX 100
#define N_PARTICLE_MAX 2000
#define MAX_LINE_LENGTH 1024
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

// Array length limit: 2e9 (range of int). So normally, it's fine to use int for array subscripts

const double Mdot_to_Me = 332949.9009605902, parsec_to_au = 206264.806, PI = 3.141592653589793, kmps_to_earth_system_speed = 0.03357294283791669, myr_to_earth_system_time = 6283185.307179586, yr_to_earth_system_time = 6.283185307179586;

// solar system planets
// const int nplanet_total = 8;
// const double mlist_me_all[] = {0.05527, 0.815, 1, 0.10745, 317.83, 95.159, 14.5, 17.204};
// const double alist_AU_all[] = {0.38709893, 0.72333199, 1.00000011, 1.52366231, 5.20336301, 9.53707032, 19.19126393, 30.06896348};
// const double elist_all[] = {0.20563069, 0.00677323, 0.01671022, 0.09341233, 0.04839266, 0.0541506, 0.04716771, 0.00858587};
// const double ilist_deg_all[] = {7.00487, 3.39471, 0.00005, 1.85061, 1.3053, 2.48446, 0.76986, 1.76917};
// const double node_deg_all[] = {48.33167, 76.68069, -11.26064, 49.57854, 100.55615, 113.71504, 74.22988, 131.72169};
// const double peri_deg_all[] = {77.45645, 131.53298, 102.94719, 336.04084, 14.75385, 92.43194, 170.96424, 44.97135};
// const double L_deg_all[] = {252.25084, 181.97973, 100.46435, 355.45332, 34.40438, 49.94432, 313.23218, 304.88003};
// const double period_in_earth_system[] = {0.03833122, 0.09790866, 0.15915494, 0.29921129, 1.88757763, 4.68870462, 13.37060677, 26.22714307}; // for dt calculation
// // apprixomate period (yr): 0.24, 0.62, 1, 1.9, 12, 29, 84, 164

// // 50-AU-Jupiter only
// const int nplanet_total = 1;
// const double mlist_me_all[] = {317.83};
// const double alist_AU_all[] = {50.0};
// const double elist_all[] = {0.04839266};
// // const double ilist_deg_all[] = {1.3053}; // 0.022689 rad
// const double ilist_deg_all[] = {0.0}; 
// const double node_deg_all[] = {100.55615};
// const double peri_deg_all[] = {14.75385};
// const double L_deg_all[] = {34.40438};
// const double period_in_earth_system[] = {56.2799};

// double mlist_me[8], alist_AU[8], elist[8], ilist_deg[8], node_deg[8], peri_deg[8], L_deg[8], mlist_mdot[8], ilist_rad[8], node_rad[8], peri_rad[8], L_rad[8];

// // kuiper (comet) data
// double m_kuiper = 0.0; // required by REBOUND test massless particle(MLP). REBOUND will NOT automatically set 0 mass to testparticles
// // double a_kuiper[] = {150., 250., 350., 450., 550., 650., 750., 850., 950., 1050};
// double a_kuiper[] = {20., 25., 30., 35., 40., 45., 55., 60., 65., 70.,
//                     75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125.,
//                     130., 135., 140., 145.,
//                     150., 175., 200., 225., 250., 275., 300., 325., 350.,
//                     375., 400., 425., 450., 475., 500., 525., 550., 575.,
//                     600., 625., 650., 675., 700., 725., 750., 775., 800.,
//                     825., 850., 875., 900., 925., 950., 975., 1000., 1025.,
//                     1050., 1075., 1100., 1125., 1150., 1175., 1200., 1225., 1250.,
//                     1275., 1300., 1325., 1350., 1375., 1400., 1425., 1450., 1475.,
//                     1500., 1525., 1550., 1575., 1600., 1625., 1650., 1675., 1700.,
//                     1725., 1750., 1775., 1800., 1825., 1850., 1875., 1900., 1925.,
//                     1950., 1975., 2000.}; // 100 positions
// double e_kuiper = 0.01;
// double i_kuiper = 0.01;
// // other three angles are random. See main program


#define NROW_LPS_HOST_MAX 10000000
// #define NROW_LPS_PTB_MAX 100000000
// struct LpsoutHostFile
// {
//     double time[NROW_LPS_HOST_MAX];
//     long int index[NROW_LPS_HOST_MAX], name[NROW_LPS_HOST_MAX];
//     double mass[NROW_LPS_HOST_MAX];
//     double x1[NROW_LPS_HOST_MAX], x2[NROW_LPS_HOST_MAX], x3[NROW_LPS_HOST_MAX];
//     double v1[NROW_LPS_HOST_MAX], v2[NROW_LPS_HOST_MAX], v3[NROW_LPS_HOST_MAX];
//     long int nline;
//     // header: "Time(Myr), index, name, mass, x1, x2, x3, v1, v2, v3"
// };

// struct LpsoutPtbFile
// {
//     double time[NROW_LPS_PTB_MAX];
//     long int index[NROW_LPS_PTB_MAX], name[NROW_LPS_PTB_MAX];
//     double mass[NROW_LPS_PTB_MAX];
//     double x1[NROW_LPS_PTB_MAX], x2[NROW_LPS_PTB_MAX], x3[NROW_LPS_PTB_MAX];
//     double v1[NROW_LPS_PTB_MAX], v2[NROW_LPS_PTB_MAX], v3[NROW_LPS_PTB_MAX];
//     long int nline;
// };

struct nbparticle
{
    double time;
    long int index, name, linen;
    double mass, x1, x2, x3, v1, v2, v3;
};

typedef struct {
    double m_me;
    double a_au;
    double e;
    double i_deg;
    double peri_deg;
    double node_deg;
    double M_deg;
} inputparticle;

// const unsigned int RANDSEEDS[] = {8297, 7474, 43265, 41182, 3158, 48012, 43358, 94073, 59305, 64405, 5536, 38289, 16674, 2837, 35065, 47955, 7825, 87285, 78460, 71837, 81540, 1182, 968, 68649, 32040, 722, 97682, 69611, 82476, 34132, 26585, 78236, 30724, 63995, 34291, 91129, 17812, 95835, 23229, 50174, 15536, 86055, 15843, 88928, 99311, 26593, 83398, 24798, 90170, 70439, 35943, 40161, 17307, 35677, 13507, 17107, 23837, 34557, 21561, 86651, 78627, 46586, 79168, 96405, 36049, 31311, 45184, 38279, 76998, 93591, 73946, 79386, 94490, 37888, 66282, 51560, 81067, 64213, 58892, 54696, 11579, 91743, 13339, 34824, 95143, 9343, 74677, 61480, 52244, 18038, 14714, 45930, 88521, 99597, 36175, 95466, 68340}; //97 numbers

double getRandom(const double min, const double max) {
    static unsigned int theseed = 1;
    srand(theseed++);
    if (theseed == UINT_MAX) theseed = 0;
    return rand() * ((max - min) / RAND_MAX) + min;
}

const double comet_density = 0.6; // g/cm3
const double jupyter_density = 0.6; // g/cm3

int small_to_big_double(const void *a, const void *b){
    return (*(double *)a - *(double *)b);
}

int get_output(char *cmd, char *result)
{
    char buf[1024];
    FILE *fp;

    if ((fp = popen(cmd, "r")) == NULL) {
        printf("Error opening pipe!\n");
        return -1;
    }

    while (fgets(buf, 1024, fp) != NULL)
        strcat(result, buf);

    if (pclose(fp)) {
        printf("Command not found or exited with error status\n");
        return -1;
    }

    return 0;
}

size_t trimwhitespace(char *out, size_t len, const char *str)
{
  if(len == 0)
    return 0;

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace((unsigned char)*str)) str++;

  if(*str == 0)  // All spaces?
  {
    *out = 0;
    return 1;
  }

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = (end - str) < len-1 ? (end - str) : len-1;

  // Copy trimmed string and add null terminator
  memcpy(out, str, out_size);
  out[out_size] = 0;

  return out_size;
}