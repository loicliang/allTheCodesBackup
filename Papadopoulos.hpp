
extern"C" {
   void F2C_Papadopoulos_POLYCRYSTAL(int ninc, int ngrains, int nsys, double sig[ngrains][ninc][6], double n[ngrains][nsys][3], double l[ngrains][nsys][3], double frac[ngrains], double alpha, double *crit, double gcrit[ngrains]);
   void F2C_Papadopoulos_Macro(int ninc, double sig[ninc][6], double alpha, double *crit);
}


