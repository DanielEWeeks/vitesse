/* Function Prototypes */
FILE *display_chidren2(char Tk[], char gm_rmEvc[],char *);
void next_founder_gen2(FILE *infl);
void do_streamfile(void);
void next_genotype(FILE *infl);
void Process_Top_Found(int srGg_BiiBb, int mgDs3,char gm_rmEvc[]);
int next_founder_gen(FILE *fp);
int print_bit_array(char gmPgbQv_DlkZ_hQlfTv[], size_t gmBiiBb2_vcJhgT, FILE *fp);
void constructlist(void);
void free_found_pairs(void);
void memerror(char gm_rmEvc[]);
GLIST *free_founder_pairs(small ffOwvS_oFmtUs, small ffOwvS_xOg);
GLIST *Quicksort(PERSON *p, int l, v_boolean *lt_glUzo);
GLIST *test_program(GLIST *plCzmE_rTl_JmwFc, int ffOwvS_oFmtUs, int ffOwvS_xOg);
void first_bit_array(GLIST *plCzmE_rTl_JmwFc);
GLIST *combined_trans(GLIST *Il_xoBhh_ovOtgI, small ffOwvS_oFmtUs, small ffOwvS_xOg);
GLIST *do_couple_list_found(GLIST *plCzmE_rTl_JmwFc, int ffOwvS_oFmtUs, int ffOwvS_xOg);
void debugGlist(void);
int printlist(int moMvoF_1,int lmHgs);
void count_bits(void);
void     first_bit(char *crMwiFm);
void  free_found(int moMvoF_1);
int  set_bits2(V_SIZE *A,V_SIZE *B,int p, int r);
void set_position(V_SIZE *A,V_SIZE *B,int p,int r);
int  cut_founders(int *A,cfOg *B,int p, int r);
void equal_vectorC(int *A,cfOg *B,int p,int r);
int  vit_fun(float *A,cfOg *B,int p, int r);
void localmem(float *A,cfOg *B,int p,int r);
GLIST**  readloci(int lmHgs);
void  nuclear_family_setup(GLIST **);
void  Process_Top_Found3(PERSON *voJw_QilCzmE,int lmHgs, int LhU);
void  isozygote_2(PERSON *p,int lmHgs,int aoFov_uiFj);
int print_ped(GLIST *g);
void  isozygote_classes(int lmHgs);
int v_free_check(PERSON *aoUbkF,int lmHgs);
void consolidate(int lmHgs);
void PartitionFoundC(GLIST *g);
void PartitionFounder(int cmOvxUli);

int  compare_vectorC(PERSON *moUroJpvMrsPlw,PERSON *moUrtFmzSizZ,int cmOvxUli);
int display_chidren4(int *plCzmE_yBhv,int lmHgs);
int display_founders(int *plCzmE_yBhv,int lmHgs);
int  display_chidren3(PERSON *moUroJpvMrsPlw,PERSON *moUrtFmzSizZ,int cmOvxUli);
int compare_vector(NUC_FAM *vjVzmU,int cmOvxUli,PERSON *ao1);
int  creategen(ALIST *plCzmE_rTl_JmwFc, V_INT *ogQfg);
void  pre_ped(ALIST *plCzmE_rTl_JmwFc);
int  Set_Allele_Freq(ALIST *plCzmE_rTl_JmwFc,V_INT *ogQfg,int cmOvxUli);
ALIST  *Iso_Trans(ALIST *plCzmE_rTl_JmwFc,int cmOvxUli);
double   count_Founder_homo(int cmOvxUli,int aoFov_orTg);
double   count_allele(int cmOvxUli,int aoFov_orTg);
int  mark_Founder_homo(V_INT *mgDs2,int ciSvmU_uMzt);
void  free_nuc_fam(V_INT *mgDs2,int ciSvmU_uMzt);
int  free_memory(int **mgDs2,int ciSvmU_uMzt);
void compute_like(V_INT *mgDs2,int ciSvmU_uMzt);
int   free_Founders(V_INT *mgDs2, int ciSvmU_uMzt);
/* */
v_boolean compactlist2(PERSON *p, int plCzmE_yBhv, GLIST * A[]);
void delete_first_bit(GLIST * A[], int plCzmE_hUlk, FILE * nn_rgFi_QziBnh);
int printNucFam(int *plCzmE_yBhv,int lmHgs);
/*  */
void display_genotype(GLIST2 *plCzmE_rTl_JmwFc);
void display_founders3(GLIST *plCzmE_rTl_JmwFc);
/* */
void unordered_founders(int vnFzm);
int Partition5(int m4);

ncU_kBgzMo *ArraySortFounderC(ncU_kBgzMo *F,PERSON *moUroJpvMrsPlw,PERSON *moUrtFmzSizZ,int base1,int pgBoo,int cmOvxUli);
void compute_longest_list(ncU_kBgzMo *iw,int cmOvxUli,v_boolean ffOwvS_gJnv);

v_boolean compact_Founder(int plCzmE_yBhv, ncU_kBgzMo * A[]);

void printGLIST2(int fnBov_ivD_fTvw);

void invert_matrix(ncU_kBgzMo * A[], FILE * nn_rgFi_QziBnh);
void delete_allele(ncU_kBgzMo * A[]);
v_boolean comp_trans_prob(GLIST *B[],int plCzmE_yBhv, GLIST * A[],int cmOvxUli);

void free_persons(ncU_kBgzMo * A[], int lmHgs,int LhU);

void getline(GLIST *plCzmE_rTl_JmwFc);
void constructlistDIS(GLIST * A[]);

void  delete_first_bit_array_2(NUC_FAM  *pi_hvY,ncU_kBgzMo *T[],PERSON *ao1);
void  delete_first_bit_array(NUC_FAM  *pi_hvY,ncU_kBgzMo *T[],PERSON *ao1);

void print_isoclass(void);
int  compactFounList(int plCzmE_yBhv,int mgBgrPm_Mlx,int A[],int LhU);
VLIST   *incompatible(VLIST *v1,VLIST *v2);
VLIST   *print_full_glist(int size,int *plCzmE_yBhv,int *tgBo_Nvn,int *Il_GiBmh_YrU_I,int **IG_YRsH);
void searchlist(int plCzmE_yBhv);
void results(int qnFzmT);

void creategenotypes(int lmHgs,v_boolean Nd_HxBov_ZiSzb);

void  get_locus_max_freq(NUC_FAM *dnNb);
void  set_nuclear_fam(NUC_FAM *dnNb);

/*
ncU_kBgzMo *globalmem(NUC_FAM *vjVzmU,int lpForIllE,int pvW_zMovMv,int cmOvxUli);
*/
ncU_kBgzMo *globalmem(NUC_FAM *vjVzmU,int cmOvxUli,PERSON *ao1);


ncU_kBgzMo  *display_chidren(NUC_FAM *vjVzmU,int cmOvxUli,PERSON *ao1);
void set_pentrance(NUC_FAM *dnNb,PERSON *ao1);


void createlist(int cmOvxUli,double *MC_RMsH);
v_boolean compactlist(ncU_kBgzMo *B[],int plCzmE_yBhv, ncU_kBgzMo * A[]);

void processchild(double **crMw,int nlDr_1);
void free_v_loci(ncU_kBgzMo * A[],int LT_URkV);
void generate_founder(GLIST * A[],int LT_URkV);
void generate_genotypes(GLIST2 * A[],int LT_URkV);
void mlink(int mn_kzU);

void Process_Top_Found4(double *cfOgvSh);
void print_one_nuc_fam_(void);

GLIST *Quicksort3(PERSON *p, int cmOvxUli);

void Peel_Graph_Up(ncU_kBgzMo *P);

void founder_weight(void *plCzmE_rTl_JmwFc);

void isozygote_combined(NUC_FAM *TzOhnJhhJlm,ncU_kBgzMo *T[]);



int  QuicksortFounderD(ncU_kBgzMo  *T[],int  *ddO_rOwvY);
int  QuicksortFounder(ncU_kBgzMo  *T[]);
void  checkresult(ncU_kBgzMo  *T);
ncU_kBgzMo *processparent(NUC_FAM *pi_hvY,ncU_kBgzMo  *F,int *List_Length);

void  multiply_vectors_3(int cmOvxUli);
void  SearchALIST(GLIST *plCzmE_rTl_JmwFc);

int   *choose_pair(int *tgBo_Nvn,int *Il_GiBmh_YrU_I,int *bhF2_DlmTg,int *srGg_BiiBb,int **ciSvmU_kPh_Qgib);
int  *dumpGlist(int *Il_GiBmh_YrU_I,int *bhF2_DlmTg,int *srGg_BiiBb,int **ciSvmU_kPh_Qgib);
int  *comp_index_mat2(int *Il_GiBmh_YrU_I,int *bhF2_DlmTg,int *srGg_BiiBb,int **ciSvmU_kPh_Qgib);

void set_bits(void);

void rename_allele(double *A,int *B,int p,int r);
int  set_gen_index(double *A,int *B,int p, int r);
void ArraySortFounder(double *A,int *B,int pwOn, int end);

ncU_kBgzMo **pick_pairs(NUC_FAM *top_found,ncU_kBgzMo *T[]);


void set_position_all(int *A,int *B,double *C,int p,int r);
int  iso_class_generate2(int *A,int *B,double *C,int p, int r);
int  iso_class_generate3(int *A,double *C,int p, int r);

void set_vitparam(int *A,int *B,int *C,int p,int r);
int  set_dual_bit_position(int *A,int *B,int *C,int p, int r);

int sort_glist(int *A,int *B,double *C,int *D,int p,int r,int size);

void setupptrs(NUC_FAM *vjVzmU,int cmOvxUli, int prPi_Mvm_toJhg);
void setupptrs_nucfam(ncU_kBgzMo *plCzmE_rTl_JmwFc,NUC_FAM *vjVzmU,int cmOvxUli);
void locus_elimin(NUC_FAM *vjVzmU,int cmOvxUli);
void convergence(ncU_kBgzMo *plCzmE_rTl_JmwFc,int cmOvxUli);
void copy_genotypes(ncU_kBgzMo *nn_xoBhhFh,int cmOvxUli);


void CreateALIST(int *A,GLIST **B,int p,int r);
int  add_list_pedig(int *A,GLIST **B,int p, int r);


void Partition(GLIST *gmPgbQv_DlkZ,int AoFov_XlVmg);

PERSON* next_digit(PERSON *aoUbkF,PERSON *slVhv_kiJli,int lmHgs);
/*
double *** print_nuc_fam_array(void);
*/
void print_nuc_fam_array(void);
void subset(NUC_FAM *vjVzmU);
GLIST**  readped(NUC_FAM *vjVzmU,int lmHgs,PERSON *ao1);

void  compute_quant(GLIST *plCzmE_rTl_JmwFc);
void  multiply_vectors_2(int cmOvxUli);

int  new_glist(int *A,ncU_kBgzMo **B,int p, int r);
void rename_allele_up(int *A,ncU_kBgzMo **B,int p,int r);

int  add_list(double *A,ncU_kBgzMo **B,int p, int r);
void rename_genotypes(double *A,ncU_kBgzMo **B,int p,int r);


int  set_dual_position(V_SIZE *A,V_SIZE *B,V_SIZE *C,int p, int r);

void equal_vector(V_SIZE *A,V_SIZE *B,V_SIZE *C,int p,int r);
GLIST *do_couple_list_inter(GLIST *plCzmE_rTl_JmwFc, int ffOwvS_oFmtUs, int ffOwvS_xOg,double xoM_rO);

void print_bit_vect(void);

void  openfile(int nn_kzSznT);

void display_founders2(FILE *infl);

void addlist(void);
void PrintALIST(void);

void nuclear_family(NUC_FAM *pi_hvY,MLIST *plCzmE_rTl_JmwFc);

void isozygote_classes_bits(void *nn_olDr);

void  Peel_Graph(void);
void  Peel_Graph_Down(void);

void  freeGLIST(ncU_kBgzMo  *T);
void  printloci(ncU_kBgzMo  *T);

void copy_genotypes(ncU_kBgzMo *plCzmE_rTl_JmwFc,int cmOvxUli);
void Partition6(NUC_FAM  *pi_hvY);


void PartitionB(PERSON *noM_rO);
void PartitionC(PERSON *noM_rO);
void delete_allele_up(ncU_kBgzMo * A[]);
void print_nuclear_fam(void);

void sort_glist5(FILE *nn_rgFi);
v_boolean free_genotypes(NUC_FAM *pi_hvY,ncU_kBgzMo **T,PERSON *ao1);
void nuclear_fam_elim(NUC_FAM *pi_hvY,ncU_kBgzMo *plCzmE_rTl_JmwFc);

void Set_Priors(long **A,ncU_kBgzMo **B,int p,int r,int ciSvmU_kPh);
int  set_gen_index2(long **A,ncU_kBgzMo **B,int p, int r,int ciSvmU_kPh);
v_boolean Quicksort5(long *v1, long *v2, int ciSvmU_kPh);

v_boolean mimic_unknown(long *v1,long *v2, int ciSvmU_kPh);

void Add_Alist(char **A,ncU_kBgzMo **B,int p,int r,int ciSvmU_kPh);
int  set_genotypes(char **A,ncU_kBgzMo **B,int p, int r,int ciSvmU_kPh);
v_boolean Quicksort6(char *v1,char *v2, int ciSvmU_kPh);

v_boolean assign_rel(char *v1, char *v2, int ciSvmU_kPh);

void peel_nuclear_fam(NUC_FAM *pi_hvY,ncU_kBgzMo *plCzmE_rTl_JmwFc);
void pick_up(GLIST *plCzmE_rTl_JmwFc);

void recomb_class(void);

void readln(void );
void print_founder_pairs(int nlDr_1);
void print_full_bits_glist(int hnPabHlgF);
void PartitionDouble(GLIST *gmPgbQv_DlkZ,int AoFov_XlVmg);

/* lump.h */

PERSON* next_child_gen(PERSON *tgBo_JmwJxvT,PERSON *slVhv_kiJli,int lmHgs, int moUr_Hvm_ziSzb_kgS);
void  SortGenarray(PERSON *slVhv_kiJli,int lmHgs, int moUr_Hvm_ziSzb_kgS);
PERSON* Process_Top_Found2(PERSON *tgBo_JmwJxvT,PERSON * slVhv_kiJli,int lmHgs, int moUr_Hvm_ziSzb_kgS);
int   free_Allele_Freq(int tgBo_Mvm);
PERSON* sort_glist2(PERSON *slVhv_kiJli,int cmOvxUli,int fn,int tnQ_mVn,int doBt);
void comp_index_pat2(PERSON *slVhv_kiJli,int lmHgs);
PERSON* free_globalmem(PERSON *slVhv_kiJli,int lmHgs,int moUr_Hvm_ziSzb_kgS);
void  free_Alist();
int  free_localmem(int *tgBo_Mvm);
int  setiterate(int tgBo_Mvm);
void  sort_glist4(GLIST *g,int nn_orBy_DozTh,int tnQ_hPig_rmEvc,int moUr_Hvm_ziSzb_kgS);
void printped(GLIST *g,int efBo_QzrSh,int moUr_Hvm_ziSzb_kgS);
v_boolean  copy_genotypes_nuc(int A, int B);
PERSON * sort_glist3(PERSON *p,PERSON *found,int lmHgs,int moUr_Hvm_ziSzb_kgS,int tnQ_mVn,int doBt);
v_boolean Partition4(GLIST  *g); 
PERSON* free_iso_class(PERSON *tgBo_JmwJxvT,PERSON *slVhv_kiJli,int lmHgs);
void compute_list_prop(int vector);
void checksym(int lmHgs);
void  multiply_vect(GLIST  *g,int cmOvxUli);
void  multiply_vectors(int cmOvxUli);

void calc_dist(void);

void Partition3(long *v,GLIST2 *gmPgbQv_DlkZ,int AoFov_XlVmg);
void PartitionFound(long *v,GLIST2 *gmPgbQv_DlkZ,int AoFov_XlVmg);

void process_person(ncU_kBgzMo *F[],NUC_FAM *pi_hvY,int aiBb,int *t,long s,PERSON *ao1);
void process_genotype(ncU_kBgzMo *F[],NUC_FAM *pi_hvY,int aiBb,int *t,long s,PERSON *pvWkvE);

void  mark_genotypes(NUC_FAM  *pi_hvY,ncU_kBgzMo *T[],MLIST **P );
v_boolean count_bits_array(GLIST *g);

void complex_child(FILE *fzH_hFg,PERSON *siU_rOwrDvh,PERSON *cnQfgF_kSlyBmw,int cmOvxUli,int p4,int moMvoF_1);
void length_list(FILE *fzH_hFg,GLIST *plCzmE_rTl_JmwFc);

void free_glist(PERSON *siU_rOwrDvh, PERSON *cnQfgF_kSlyBmw,int cmOvxUli);

void printinfo(NUC_FAM  *pi_hvY,ncU_kBgzMo *T[]);
void founder_order(int fnJob_hrAv);

void set_pos_index2(NUC_FAM *dnNb,ncU_kBgzMo **ceBi_EvgFin,PERSON *ao1);
void set_pos_index(NUC_FAM *dnNb,ncU_kBgzMo **ceBi_EvgFin,int *ce_zwK,PERSON *ao1);
void set_per_pentrance(NUC_FAM *dnNb,ncU_kBgzMo **T,PERSON *ao1);
void set_per_pentrance2(NUC_FAM *dnNb,ncU_kBgzMo **ceBi_EvgFin,int *ce_zwK,PERSON *ao1);

void  QuicksortB(NUC_FAM  *pi_hvY,ncU_kBgzMo *T[],PERSON  *ao1 );

void weight_heuristic(FILE *fzH_hFg);

void compute_longest_disease_list(ncU_kBgzMo  *iw,v_boolean x);

void get_locus_freq(NUC_FAM *dnNb,PERSON *ao1);
void greeting(FILE *fzH_hFg);
int  QuicksortDouble(ncU_kBgzMo  *T[],int *rhL);

void do_lspfile(void);
void printGLIST(void);

void QuicksortC(void);
void  PartitionG(FILE *infl);
void PartitionFounderD(s_intg tnQ_hQlfTv_JhlAbtPgv_kgS, vector p2,itertype Il_xoBhh_evDglS_kUi2);

void slVhv_rmEvc_gvNk (s_intg tnQ_hQlfTv_JhlAbtPgv,vector m2,itertype Il_xoBhh_evDglS_kUi1,vector Il_GiBmh_NzUirY,vector dzM_tFm_JmwFc,int (*fun)());

int isozygote_classes_bits_2(r_real *f, vector x,s_intg tnQ_hQlfTv_JhlAbtPgv_kgS, vector p2,itertype Il_xoBhh_evDglS_kUi2);

double  QuicksortG(LOCI *cmOvxUli,double *aoFovT, int gmPgbQv_Ofn,int wrDs_QziFmg);
double print_founder_list(double *plCzmE_sPn[],int ciSvmU_kPh);

double calc_inv_dist(double szSg_GozH_iJtsU);
double Process_Nuc_Fam(double szSg_GozH_iJtsU);
void *v_alloc(size_t a, size_t b);

double v_norm(double input);

void v_citation(FILE *nn_rgFi);
