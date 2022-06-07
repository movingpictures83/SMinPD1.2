#ifndef SMINPD12PLUGIN_H
#define SMINPD12PLUGIN_H

#include  <string>
using std::string;

#define INITIAL_CODON_LIST_SIZE 20
#define CODON_LIST_SIZE_MULTIPLIER 2
#define DEFAULT_OFFSET 0
#define INITIAL_LINE_SIZE 500
#define LINE_SIZE_MULTIPLIER 2
#define MARKER_SEPARATORS ","
#define CODON_SIZE 3
/* Period parameters */  
#define TOLERANCE 1.0e-10
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)
#include "Plugin.h"
class SMinPD12Plugin : public Plugin {
   public:
	   void input(std::string);;
	   void run();;
	   void output(std::string);;

   private:
	           FILE    *fp1;
        char    *cc; /*BC: removed file name declarations from this line*/
        char    temp[LINELIMIT], temp2[MAXLENGHT], timechar[4]; /*BC: timechar[] should have length 4*/
        char* ch = NULL;
        char* ch2 = NULL;
        int             i, k, n,j,win_Count, fr_Dim, align, opt, res_no =0;
        double  **dist, ***dist_frags;
        double  GOP, GEP, match, mismatch,  mutrate;
        unsigned int    maxDim;
        Mintaxa *minresults[NUMSEQS];
        Fasta **temp_seqs;
        int             wSize=0, stSize,numWins;
        int old_n;
        short int bootreps=100;
        int closestRel=0;
        int output_file = 0;
        char* Cluster_Suffix = ".clust.txt"; /*BC clusters suffix*/
        char* Distances_Suffix = ".d.txt"; /*BC distance suffix*/
        char* Transitions_Suffix = ".transitions.txt"; /*BC distance suffix*/
        char ClusterFileName[ MAX_FILE_NAME_LENGTH + 1000 ]; /*BC: The file name for the cluster file */
        char DistanceFileName[ MAX_FILE_NAME_LENGTH + 1000 ]; /*BC: The file name for the distance file */
        char InputFileName[ MAX_FILE_NAME_LENGTH]; /*BC: The file name for the input file */
        char OutputFileName[ MAX_FILE_NAME_LENGTH + 1000 ]; /*BC: The file name for the transition file */

int counter = 0;

double PCCThreshold=0;
int model=TN93;
int seed=-3;
double alpha=0.5;
int fullBootscan=false;
int distPenalty=0;
int BootThreshold=96;
int BootTieBreaker=1;
int crossOpt=1;
int Bootstrap=0;
int recOn = 1;
int reportDistances=0;
double gapPenalty = 1;/* 1 is default: means gap columns are ignored are not counted*/
int printBoot = 0;
char codonFile[1000];
char output_directory[1000];
int clustBoot=1;
int *codonList;
double threshold = 0.001; /*BC*/

Fasta *seqs[NUMSEQS]; /* Declare an array of Fasta structures */ 

/*BC: begin variables to support amino acid clustering*/
int RNA = 0;
int ClusterByAminoAcid = 0;
int NumberOfTriplets = 72; /*64 normal + 8 with gaps*/
char* DNATriplets[72] = {
	"ATT", "ATC", "ATA",
	"CTT", "CTC", "CTA", "CTG", "TTA", "TTG", "CT-",
	"GTT", "GTC", "GTA", "GTG", "GT-",
	"TTT", "TTC",
	"ATG",
	"TGT", "TGC",
	"GCT", "GCC", "GCA", "GCG", "GC-",
	"GGT", "GGC", "GGA", "GGG", "GG-",
	"CCT", "CCC", "CCA", "CCG", "CC-",
	"ACT", "ACC", "ACA", "ACG", "AC-",
	"TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "TC-",
	"TAT", "TAC",
	"TGG",
	"CAA", "CAG",
	"AAT", "AAC",
	"CAT", "CAC",
	"GAA", "GAG",
	"GAT", "GAC",
	"AAA", "AAG",
	"CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "CG-",
	"TAA", "TAG", "TGA"
};

char DNAAminoAcids[73] = {
	'I','I','I',
	'L','L','L','L','L','L','L',
	'V','V','V','V','V',
	'F','F',
	'M',
	'C','C',
	'A','A','A','A','A',
	'G','G','G','G','G',
	'P','P','P','P','P',
	'T','T','T','T','T',
	'S','S','S','S','S','S','S',
	'Y','Y',
	'W','W',
	'Q','Q',
	'N','N',
	'H','H',
	'E','E',
	'D','D',
	'K','K',
	'R','R','R','R','R','R','R',
	'Z','Z','Z'
};

FILE* ClusterFile = NULL;
/*BC: end variables to support amino acid clustering*/

/*BC: begin variables for preclustering*/
int Cluster = 0;
int PreCluster = 0;
int old_opt;
int old_recOn;

int z_rndu=137;
unsigned w_rndu=13757;

unsigned long mt[624]; /* the array for the state vector  */
int mti=625; /* mti==N+1 means mt[N] is not initialized */

int is_break( char ch );
void freeBKPNode(BKPNode *bkpLL);
int BitCriteria(int Num, int Bit);/* Is bit set in Num */
int TwoCrossovers(int Num, int MaxBit, int *Begin, int *End);/* Is bit set in Num */
double SumofSquaredValues(double ***values, int **values2, int points, int data, int s,double *m_x);
char ** MatrixInit(int s, int dim1, int dim2);
double JC69distance(int gaps, int matches, int al_len, double alpha);
double K2P80distance(int p1, int p2, int gaps, int matches, int al_len, double alpha);
double TN93distance(int pT, int pC, int pG, int pA, int p1, int p2, int gaps,int matches, int al_len, double alpha, int s, int t);
double PairwiseDistance(char **seqsChars,int t, int s, int windowSize, int model);
int DistOnly(Fasta **seqs, double **dist, double ***dist_frags, int maxDim, int n,int win_Count, int model, int windowSize, int step, int numWins);
int SaveMinResults(Mintaxa **minresults, int *resno, int a, int c, double dist, double div);
void OptionalRecPrintout(FILE *fp4,  RecRes *recRes,Fasta **seqs,int *Frag_Seq,FDisNode (*recSolutions)[3][2], int i, int l, int r, int p, int q, int f, int index, int last, int BKP, int win_Count,int align,int fr_size);
double FillArrayInOrder(int *arrayInOrder, double ***dist_frags, double **bootVals, int  *Frag_Seq, int s, int k, int noSeqs, int bootOpt);
BKPNode * ArrayToBKPNode( int  *F, int win_Count, int fr_size, int align, int bootOpt, int step);
BKPNode * GetRecSolutionsII( double ***dist_frags, double **bootVals, Fasta **seqs, int  *Frag_Seq, int noSeqs, int s, int win_Count, int min_seq, int fr_size, double min_seq_dist, int align, int bootOpt, int step, int *avgBS);
double FindBestinStretch( double ***dist_frags, double **bootVals, int  *Frag_Seq,  int bootOpt,int noSeqs, int s, int win_Count, int candidate,int lb,int ub,double largest);
BKPNode * GetRecSolutionsI(double ***dist_frags, double **bootVals, Fasta **seqs, int  *Frag_Seq, int noSeqs, int s, int win_Count, int min_seq, int fr_size, double min_seq_dist, int align, int bootOpt, int step,int *avgBS);
int PickSeqofLargerAvgDis(double ***dist_frags, int k_idx, int f_idx, int min_idx,  int s, int win_Count, int k, int f,double min_d);
int PickSeqMinDist(double ***dist_frags, int  *Frag_Seq, int min_idx,int noSeqs, int k_idx, int f_idx, int s, int win_Count, int k, int f,double min_d);
void PrintFragments(FILE *fp3, int  *Frag_Seq, Fasta **seqs, double **dist, double ***dist_frags, int win_Count, int i, int s, int align );
int PickCandidates( Fasta **seqs, int  *Frag_Seq, double ***dist_frags, int win_Count, int min_idx, int s, int startAS, double min_d);
BKPNode * Check4Recombination(int  *Frag_Seq,  Fasta **seqs, double **dist, double ***dist_frags, int startAS, int s, int win_Count,  int min_idx,  int align, int maxDim, int opt, int step, int windowSize , short int bootReps, double min_d);
double GetBootstrapValues(Fasta **seqs, double ***bootVals,  int s, int a, int bootreps);
double GetClusteredBootstrapValues(Mintaxa *mintaxa, int max_times, Fasta **seqs, double ***bootVals,  int s, int min_idx, int bootreps);
char* RNACodonToDNACodon(char* RNACodon);
char CodonToAminoAcid(char* Codon);
int CheckClustering(Mintaxa *mintaxa,int max_times,Fasta **seqs,int s,int min_idx);
void StoreBootstrapValues(int *boots, Mintaxa **minresults, int resno, int start);
int OutputResults12(int win_Count, int s, int startAS, double iBase, double min_d, int min_idx, int *times_max, double *div, BKPNode *bkpLL, double ***bootVals, Mintaxa *mintaxa,Mintaxa **minresults, Fasta **seqs, double **dist, int maxDim, int n, FILE *fp1, double mutrate, int fr_size, int *resno, int align, int avgBS,int bootreps);
int GetMinDist12(Mintaxa **minresults, Fasta **seqs, double **dist, double ***dist_frags, int maxDim, int n, char *File1, char *File2, double mutrate, int win_Count, int *resno, int align, int opt, int step, int winSize);
void BuildPartialNJTrees(Fasta **seqs, Mintaxa **minresults, int resno, double **dist, char *File1, int n);
int PrepareWeights(int *wmod,int bootreps, int windowSize);
int BootRip(short int bootreps,  int model, double alpha, int noSeqs, int windowSize, char *seqs,  double *dist,  int *wmod);
int GetBootNJdistances(int bootreps,int noSeqs, double *dstMat,Fasta **seqs, int distPenalty);
int DoBootscan(double ***dist_frags,int **plotVal,int numWins, Fasta **seqsBS,  int maxDim, int noSeqs,int stepSize, int windowSize,int closestRel, int opt, short int bootreps, int fullBootscan);
 int BootscanShrinkPool(Fasta **seqs, int  *Frag_Seq, int step, int windowSize , short int bootReps, int numSeqs, int s, int numWins, int maxDim, int bootThreshold, int spikeLen);
BKPNode * BootscanCandidates(Fasta **seqs, int  *Frag_Seq, int step, int windowSize , short int bootReps,int numSeqs, int s, int numWins, int min_seq, double min_seq_dist, int align, int maxDim, int opt, int *avgBS);
int CopyToSeqsBootscan(Fasta **seqsBS, Fasta **seqs,int numSeqs, int seqNo, int *maxBSDim);
int GetParentsFromAnswer(Fasta **seqs,char *temp,int *parents, int n);
int Parents_Missed(int *Frag_Seq,int numSeqs,int *parents,int p);
int GetAvgBootstrapDist(Fasta **seqs,double **dist, double ***bootVals, int n, int bootreps);
int Bootscan(char *File1, char *File2, double ***dist_frags, int numWins, Mintaxa **minresults,Fasta **seqs, double **dist, int maxDim, int n,int step, int windowSize, int closestRel, int *resno, short int bootreps, int align);
void ReadUntil(FILE *fv, char stopChar, char *what);
void PrintUsage();
void ReadParams12( char *inputFile, char *outputFile, int *opt, int *align,  int *wSize, int *stSize, char* ifile);
void CheckTimes(Fasta **seqs, int n);
void heapify_mintaxa ( Mintaxa **list , int newnode );
void heapsort_mintaxa ( Mintaxa **list, int last );

void heapify_Fasta ( Fasta **list , int newnode );
void heapsort_Fasta ( Fasta **list, int last );
void PrintTree(FILE *fv, TNode *node, Fasta **seqs,int *boots, int printboot );
void AddTipOrNode(int *w,int *p,int *max_w, int chidx, double brlen,TNode **NodeStorage, int Join);
void NJTree(TNode **NodeStorage, double **DistMatrix, int UBound, int OutgroupIn0, int *seqIds, int *max_w, int *w);
int MatricizeTopology(TNode *node, int *topoMatrix, int n, int **dist, double avgBrLen,int *max);
double GetAvgBrLen(TNode *node, int *nodes);


TNode* Reroot( TNode* node, int id, int Branch_Zero, Fasta** seqs );
int remove_old_root( TNode* root_node, TNode* old_root );
void reverse( TNode* old_root, TNode* node, TNode* old_child, double new_length );
TNode* find_parent( TNode* root_node, TNode* child_node );
TNode* find_node( TNode* node, int id );
void print_nodes( TNode* node, int level, Fasta** seqs );
int parse_markers_file( char* file_name, int** codon_list, int* number_of_codons );
void SetSeed (unsigned long seed);
double rndu (void);

double rndgamma (double s);
double zeroin(double ax, double bx, double (*f)(double x));
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);
};

#endif
