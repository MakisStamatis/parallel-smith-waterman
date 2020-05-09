#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <ctype.h>
#include <cstring>
#include <stdarg.h>
#include <sys/time.h>
#include <omp.h>

#define THREADS 4
#define COUT 1

/* Size of Blocks */
//#define BLOCK_J 200
//#define BLOCK_I 200

int BLOCK_I ;
int BLOCK_J ;
FILE * f;

/* Variables that used for block status*/
#define VOID 0
#define CREATED 1
#define DONE 2

using namespace std;

/* We split in functions to view the precent of time of each funtions 
    all function are declared inline */
float  similarity_score(char a, char b);
float  find_array_max(float array[], int length);
void   checkfile(int open, char filename[]);
string read_sequence(ifstream& f);
void   my_H_malloc(void);
void   my_B_malloc(void);
void   my_Ii_malloc(void);
void   my_Ij_malloc(void);
void   my_H_findmax(void);
void   my_actual_alg(int i_block, int j_block, int max_i_in_block, int max_j_in_block);
void   my_backtracking(void);
void   my_H_print(void);
int    main_function(int argc, char **argv);

/* The number of blocks in Horizontal Axis and in Vertical Axis*/
int num_blocks_in_j=0;
int num_blocks_in_i=0;

/* This is a 2d table with the size of all blocks 
    and contains information for each block :
     VOID
     CREATED
     DONE
*/
short** block_status;

/* Maximun rows and collumns of each block */
int max_i_in_Block;
int max_j_in_Block;

string seq_a,seq_b;

/* Size of Vertical and Horizontal Axis of table H which contains aligment socres */
int N_i ;
int N_j ;


unsigned short **I_i, **I_j;
float mu, delta, **H;

float H_max = 0.;
int i_max=0,j_max=0;

int current_i , current_j ;
int next_i ;
int next_j ;
int tick = 0;

struct timeval	StartTime, EndTime;



int main(int argc, char** argv){
    /* log file */
    f = fopen("version3.2.log","a+");
    BLOCK_I = 200;
    BLOCK_J = 200;

    /* Multiple execution for testing and debuging */
//    for(int i=0;i<5;i++){
        main_function(argc,argv);
//    }


}


int main_function(int argc, char** argv) {

    #ifdef _OPENMP
    omp_set_num_threads((int)THREADS);
    cout<<"OPENMP defined OK"<<endl;
    #endif

    // read info from arguments
    if(argc!=6) {
        #if COUT
        cout<<"Give me the propen number of input arguments:"<<endl<<"1 : mu"<<endl;
        cout<<"2 : delta"<<endl<<"3 : filename sequence A"<<endl<<"4 : filename sequence B"<<endl;
        cout<<"5 : maximal length N of sequences"<<endl;
        #endif
        exit(1);
    }

    mu    = atof(argv[1]);
    delta = atof(argv[2]);

    char *nameof_seq_a = argv[3];
    char *nameof_seq_b = argv[4];

    // read the sequences into two vectors:
    ifstream stream_seq_a;
    stream_seq_a.open(nameof_seq_a);
    checkfile(! stream_seq_a,nameof_seq_a);
    #if COUT
    cout << "Reading file \"" << nameof_seq_a << "\"\n";
    #endif
    seq_a = read_sequence(stream_seq_a);
    #if COUT
    cout << "File \"" << nameof_seq_a << "\" read\n\n";
    #endif

    ifstream stream_seq_b;
    stream_seq_b.open(nameof_seq_b);
    checkfile(! stream_seq_b,nameof_seq_b);
    #if COUT
    cout << "Reading file \"" << nameof_seq_b << "\"\n";
    #endif
    seq_b = read_sequence(stream_seq_b);
    #if COUT
    cout << "File \"" << nameof_seq_b << "\" read\n\n";
    #endif
    N_i = seq_a.length();
    N_j = seq_b.length();
    #if COUT
    cout << "First sequence has length  : " << setw(6) << N_i <<endl;
    cout << "Second sequence has length : " << setw(6) << N_j << endl << endl;
    cout << "Allocating memory for matrix H\n";
    #endif

    /* Compute the number of blocks */
    num_blocks_in_j = (N_j - 1) / BLOCK_J;
    num_blocks_in_i = (N_i - 1) / BLOCK_I;

    /* Compute the Rows and Collumns in case of smaller files than BLOCK_J and BLOCK_I  */
    max_i_in_Block = num_blocks_in_i > 1 ? BLOCK_I : N_i - 1;
    max_j_in_Block = num_blocks_in_j > 1 ? BLOCK_J : N_j - 1;

    /* Function for allocation and initialization of block_status table */
    my_B_malloc();

    /* Function for allocation and initialization of H table */
    my_H_malloc();

    /* Function for allocation and initialization of Ii table */
    my_Ii_malloc();

    /* Function for allocation and initialization of Ij table */
    my_Ij_malloc();

    /* Start the timer */
    gettimeofday(&StartTime, NULL);

    /* Starting the parallel secion */
    #pragma omp parallel
    {
        /* The First call of actual algorithm going to be executed in the master thread */
        #pragma omp master
        my_actual_alg(0,0,max_i_in_Block,max_j_in_Block);
        #pragma omp taskwait
    }

    /* Function that finds the maximum score in table H */
    my_H_findmax();

    /* Function that fill the sequense with back tracking */
    my_backtracking();


    /* Conversion of time in secs */
    if (EndTime.tv_usec < StartTime.tv_usec) {
        int nsec = (StartTime.tv_usec - EndTime.tv_usec) / 1000000 + 1;
        StartTime.tv_usec -= 1000000 * nsec;
        StartTime.tv_sec += nsec;
    }
    if (EndTime.tv_usec - StartTime.tv_usec > 1000000) {
        int nsec = (EndTime.tv_usec - StartTime.tv_usec) / 1000000;
        StartTime.tv_usec += 1000000 * nsec;
        StartTime.tv_sec -= nsec;
    }

    printf("\n\nParallel calculation time: %ld.%.6ld seconds\n", EndTime.tv_sec  - StartTime.tv_sec, EndTime.tv_usec - StartTime.tv_usec);
    
    /* Prints for Debuging */
    printf("\nThe H[%d][%d] = %f\n",i_max,j_max,H_max);
    printf("\nThe N _a= %d\t N_b=%d\n",N_i,N_j);
    printf("\nThe B_i= %d\t B_j=%d\n",BLOCK_I,BLOCK_J);
    fprintf(f,"%ld.%.6ld\t%d\t%d\n",(EndTime.tv_sec  - StartTime.tv_sec), (EndTime.tv_usec - StartTime.tv_usec), BLOCK_I, BLOCK_J);

    return 0;
}

/******************************************************************************/
/* auxiliary functions used by main                                          */
/******************************************************************************/

void checkfile(int open, char filename[]) {

    if (open) {
        cout << "Error: Can't open the file "<<filename<<endl;
        exit(1);
    }
    else cout<<"Opened file \"" << filename << "\"\n";
}

/******************************************************************************/

inline float similarity_score(char a,char b) {

    float result;
    if(a==b) {
        result=1.;
    }
    else {
        result=-mu;
    }
    return result;
}

/******************************************************************************/

string read_sequence(ifstream& f){
    string seq;
    char line[20000];
    while( f.good() )
    {
        f.getline(line,20000);
        if( line[0] == 0 || line[0]=='#' )
            continue;
        for(int i = 0; line[i] != 0; ++i)
        {
            int c = toupper(line[i]);
            if( c != 'A' && c != 'G' && c != 'C' && c != 'T' )
                continue;

            seq.push_back(char(c));
        }
    }
    return seq;
}

/******************************************************************************/
inline void my_H_malloc(){
    H = (float **)malloc((N_i + 1) * sizeof(float *));
    if (H == NULL) {
        #if COUT
        cout << "Could not allocate memory for matrix H\n";
        #endif
        exit(1);
    }
    #pragma omp parallel for
    for (int i = 0; i < (N_i + 1); i++) {
        H[i] = (float *)malloc((N_j + 1) * sizeof(float));
        if (H[i] == NULL) {
            #if COUT
            cout << "Could not allocate memory for matrix H[" << setw(6) << i << "]\n";
            #endif
            exit(1);
        }
    }
    #if COUT
    cout << "Memory for matrix H allocated\n\n";

    cout << "Initializing matrix H\n";
    #endif

    for(int i=0; i<=N_i; i++) {
        #pragma omp parallel for
        for(int j=0; j<=N_j; j++) {
            H[i][j]=0.;
        }
    }
    #if COUT
    cout << "Matrix H initialized\n\n";
    #endif
    return;
}

/******************************************************************************/
inline void my_Ii_malloc(){
    I_i = (unsigned short **)malloc((N_i + 1) * sizeof(unsigned short *));
    if (I_i == NULL) {
    #if COUT
        cout << "Could not allocate memory for matrix I_i\n";
    #endif
        exit(1);
    }
    #pragma omp parallel for
    for (int i = 0; i < (N_i + 1); i++) {
        I_i[i] = (unsigned short *)malloc((N_j + 1) * sizeof(unsigned short));
        if (I_i[i] == NULL) {
        #if COUT
            cout << "Could not allocate memory for matrix I_i[" << setw(6) << i << "]\n";
        #endif
            exit(1);
        }
    }
    #if COUT
    cout << "Memory for matrix I_i allocated\n\n";
    #endif
    
    return;
}

inline void my_B_malloc(){
    block_status = (short**)malloc(num_blocks_in_i * sizeof(short*));
    if (block_status == NULL) {
        cout << "Could not allocate memory for matrix I_i\n";
        exit(1);
    }
    #pragma omp parallel for
    for(int i=0 ; i< (num_blocks_in_i) ; i++) {
        block_status[i] = (short*) malloc (num_blocks_in_j * sizeof(short));
        if (block_status[i] == NULL) {
            cout << "Could not allocate memory for matrix I_i[" << setw(6) << i << "]\n";
            exit(1);
        }
    }
    for(int i=0; i<(num_blocks_in_i); i++)
    {
     #pragma omp parallel for
        for(int j=0;j<(num_blocks_in_j);j++){
         block_status[i][j]=VOID;
        }
     }
}

/******************************************************************************/
inline void my_Ij_malloc(){
    I_j = (unsigned short **)malloc((N_i + 1) * sizeof(unsigned short *));
    #if COUT
    cout << "Allocating memory for matrix I_j\n";
    #endif
    if (I_j == NULL) {
        #if COUT
        cout << "Could not allocate memory for matrix I_j\n";
        #endif
        exit(1);
    }
    #pragma omp parallel for
    for (int i = 0; i < (N_i + 1); i++) {
        I_j[i] = (unsigned short *)malloc((N_j + 1) * sizeof(unsigned short));
        if (I_j[i] == NULL) {
            #if COUT
            cout << "Could not allocate memory for matrix I_j[" << setw(6) << i << "]\n";
            #endif
            exit(1);
        }
    }

    #if COUT
    cout << "Memory for matrix I_j allocated\n\n";
    #endif
    
    return;
}

/******************************************************************************/
inline void my_H_findmax(){
    for(int i=1; i<=N_i; i++) {
        #pragma omp parallel for
        for(int j=1; j<=N_j; j++) {
            if(H[i][j]>H_max) {
                H_max = H[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }
    return;
}

/******************************************************************************/
inline void my_H_print(){
    #if COUT
    // Print the matrix H to the console
    cout<<"**********************************************"<<endl;
    cout<<"The scoring matrix is given by  "<<endl<<endl;
    for(int i=1; i<=N_i; i++) {
        for(int j=1; j<=N_j; j++) {
            cout<<H[i][j]<<" ";
        }
        cout<<endl;
    }
    #endif
}

/******************************************************************************/
inline void my_actual_alg(int i_block, int j_block, int max_i_in_block, int max_j_in_block){
    
    /* This are variables for finding the max value between 4 values in temp */
    float temp[4];
    float max;            
    int ind;

    /*Compute of starting and ending Indexes for Rows in block*/
    int BlockStartI=i_block * BLOCK_I + 1;
    int BlockEndI=BlockStartI + max_i_in_block - 1;

    /*Compute of starting and ending Indexes for Columns in block*/
    int BlockStartJ=j_block * BLOCK_J + 1;
    int BlockEndJ=BlockStartJ + max_j_in_block - 1;;

    /* Access into each Block */
    for(int i=BlockStartI; i<=BlockEndI; i++) {
        for(int j=BlockStartJ; j<=BlockEndJ; j++) {

            temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]);
            temp[1] = H[i-1][j]-delta;
            temp[2] = H[i][j-1]-delta;
            temp[3] = 0.;
            max = temp[0];
            ind = 0;
            /* Finding maximum value from temp vector */

            for(int x = 1; x<4; x++) {
                if(temp[x] > max) {
                    max = temp[x];
                    ind = x;
                }
            }
            H[i][j] = max;

            switch(ind) {
            case 0:                                  // score in (i,j) stems from a match/mismatch
                I_i[i][j] = i-1;
                I_j[i][j] = j-1;
                break;
            case 1:                                  // score in (i,j) stems from a deletion in sequence A
                I_i[i][j] = i-1;
                I_j[i][j] = j;
                break;
            case 2:                                  // score in (i,j) stems from a deletion in sequence B
                I_i[i][j] = i;
                I_j[i][j] = j-1;
                break;
            case 3:                                  // (i,j) is the beginning of a subsequence
                I_i[i][j] = i;
                I_j[i][j] = j;
                break;
            }
        }
    }
    
    /* When the access in block ends only one thread
       has the right to set it DONE in block status table */
    #pragma omp critical
    {
        block_status[i_block][j_block] = DONE;
    }

    /* Check if we can start right block */

    if (j_block + 1 < num_blocks_in_j) {
        short prev_i_BlockStatus = DONE;
        short prev_j_BlockStatus = DONE;
        short prev_diag_BlockStatus = DONE;

        if((i_block - 1 >= 0) && (j_block - 1 >= 0)){
            prev_i_BlockStatus = block_status[i_block-1][j_block+1];
            prev_j_BlockStatus = block_status[i_block][j_block];
            prev_diag_BlockStatus = block_status[i_block-1][j_block];
        }

        /* Checking if the 3 blocks that are necessary for the block is DONE */
        if((prev_i_BlockStatus == DONE) && (prev_j_BlockStatus == DONE) && (prev_diag_BlockStatus == DONE)){
            short rightBlockStatus;

            /* When a block is CREATED only one thread must
                have access to it */
            #pragma omp critical
            {
                rightBlockStatus = block_status[i_block][j_block+1];
                if(rightBlockStatus == VOID) {
                    block_status[i_block][j_block+1] = CREATED;
                }
            }

            if(rightBlockStatus == VOID) {
                int new_max_j_in_block = max_j_in_block;

                if(j_block + 1 == num_blocks_in_j -1) {
                    new_max_j_in_block += (N_j) % BLOCK_J;
                }

                /* The nested call for the next block computation is a new task 
                    when a thread is avoilabe take this task */
                #pragma omp task private(ind,max,temp)
                my_actual_alg(i_block,j_block+1,max_i_in_block,new_max_j_in_block);
            }
        }
    }

    //Check if we need to start block below
    if(i_block + 1 < num_blocks_in_i) {
        short prev_j_BlockStatus = DONE;
        short prev_i_BlockStatus = DONE;
        short prev_diag_BlockStatus = DONE;

        if((i_block - 1 >= 0) && (j_block - 1 >= 0)){
            prev_j_BlockStatus = block_status[i_block+1][j_block-1];
            prev_i_BlockStatus = block_status[i_block][j_block];
            prev_diag_BlockStatus = block_status[i_block][j_block-1];
        }

        /* Checking if the 3 blocks that are necessary for the block is DONE */
        if((prev_i_BlockStatus == DONE) && (prev_j_BlockStatus == DONE) && (prev_diag_BlockStatus == DONE)){
            short belowBlockStatus;

            /* When a block is CREATED only one thread must
                have access to it */
            #pragma omp critical
            {
                belowBlockStatus = block_status[i_block+1][j_block];
                if(belowBlockStatus == VOID) {
                    block_status[i_block+1][j_block] = CREATED;
                }
            }

            if(belowBlockStatus == VOID) {
                int new_max_i_in_block = max_i_in_block;

                if(i_block + 1 == num_blocks_in_i - 1) {
                    new_max_i_in_block += (N_i ) % BLOCK_I;
                }

                /* The nested call for the next block computation is a new task 
                    when a thread is avoilabe take this task */                
                #pragma omp task private(ind,max,temp)
                my_actual_alg(i_block+1,j_block,new_max_i_in_block,max_j_in_block);
            }
        }
    }
}

/******************************************************************************/
inline void my_backtracking(){
    // Backtracking from H_max
    current_i = i_max, current_j = j_max;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    tick = 0;

    char consensus_a[N_i+N_j+2], consensus_b[N_i+N_j+2];

    while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)) {

        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
        else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
        else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

        current_i = next_i;
        current_j = next_j;
        next_i = I_i[current_i][current_j];
        next_j = I_j[current_i][current_j];
        tick++;
    }

    /* End of Timer */
    gettimeofday(&EndTime, NULL);
    #if COUT
    cout<<endl<<"***********************************************"<<endl;
    cout<<"The alignment of the sequences"<<endl<<endl;
    for(int i=0; i<N_i; i++) {
        cout<<seq_a[i];
    };
    cout<<"  and"<<endl;
    for(int i=0; i<N_j; i++) {
        cout<<seq_b[i];
    };
    cout<<endl<<endl;
    cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<endl<<endl;
    for(int i=tick-1; i>=0; i--) cout<<consensus_a[i];
    cout<<endl;
    for(int j=tick-1; j>=0; j--) cout<<consensus_b[j];
    cout<<endl;
    #endif
}

/******************************************************************************/
