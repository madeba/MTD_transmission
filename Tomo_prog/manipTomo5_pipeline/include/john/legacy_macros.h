
/*******************************************************************************
 Memory-allocation related: MATRIX
*******************************************************************************/


// section deprecated, kept for legacy AIR_Volume
// WAY easier to use 1D arrays


/* allocates to given MAT pointer enough memory for a matrix of 
   X * Y elements of type DATATYPE. Initialization to 0.
   - Requires <assert.h> */
#define MATRIX_ALLOC(MAT, X, Y, DATATYPE) \
assert(MAT = (DATATYPE **) malloc(X * sizeof(DATATYPE*))); \
{ \
  int CTR; \
  for (CTR = 0; CTR < X; CTR++) \
    {assert(MAT[CTR] = (DATATYPE *) calloc(Y, sizeof(DATATYPE)));} \
} 
/* modified from img.h: 
   calloc used: elements initialized to 0
   CTR declared on heap and not inherited 
   X and Y reverted: now mat[X-1][Y-1] exists */


/* frees memory allocated by MATRIX_ALLOC */
#define MATRIX_FREE(MAT, X) \
{ \
  int I; \
  for (I = 0; I < X; I++) \
    free(MAT[I]); \
  free(MAT); \
} 
/* this could be a (inline?) function since type-independent*/


/* allocates to given VOL pointer enough memory for a volume of 
   X * Y * Z elements of type DATATYPE. Initialization to 0.
   - Requires <assert.h> */
#define VOLUME_ALLOC(VOL, X, Y, Z, DATATYPE) \
assert(VOL = (DATATYPE ***) malloc(X * sizeof(DATATYPE**))); \
{ \
  int I, J; \
  for (I = 0; I < X; I++) \
    { \
      assert(VOL[I] = (DATATYPE **) calloc(Y, sizeof(DATATYPE*))); \
      for (J = 0; J < Y; J++) \
        assert(VOL[I][J] = (DATATYPE *) calloc(Z, sizeof(DATATYPE))); \
    } \
} 
/* to declare a volume in SEG convention:
   VOLUME_ALLOC(vol, dimz, dimx, dimy);
   vol[z][x][y]; */


/** a verif sur poly michel, je crois que type facultatif */
#define MATRIX_ALLOC_2(MYMAT, X, Y, DATATYPE) \
{ \
  int I; \
  DATATYPE **MAT; \
  \
  assert(MAT = (int **) calloc(X, sizeof(*MAT))); \
  assert(*MAT=(int *) calloc (X * Y, sizeof(**MAT))); \
  for (I = 1; I < X; I++) \
    MAT[I] = MAT[I-1] + Y; \
  \
  MYMAT = **MAT; \
}


/* frees memory allocated by VOLUME_ALLOC */
#define VOLUME_FREE(VOL, X, Y) \
{ \
  int I, J; \
  for (I = 0; I < X; I++) \
    { \
      for (J = 0; J < Y; J++) \
        free(VOL[I][J]); \
    } \
} 
/* this could be a (inline?) function since type-independent*/


/* allocates to given MAT pointer enough memory for a matrix of 
   y * x elements of type DATATYPE.
   - Requires extra integer variable INT_CTR.
   - Requires <assert.h> */
#define REV_MATRIX_ALLOC(MAT, X, Y, DATATYPE, INT_CTR) \
assert(MAT = (DATATYPE **) malloc(Y * sizeof(DATATYPE*))); \
for (INT_CTR = 0; INT_CTR < Y; INT_CTR++) \
    {assert(MAT[INT_CTR] = (DATATYPE *) malloc(X * sizeof(DATATYPE)));} 
/* Attention: colle probablement pas avec SEG */



/*******************************************************************************
 VOLUME
*******************************************************************************/



/* initializes volume to given value */
#define VOLUME_INIT(MAT, X, Y, Z, VAL) \
{ \
  int T_I, T_J, T_K; \
  for (T_I = 0; T_I < X; T_I++) \
    for (T_J = 0; T_J < Y; T_J++) \
      for (T_K = 0; T_K < Z; T_K++) \
        MAT[T_I][T_J][T_K] = VAL; \
}


/* transfer a whole slice (possibly of a volume, ie calling vol[k]) */
#define TRANSFER_SLICE(slice_in, slice_out, x, y) \
{ \
  int a, b; \
  for (a = 0; a < x; a++) \
    for (b = 0; b < y; b++) \
      slice_out[a][b] = slice_in[a][b]; \
}



/*******************************************************************************
 Iterators
*******************************************************************************/


/* sais plus a quoi, mais ca sert */
#define WITH_BREAK(CODE_EXEC)			\
{ \
  int I; \
  for (I = 0; I < 1; I++) \
    { \
	CODE_EXEC; \
    } \
}



/*******************************************************************************
 Time Measurement - C - stand alone
*******************************************************************************/


/** times (in seconds) the execution of given code. Displays given message 
    at start and end of execution.
    Time measured is actual time from system clock: it accounts for other 
    running processes (time()). See clock() for elapsed cycles 
    WARNING: do not nest calls. 
    WARNING: variables declared within CODE_EXEC will 
*/
#define TIME(MSG, CODE_EXEC) \
{ \
  time_t t_start; \
  time_t t_end; \
  double seconds_elapsed; \
  \
  t_start = time(NULL); \
  printf("\n" MSG ": START "); \
  \
  CODE_EXEC; \
  \
  t_end = time(NULL); \
  seconds_elapsed = difftime(t_end, t_start); \
  printf("\n" MSG ": COMPLETE in %5.1f seconds \n", seconds_elapsed); \
}
/* usage:
TIME("my message",
{
  code to time;
});
*/ 

/* same as before, but with two separate calls and without context */
#define TIME_START(NAME, MSG)	                                    \
  				                                    \
  time_t t_start_s_##NAME, t_end_s_##NAME;	                    \
  clock_t t_start_c_##NAME, t_end_c_##NAME;	                    \
  double seconds_elapsed_s_##NAME, seconds_elapsed_c_##NAME;        \
  				                                    \
  t_start_s_##NAME = time(NULL);	                            \
  t_start_c_##NAME = clock();	                                    \
  printf("\n" MSG ": START \n");	                            \


#define TIME_END(NAME, MSG)                                                  \
  t_end_s_##NAME = time(NULL);						     \
  t_end_c_##NAME = clock();						     \
  seconds_elapsed_s_##NAME = difftime(t_end_s_##NAME, t_start_s_##NAME);     \
  seconds_elapsed_c_##NAME = (double)(t_end_c_##NAME - t_start_c_##NAME) / CLOCKS_PER_SEC; \
  printf("\n" MSG ": COMPLETE in %f secs sys. time and %f secs cpu cycles\n", seconds_elapsed_s_##NAME, seconds_elapsed_c_##NAME); \


/* usage:
TIME_START(CHRONO1, "message");
// random code
TIME_END(CHRONO1, "message");
*/

#define CHRONO_SET(NAME)	                                    \
  time_t t_start_s_##NAME, t_end_s_##NAME;                          \
  clock_t t_start_c_##NAME, t_end_c_##NAME;			    \
  double t_acc_s_##NAME = 0, t_acc_c_##NAME = 0;                    \

#define CHRONO_START(NAME)					    \
  t_start_s_##NAME = time(NULL);	                            \
  t_start_c_##NAME = clock();	                                    \

#define CHRONO_STOP(NAME)					                           \
  t_end_s_##NAME = time(NULL);						                   \
  t_end_c_##NAME = clock();						                   \
  t_acc_s_##NAME += difftime(t_end_s_##NAME, t_start_s_##NAME);                            \
  t_acc_c_##NAME += t_end_c_##NAME - t_start_c_##NAME                                      \

#define CHRONO_END(NAME, MSG)                                                              \
  printf("\n" MSG ": COMPLETE in %f secs sys. time and %f secs cpu cycles\n", t_acc_s_##NAME, (double)(t_acc_c_##NAME) / CLOCKS_PER_SEC); \


/* variant with generic output */
#define CHRONO_END_OUT(NAME, MSG, STREAM)					\
  fprintf(STREAM, "\n" MSG ": COMPLETE in %f secs sys. time and %f secs cpu cycles\n", t_acc_s_##NAME, (double)(t_acc_c_##NAME) / CLOCKS_PER_SEC); \



/* same as previous but shows both system and cycle times*/
#define TIMES(MSG, CODE_EXEC)  \
{                              \
  time_t t_start_s, t_end_s;   \
  clock_t t_start_c, t_end_c;  \
  double seconds_elapsed_s, seconds_elapsed_c;	\
  \
  t_start_s = time(NULL);      \
  t_start_c = clock();         \
  printf("\n" MSG ": START "); \
  \
  CODE_EXEC; \
  \
  t_end_s = time(NULL); \
  t_end_c = clock();    \
  seconds_elapsed_s = difftime(t_end_s, t_start_s);           \
  seconds_elapsed_c = (double)(t_end_c - t_start_c) / CLOCKS_PER_SEC ; \
  printf("\n" MSG ": COMPLETE in %f secs sys. time and %f secs cpu cycles\n", seconds_elapsed_s, seconds_elapsed_c); \
}


/* variant of TIME_END giving results in deci-ms (10-4 secs) for accurate measurements */
#define TIME_DMS_END(NAME, MSG)                                                  \
  t_end_s_##NAME = time(NULL);						     \
  t_end_c_##NAME = clock();						     \
  seconds_elapsed_s_##NAME = difftime(t_end_s_##NAME, t_start_s_##NAME);     \
  seconds_elapsed_c_##NAME = (double)(t_end_c_##NAME - t_start_c_##NAME) / ( CLOCKS_PER_SEC / 10000 ); \
  printf("\n" MSG ": COMPLETE in %f secs sys. time and %f 10-4 secs cpu cycles\n", seconds_elapsed_s_##NAME, seconds_elapsed_c_##NAME); \




/*******************************************************************************
 Booleans
*******************************************************************************/


#define BOOL unsigned short int
#define FALSE 0
#define TRUE 1
