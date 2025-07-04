/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */


#ifndef _VECTRA_MACROS_
#define _VECTRA_MACROS_


#include <time.h>
//#include <assert.h>
#include <stdlib.h>


/** warning: numerical expressions are subject to syntax errors or
    operators precedence problems, unless the parameters are quoted in
    the body inside parenthesis */


/*******************************************************************************
| 
*******************************************************************************/

/*
  typedef struct vectra_tuple{ 
  float car; 
  float cdr; 
  } v_tuple;
*/

/*******************************************************************************
| ASSERT
*******************************************************************************/

/* 
   quand on compile avec:
   g++ -D NDEBUG 
   on définit la variable NDEBUG, et par conséquent tous les appels de assert() sont désactives:

   #ifndef NDEBUG
   #   define ASSERT(condition, message) \
   do { \
   if (! (condition)) { \
   std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
   << " line " << __LINE__ << ": " << message << std::endl; \
   std::exit(EXIT_FAILURE); \
   } \
   } while (false)
   #else
   #   define ASSERT(condition, message) do { } while (false)
   #endif

   On doit donc définir notre propre my_assert() si on veut que ca reste!
*/



#define MSG_ASSERT( CONDITION, MSG ) {					\
  std::cout.flush(); std::cerr.flush();					\
  if (! (CONDITION))							\
    {									\
      std::cerr << std::endl << "Dynamic assertion `" #CONDITION "` failed in " << __FILE__ \
		<< " line " << __LINE__ << ": <" << MSG << ">" << std::endl; \
      std::cerr.flush();						\
      exit( 1 );							\
    } \ 
}


#define ASSERT( CONDITION )			\
  MSG_ASSERT( CONDITION, " " )


/* if (!condition) */
/*     throw std::logic_error("Ach! Große malheur"); */


/*******************************************************************************
 Class Helper Macros
*******************************************************************************/


/* automatically define inline accessors in class definition */

#define READER(TYPE, VARNAME)			\
  TYPE get_##VARNAME() const			\
  { return this -> VARNAME; };

/* when the accessed variable is const */
#define CREADER(TYPE, VARNAME)			\
  const TYPE get_##VARNAME() const		\
  { return this -> VARNAME; };

#define WRITER(TYPE, VARNAME)			\
  void set_##VARNAME(const TYPE &value)		\
  { this -> VARNAME = value; };

#define READWRITER(TYPE, VARNAME)		\
  READER(TYPE, VARNAME)				\
  WRITER(TYPE, VARNAME) 


#define NREADER(TYPE, VARNAME, FUNNAME)		\
  TYPE FUNNAME() const				\
  { return this -> VARNAME; }			\


#define NWRITER(TYPE, VARNAME, FUNNAME)		\
  void FUNNAME(const TYPE &value)		\
  { this -> VARNAME = value; }			\

#define NREADWRITER(TYPE, VARNAME, SUFFIXNAME)	\
  NREADER(TYPE, VARNAME, get_##SUFFIXNAME)	\
  NWRITER(TYPE, VARNAME, set_##SUFFIXNAME)	\


/*******************************************************************************
 Language modifiers
*******************************************************************************/


#define SWITCH_STR_BEGIN(STR_VAR)		\
  {						\
  char* local_string = STR_VAR;			\
  if (STR_VAR == NULL)				\
    { fprintf(stderr, "\n null string read"); }	\


#define CASE_STR(STR_VAL)			\
  else if (!strcmp(local_string, STR_VAL))	\


#define SWITCH_STR_END				\
  else						\
    { fprintf(stderr, "unknown value"); }	\
  }						\


#define SWITCH_STR_ELSE_END( ELSE_CODE )	\
  else						\
    { ELSE_CODE }				\
  }						\


/*
  SWITCH_STR_BEGIN( argv[1] )
  CASE_STR("toto")
  { cout << "toto"; }
  CASE_STR("titi")
  { cout << "titi"; }
  SWITCH_STR_END;

  SWITCH_STR_BEGIN( argv[1] )
  CASE_STR("toto")
  { cout << "toto"; }
  SWITCH_STR_ELSE_END
  ( fprintf(stderr, "unknown value"); );
*/


/*******************************************************************************
AIR VOL HELPERS
*******************************************************************************/


// permet de charger en mémoire (row_major order) les données 
#define AIR_LOAD_LINEAR( _AIRVOL_PTR )		\
  _AIRVOL_PTR -> set_data_linear_mode(true);	\
  _AIRVOL_PTR -> update_definitions();		\
  _AIRVOL_PTR -> allocate();			\
  _AIRVOL_PTR -> update_data();


/*******************************************************************************
 CUDA HELPERS
*******************************************************************************/


#define CUDACALL( code )				\
  MSG_ASSERT( code == cudaSuccess , "failed CUDACALL");

// allocate data on the GPU memory, unpinned
#define CUDALLOC_GPU( _TAB, _DIM, _DATATYPE )			\
  MSG_ASSERT(							\
  cudaMalloc( (void**) &_TAB, _DIM * sizeof( _DATATYPE) ) \ 
== cudaSuccess , "failed CUDALLOC" );

// copy a CPU array to a GPU twin
#define CUDAPUSH( _TAB_CPU, _TAB_GPU, _DIM, _DATATYPE )			\
  MSG_ASSERT(								\
  cudaMemcpy( _TAB_GPU, _TAB_CPU, _DIM * sizeof( _DATATYPE), cudaMemcpyHostToDevice ) \ 
== cudaSuccess , "failed to push data to CUDA device");

// copy back a GPU array to its CPU twin
#define CUDAPULL( _TAB_CPU, _TAB_GPU, _DIM, _DATATYPE )			\
  MSG_ASSERT(								\
  cudaMemcpy( _TAB_CPU, _TAB_GPU, _DIM * sizeof( _DATATYPE), cudaMemcpyDeviceToHost ) \ 
== cudaSuccess , "failed to retrieve data from CUDA device"); 

#define CUDACOPY( _TAB_GPU_S, _TAB_GPU_D, _DIM, _DATATYPE )		\
  MSG_ASSERT(								\
  cudaMemcpy( _TAB_GPU_D, _TAB_GPU_S, _DIM * sizeof( _DATATYPE), cudaMemcpyDeviceToDevice ) \ 
== cudaSuccess , "failed to copy data inside CUDA device"); 



// allocate pinned memory in CPU memory, accelerates in/outbound traffic
// since the memory pinned will not be remapped
#define CUDALLOC_HP( _TAB, _DIM, _DATATYPE )				\
  MSG_ASSERT(								\
  cudaHostAlloc( (void**) &_TAB, _DIM * sizeof( _DATATYPE), cudaHostAllocDefault ) \ 
== cudaSuccess , "failed CUDALLOC for pinned mem" );				\


// allocate pinned memory in CPU memory for zero copy
// ie, implicit direct access from the GPU without CUDAPUSH
#define CUDALLOC_ZC( _TAB, _DIM, _DATATYPE )				\
  MSG_ASSERT(								\
  cudaHostAlloc( (void**) &_TAB, _DIM * sizeof( _DATATYPE), cudaHostAllocMapped ) \ 
== cudaSuccess , "failed " );				\


// when CUDALLOC_ZC was invoked on _TAB_CPU, this macro sets
// null-created _TAB_GPU pointer to a device-compatible value
#define CUDABIND_ZC( _TAB_CPU, _TAB_GPU )			\
  MSG_ASSERT(							\
	     cudaHostGetDevicePointer( &_TAB_GPU, _TAB_CPU, 0 )	\
	     == cudaSuccess , "failed binding for zero copy");	\
	
								
#define CUDABARRIER							\
  MSG_ASSERT(cudaThreadSynchronize() == cudaSuccess, "failed CUDA barrier");



/* affiche dans le flux spécifié (FILE*, stdin) la quantité de mémoire GPU disponible sur une carte CUDA */
#define CUDA_MEM_DISPLAY(__STREAM)					\
  {									\
    unsigned int free_mem, total_mem, used_mem;				\
    cuMemGetInfo( &free_mem, &total_mem );				\
    used_mem = total_mem-free_mem;					\
    fprintf(__STREAM, "\n total mem: %0.3f MB, free: %0.3f MB, used : %0.3f MB\n", \
	    ((double)total_mem)/1024.0/1024.0, ((double)free_mem )/1024.0/1024.0, \
	    ((double)used_mem )/1024.0/1024.0 );			\
  }									\


/*******************************************************************************
 Time Measurement - C++ - depend on Boost
*******************************************************************************/

/* don't forget to make: 
   #include "vChrono.h"
   before other includes

   link with -lboost_system -lboost_chrono
*/

// c'est que de la MERDE, passer à vChronos


#define vCHRONO_SET(NAME, MSG)				\
  vChrono<boost::chrono::system_clock> vC_##NAME;	\
  const char* vc_msg_##NAME = MSG;			\

#define vCHRONO_START(NAME)			\
  vC_##NAME.reset();				\
  cout << "\n::: start " << vc_msg_##NAME;	\
  cout.flush();					\

#define vCHRONO_STOP(NAME)				\
  cout << "\n::: stop " << vc_msg_##NAME << " in "	\
  << vC_##NAME.seconds() << "(s)";			\
  cout.flush();						\


// shorter:

#define vCHRONO_ASTART(NAME, MSG)			\
  vChrono<boost::chrono::system_clock> vC_##NAME;	\
  const char* vc_msg_##NAME = MSG;			\
  vC_##NAME.reset();					\
  cout << "\n::: start " << vc_msg_##NAME;		\
  cout.flush();						\


/* exemple d'utilisation
   vCHRONO_SET(cpu, "cpu call"); vCHRONO_START(cpu); 
   vecMoy1_h(tab_src_h, tab_count_h, tab_dst_h, N);
   vCHRONO_STOP(cpu);

   ou bien:
   vCHRONO_ASTART(cpu, "cpu call");
   // code to bench
   vCHRONO_STOP(cpu);

*/


// cela devient plus simple avec vChronos.h
#define TIMECODE(VCHRONOS, MSG, CODE_EXEC)	\
  {						\
    VCHRONOS.clear();				\
    VCHRONOS.change_message( MSG );		\
    VCHRONOS.start();				\
    CODE_EXEC;					\
    VCHRONOS.stop();				\
  }						\
    

/* #include <boost/thread.hpp>   */
/* #include <boost/date_time.hpp> */
// permet de paralléliser une attente avec l'exécution de code autre
#define CODE_EXEC_WAIT(VCHRONOS, T_MS, CODE_EXEC)	\
  {						\
    VCHRONOS.clear();				\
    VCHRONOS.set_mute();			\
    VCHRONOS.start();				\
    CODE_EXEC;					\
    VCHRONOS.stop();				\
    if ( VCHRONOS.get_elapsed_time_ms() < T_MS )	\
      {									\
    double remaining = double(T_MS) - VCHRONOS.get_elapsed_time_ms();	\
    boost::posix_time::milliseconds delay( remaining );			\
    boost::this_thread::sleep( delay );					\
      }									\
  }


/*******************************************************************************
 Math
*******************************************************************************/


/* REM: passer en inline */
/* tests if DOW <= X <= UP, since this expression is not evaluated by
   egcs as we could expect */
#define BETWEEN(X, DOW, UP)			\
  (((DOW) <= (X)) && ((X) <= (UP)))


/* tests if X is between range [mean - sd, mean + sd] */
#define BETWEEN_SD(X, MEAN, SD)				\
  BETWEEN((X), ((MEAN) - (SD)), ((MEAN) + (SD)))


/* */
#define add_limited(NUM, ADD, LIM)				\
  ((NUM) + (ADD) > (LIM)) ? (NUM) = (LIM) : (NUM) += (ADD); 


/* */
#define sub_limited(NUM, SUB, LIM)				\
  ((NUM) - (SUB) <= (LIM)) ? (NUM) = (LIM) : (NUM) -= (SUB); 


/* */
#define sub_positive(NUM, SUB)			\
  sub_limited(NUM, SUB, 0)


/* determines the square of a number (genericity needed) */
#define SQR(EXP)				\
  ((EXP) * (EXP))
/* NB: parenthesis enable use of a numerical expression without
   precedence operator problems */


/** notice expressions are of form ( (EXP1) )*/


/* determines the absolute value of a number */
/* works for expressions but with multiple evals */
#define ABS(EXP)				\
  (((EXP) > 0) ? (EXP) : -(EXP))


/* return the minimum value of both expressions values*/
#define MIN(EXP1, EXP2)				\
  (((EXP1) > (EXP2)) ? (EXP2) : (EXP1))


/* return the maximum value of both expressions values*/
#define MAX(EXP1, EXP2)				\
  (((EXP1) > (EXP2)) ? (EXP1) : (EXP2))

/*
  #if defined(_MSC_VER) 
  #  define MIN(a,b) std::_cpp_min(a,b)
  #  define MAX(a,b) std::_cpp_max(a,b)
  #else
  #  define MIN(a,b) std::min(a,b)
  #  define MAX(a,b) std::max(a,b)
  #endif
*/

/* ensure VAR stays within DOWN =< VAR =< UP, set to boundary value
   otherwise */
#define FRAME_VALUE(VAR, DOWN, UP)		\
  if (VAR < DOWN)				\
    VAR = DOWN;					\
  else if (VAR > UP)				\
    VAR = UP;


/* */
#define BOUNDEDP(VAR, DOWN, UP)			\
  ((VAR > DOWN) && (VAR < UP))

#define BOUNDED_EQ_P(VAR, DOWN, UP)		\
  ((VAR >= DOWN) && (VAR <= UP))


/*******************************************************************************
 Memory-allocation related: ARRAYS
*******************************************************************************/


#define ARRAY_ALLOC(TAB, DIM, DATATYPE)		       \ 
MSG_ASSERT(TAB = (DATATYPE *) calloc(DIM, sizeof(DATATYPE)),	\
	   "failed array alloc" ); 


#define ARRAY_DEC_ALLOC(TAB, DIM, DATATYPE)			\
  DATATYPE* TAB;						\
  MSG_ASSERT(TAB = (DATATYPE *) calloc(DIM, sizeof(DATATYPE)),	\
	     "array was declared but allocation failed"); 


/* présuppose
#include <xmmintrin.h>
//AVX
#include <immintrin.h>
//SSE4 (?)
#include <smmintrin.h>
// pour madvise
#include <sys/mman.h>
*/

#define ARRAY_SSE_ALLOC(TAB, DIM, ALIGN, DATATYPE)		       \ 
MSG_ASSERT(TAB = (DATATYPE *) _mm_malloc(DIM * sizeof(DATATYPE), ALIGN ), \
	   "failed sse array alloc" );					\
madvise(TAB, DIM * sizeof(DATATYPE), MADV_HUGEPAGE);			\ 



/*				     
#include <assert.h>
#define ARRAY_ALLOC(TAB, DIM, DATATYPE)		       \ 
assert( TAB = (DATATYPE *) calloc(DIM, sizeof(DATATYPE)) );
 
 #define ARRAY_DEC_ALLOC(TAB, DIM, DATATYPE)			\
 DATATYPE* TAB;						\
 assert(TAB = (DATATYPE *) calloc(DIM, sizeof(DATATYPE)));
*/


#define ARRAY_NEW(TAB, DIM, DATATYPE)		\
  try						\
  { int m_dim = (int) DIM;			\
    TAB = new DATATYPE[ m_dim ]; }		\
  catch(std::bad_alloc& ex)			\
  { cerr << ex.what() << endl;			\
    exit (EXIT_FAILURE); }			\


// préferer memset en c++
#define ARRAY_SET(TAB, DIM, VAL)		\
  for (int i = 0; i < DIM; i++)			\
    {						\
      TAB[DIM] = VAL;				\
    }
/* <=> memset( TAB, VAL, DIM * sizeof( le_type )) */


// préférer memcpy
#define ARRAY_CPY(TAB_SRC, TAB_DST, I, J)	\
  for (int n = I; n <= J; n++)			\
    {						\
      TAB_DST[n] = TAB_SRC[n];			\
    }

/* <=> memcpy( TAB_DST, TAB_SRC, n * sizeof( le_type )) */


#define ARRAY_REALLOC(TAB, NEWDIM, DATATYPE)				\
  MSG_ASSERT(TAB = (DATATYPE *) realloc(TAB, NEWDIM * sizeof(DATATYPE)), "reallocation failed" ); 


#define ARRAY_I_DISPLAY(array, size)		\
  {						\
    for (size_t i = 0; i < size; i++)		\
      printf("%d | ", array[i]);		\
  }



/*******************************************************************************
 n-arrays operations										
*******************************************************************************/


/* Reverses (destructively) the N first elements of given array */
/* déconne un brin */
#define ARRAY_REVERSE(TAB, N, SWP)			\
  {							\
    int I;						\
    int STOP = (int)N / 2;				\
    for (I = 0; I < STOP; I++)				\
      {REVERSE_VALUES(TAB[I], TAB[N - 1 - I], SWP);}	\
  }


/* initializes matrix to given value */
#define MATRIX_INIT(MAT, X, Y, VAL)		\
  {						\
    int T_I, T_J;				\
    for (T_I = 0; T_I < X; T_I++)		\
      for (T_J = 0; T_J < Y; T_J++)		\
	MAT[T_I][T_J] = VAL;			\
  }

/* copy matrix MAT1 to MAT2 */
#define MATRIX_COPY(MAT1, MAT2, X, Y)		\
  {						\
    int T_I, T_J;				\
    for (T_I = 0; T_I < X; T_I++)		\
      for (T_J = 0; T_J < Y; T_J++)		\
	MAT2[T_I][T_J] = MAT1[T_I][T_J];	\
  }


/* initializes matrix to given file content. TYPECODE = "f" for float,
   "d" for integer, "lf" for double, etc */
#define MATRIX_FILE_TXT(MAT, X, Y, FILENAME, TYPECODE)			\
  {									\
    FILE* fp;								\
    int T_I, T_J;							\
									\
    ASSERT(fp = fopen(FILENAME, "r"));					\
									\
    for (T_I = 0; T_I < X; T_I++)					\
      for (T_J = 0; T_J < Y; T_J++)					\
	ASSERT(fscanf(fp, "%" TYPECODE, &MAT[T_I][T_J]) != EOF);	\
									\
    fclose(fp);								\
  }


/* initializes matrix to given file content. TYPECODE = "f" for float,
   "d" for integer, etc */
#define MATRIX_FILE_RAW(MAT, X, Y, FILENAME, D_TYPE)			\
  {									\
    FILE* fp;								\
    int T_I, T_J;							\
									\
    ASSERT(fp = fopen(FILENAME, "rb"));					\
									\
    for (T_I = 0; T_I < X; T_I++)					\
      for (T_J = 0; T_J < Y; T_J++)					\
	ASSERT(fread(&MAT[T_I][T_J], sizeof(D_TYPE), 1, fp) != EOF);	\
									\
    fclose(fp);								\
  }


#define FILE_MATRIX_TXT(MAT, X, Y, FILENAME, TYPECODE)	\
  {							\
    FILE* fp;						\
    int T_I, T_J;					\
							\
    ASSERT(fp = fopen(FILENAME, "w"));			\
							\
    for (T_I = 0; T_I < X; T_I++)			\
      {							\
	for (T_J = 0; T_J < Y; T_J++)			\
	  fprintf(fp, "%" TYPECODE " ", MAT[T_I][T_J]); \
	fprintf(fp, "\n");				\
      }							\
							\
    fclose(fp);						\
  }



#define MATRIX_DISPLAY(MAT, X, Y, TYPECODE)		\
  {							\
    int T_I, T_J;					\
							\
    for (T_I = 0; T_I < X; T_I++)			\
      {							\
	printf("\n");					\
	for (T_J = 0; T_J < Y; T_J++)			\
	  printf("%" TYPECODE " ", MAT[T_I][T_J]);	\
      }							\
  }



/*******************************************************************************
 Connexity Iterators (img.h)
*******************************************************************************/


/* calls given function of arity 2 on every pixel coordinates around I, J 
   in 2D 4 connexity */
#define AROUND_4_CONNEX(I, J, FUN_OF_IJ)	\
  FUN_OF_IJ(I + 1, J);				\
  FUN_OF_IJ(I, J + 1);				\
  FUN_OF_IJ(I - 1, J);				\
  FUN_OF_IJ(I, J - 1);

/* 8 connexity without 4 connexity cases */
#define AROUND_8_ONLY_CONNEX(I, J, FUN_OF_IJ)	\
  FUN_OF_IJ(I + 1, J + 1);			\
  FUN_OF_IJ(I - 1, J + 1);			\
  FUN_OF_IJ(I - 1, J - 1);			\
  FUN_OF_IJ(I + 1, J - 1); 

/* calls given function of arity 2 on every pixel coordinates around I, J 
   in 2D 8 connexity. Order comply with Freeman code run. */
#define AROUND_8_FREEMAN_CONNEX(I, J, FUN_OF_IJ)	\
  FUN_OF_IJ(I + 1, J);					\
  FUN_OF_IJ(I + 1, J + 1);				\
  FUN_OF_IJ(I, J + 1);					\
  FUN_OF_IJ(I - 1, J + 1);				\
  FUN_OF_IJ(I - 1, J);					\
  FUN_OF_IJ(I - 1, J - 1);				\
  FUN_OF_IJ(I, J - 1);					\
  FUN_OF_IJ(I + 1, J - 1); 

/* */
#define AROUND_8_CONNEX(I, J, FUN_OF_IJ)	\
  AROUND_4_CONNEX(I, J, FUN_OF_IJ)		\
  AROUND_8_ONLY_CONNEX(I, J, FUN_OF_IJ)


/* the same , but if local boolean WITH_HALT is set to
   TRUE, iteration stops. */
#define AROUND_COND_8_CONNEX(I, J, FUN_OF_IJ)		\
  WITH_BREAK({						\
      {							\
	BOOL WITH_HALT = FALSE;				\
	FUN_OF_IJ(I + 1, J);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I, J + 1);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I - 1, J);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I, J - 1);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I + 1, J + 1); if (WITH_HALT) break;	\
	FUN_OF_IJ(I - 1, J + 1); if (WITH_HALT) break;	\
	FUN_OF_IJ(I - 1, J - 1); if (WITH_HALT) break;	\
	FUN_OF_IJ(I + 1, J - 1); if (WITH_HALT) break;	\
      }});


/* the same in Freeman order, but if local boolean WITH_HALT is set to
   TRUE, iteration stops. */
#define AROUND_COND_8_FREEMAN_CONNEX(I, J, FUN_OF_IJ)	\
  WITH_BREAK({						\
      {							\
	BOOL WITH_HALT = FALSE;				\
	FUN_OF_IJ(I + 1, J);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I + 1, J + 1); if (WITH_HALT) break;	\
	FUN_OF_IJ(I, J + 1);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I - 1, J + 1); if (WITH_HALT) break;	\
	FUN_OF_IJ(I - 1, J);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I - 1, J - 1); if (WITH_HALT) break;	\
	FUN_OF_IJ(I, J - 1);     if (WITH_HALT) break;	\
	FUN_OF_IJ(I + 1, J - 1); if (WITH_HALT) break;	\
      }});


/* calls given function of arity 3 on every voxel coordinates around I, J, K 
   in 3D 6 connexity */
#define AROUND_6_CONNEX(I, J, K, FUN_OF_IJK)	\
  FUN_OF_IJK(I + 1, J, K);			\
  FUN_OF_IJK(I - 1, J, K);			\
  FUN_OF_IJK(I, J + 1, K);			\
  FUN_OF_IJK(I, J - 1, K);			\
  FUN_OF_IJK(I, J, K + 1);			\
  FUN_OF_IJK(I, J, K - 1);


/* calls given function of arity 3 on every voxel coordinates around I, J, K 
   in 3D 6 connexity */
#define AROUND_26_CONNEX(I, J, K, FUN_OF_IJK)	\
  /* first K slide around voxel */		\
  FUN_OF_IJK(I + 1, J, K);			\
  FUN_OF_IJK(I - 1, J, K);			\
  FUN_OF_IJK(I + 1, J + 1, K);			\
  FUN_OF_IJK(I + 1, J - 1, K);			\
  FUN_OF_IJK(I - 1, J + 1, K);			\
  FUN_OF_IJK(I - 1, J - 1, K);			\
  FUN_OF_IJK(I, J + 1, K);			\
  FUN_OF_IJK(I, J - 1, K);			\
  /* the same at K+1 (one extra voxel) */	\
  FUN_OF_IJK(I + 1, J, K + 1);			\
  FUN_OF_IJK(I - 1, J, K + 1);			\
  FUN_OF_IJK(I + 1, J + 1, K + 1);		\
  FUN_OF_IJK(I + 1, J - 1, K + 1);		\
  FUN_OF_IJK(I - 1, J + 1, K + 1);		\
  FUN_OF_IJK(I - 1, J - 1, K + 1);		\
  FUN_OF_IJK(I, J + 1, K + 1);			\
  FUN_OF_IJK(I, J - 1, K + 1);			\
  FUN_OF_IJK(I, J, K + 1);			\
  /* the same at K-1 */				\
  FUN_OF_IJK(I + 1, J, K - 1);			\
  FUN_OF_IJK(I - 1, J, K - 1);			\
  FUN_OF_IJK(I + 1, J + 1, K - 1);		\
  FUN_OF_IJK(I + 1, J - 1, K - 1);		\
  FUN_OF_IJK(I - 1, J + 1, K - 1);		\
  FUN_OF_IJK(I - 1, J - 1, K - 1);		\
  FUN_OF_IJK(I, J + 1, K - 1);			\
  FUN_OF_IJK(I, J - 1, K - 1);			\
  FUN_OF_IJK(I, J, K - 1); 


/*******************************************************************************
 Generic 
*******************************************************************************/



/* execute given code unless condition is true */
#define UNLESS(COND, CODE_EXEC)			\
  if (! (COND))					\
  CODE_EXEC
/* usage:
   UNLESS(t == -1),
   {
   blah
   })*/


/* reverse the values of A and B 
   requires extra variable of same type for value swapping */
#define REVERSE_VALUES(A, B, SWP)		\
  SWP = B; B = A; A = SWP;


/* display an integer in a lazy way */
#define DISP_INT(NUM)				\
  printf("\n%d", NUM);



// en C++, préférer les template functions (qui sont également des macros)
#define DEFUN_MAX_FUNCTION(TYPE)		\
  TYPE						\
  max_##TYPE (TYPE* array, int start, int card) \
  {						\
    int i, end = start + card;			\
    TYPE max = array[start++];			\
						\
    for (i = start; i < end; i++)		\
      {						\
	if (array[i] > max) max = array[i];	\
      }						\
						\
    return max;					\
  }
/* declares:
   TYPE max_TYPE(TYPE* array, int start, int card);
   
   which returns maximum value in array from a range of card elements
   starting from given index */


/*******************************************************************************
ASM
*******************************************************************************/

//volatile: interdire au compilo d'optimiser l'intérieur

#define ASM_INSERT( CODE, FOOTER ) {			\
    __asm__ 	__volatile__				\
      (							\
       ".intel_syntax noprefix \n\t"			\
       CODE						\
       ".att_syntax noprefix  \n\t"			\
       FOOTER						\
						);	\
  }							\


// #define ASM2(code) code
#define ASM( LINESTART, LINEEND ) " " #LINESTART "," #LINEEND "\n\t"


/* example */
/*
  int out = 1, in1 = 4, in2 = 5;

  //#define MYFOOTER				\
  // output    \
  :"=r"(out)   \    
  // input     \ 
  :"r"(in1), "r"(in2) \
	
  ASM_INSERT(
  ASM(mov eax, %1)
  ASM(mov ebx, %2)
  ASM(add eax, ebx)
  ASM(mov %0, eax)
  ,
  MYFOOTER
  );

*/
 
/*******************************************************************************
ASM END
*******************************************************************************/


#include "legacy_macros.h"


#endif /* _VECTRA_MACROS_ */



/*
#define xstr(s) str(s)
     #define str(s) #s
     #define foo 4
     str (foo)
          ==> "foo"
     xstr (foo)
          ==> xstr (4)
          ==> str (4)
          ==> "4"
*/
