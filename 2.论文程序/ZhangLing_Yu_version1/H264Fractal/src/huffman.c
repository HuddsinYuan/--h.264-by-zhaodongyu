#include "huffman.h"
#include "global.h"

extern long pack(int option, long data, int bits, FILE * output); 
void CreateHuffmanCodeBook(Book,N,Prob)
HCODE*Book; /* Huffman codebook pointer */
int N; /* number of symbols */
long * Prob; /* input probability array */
{
  int i,k;
  int s1,s2; /* active symbols */
  long p1,p2; /* active probabilities */
  int *merge1, *merge2;

  merge1 = (int*) malloc((N-1)*sizeof(int));
  merge2 = (int*) malloc((N-1)*sizeof(int));

  /* Set Active flags using Book->code buffer */

  for(i=N-1;i>=0;i--) Book[i].code = 1; /* all accessible for merging */

  /* Loop of merging symbols */

  for(k=N-2;k>=0;k--)  /* merge k & k+1 symbols to k */
  {
	/* choose s1 and s2 that have the smallest prob p1 and p2 */
	p1 = (1<<30);
	p2 = (1<<30);
	s1 = 0;
	s2 = 0;
	for(i=N-1;i>=0;i--)
	{
	  if(Book[i].code == 1) /* means that the code is still accessible for merging */
	  { 
	    if(Prob[i]<=p2)
		{
		  if(Prob[i]<=p1)
		  { //(p1 < p2)
		    p2 = p1;
		    s2 = s1;
		    s1 = i;
		    p1 = Prob[s1];
		  }
		  else
		  {
		    s2 = i;
		    p2 = Prob[s2];
		  } /* At the end, p1 and p2 have the smallest probs */
		} /* The 0-prob symboles at the end of the list are the first merged */ 
	  }
	}
    /* merge合并 s2 to s1 */
    merge1[k] = s1; /* Store the indices of the smallest probabilities for each iteration (merge) */
    merge2[k] = s2;
//     printf("Merger %d:(%d,%d)\n",k,merge1[k],merge2[k]);
    Prob[s1] += Prob[s2];/* p1 becomes the sum and p2 is no longer available for merging */
    Book[s2].code = 0; 
  }

  /* Set initial code for the last remaining symbol */
  Book[s1].code = Book[s1].bits = 0; /* s1 is the Root of the tree at the end */

  /* Loop for creating codes, starting from the end (Root)            */
  /* We go from right to left and create the tree (no code, bits = 0) */
  for(k=0;k<N-1;k++) 
  {
	s1 = merge1[k]; /* We start with the strongest probs */
	s2 = merge2[k];
	if(Book[s1].bits<MAX_CL)
	{
	  Book[s2].code = Book[s1].code | (1<<Book[s1].bits); /* 1 is for the strongest prob */
	  Book[s2].bits = (++Book[s1].bits); /* 0 for the smallest and the code doesnt change. But, we add one bit each. */
	  /* Remark: if the number of bits excedes the maximum, the symbol must be coded         */
	  /* using one of the longuest codes. (We must send the symbol with that longuest code.) */
	}
	else
	{ /* at MAX_CL bits, the code is not extended so it stays the same */
	  Book[s2].code = Book[s1].code;
	  Book[s2].bits = Book[s1].bits; /* number of bits already at MAX_CL */
	}
  }
  free(merge1);
  merge1= NULL;
  free(merge2);
  merge2= NULL;
}


void HuffmanEncoder(Code,Data_length, Data, Book, Code_length, symbol_length)
unsigned char *Code; /* output code array */
int Data_length; /* input data length */
int * Data;	/* input data buffer */
HCODE * Book; /* Huffman code book */
long  *Code_length; /* output code length in bits */
int symbol_length; /* symbol length needed for extended code book */
{
  //long k,c,i;
  long k,c;
  unsigned char *codptr; /* output code cursor */
  unsigned char bitoff;	/* bit offset cursor */

  codptr = Code;
  bitoff = 0;
  for(k=0;k<Data_length;k++)//5338
  {
	/* Now would be a good time to include the symbol of the extended code book before encoding it */
	if(Book[Data[k]].bits == MAX_CL)
	{
	  c = (Data[k]<<(MAX_CL+bitoff))|(Book[Data[k]].code<<bitoff); /* symbol the value */
	  bitoff += (MAX_CL + symbol_length);
	  //printf("Using extended codebook with c = %ld\n",c);
	}
	else
	{
	  c = (Book[Data[k]].code<<bitoff); /* symbol */
	  bitoff += Book[Data[k]].bits;
	}

	codptr[0] |= (c&255);
	codptr[1] = ((c>>8)&255);
	codptr[2] = ((c>>16)&255);
	codptr[3] = ((c>>24)&255);

	while(bitoff >=8)
	{
	  bitoff -= 8;
	  codptr++;
	}
  }
  *Code_length = ((codptr - Code)<<3) + bitoff;
}


void GetInverseCodeBook(N, Book, InvBook)
int N; /* Number of symbols */
HCODE * Book; /* Huffman code book */
int * InvBook;
{
  long k,b,c;
        
  for(k=N-1;k>=0;k--) /* get all N codes and create an inverse code for each symbol */
  {
	b=(1<<Book[k].bits); /* position a '1' ahead of the code to be read */
	for(c=Book[k].code;c<=MASK;c+=b)
		InvBook[c]=k;
	/* The InvBook index begins with the actual code, then it is increased by     */
    /* positioning a '1' in front of it, then '10' then '11' and so forth until   */
	/* the maximum number of bits for a code is reached. All of these codes c now */
	/* correspond to the symbole k so we can decode the entire word instanta-     */
	/* neously. Only one code will be recognizable as the longuest code for which */
	/* we must read the symbol sent with it.                                      */ 
  }
}


void HuffmanDecoder(Data,Book, InvBook, Code, Data_length, symbol_length)
int * Data; /* output decoded data buffer */
HCODE * Book; /* Huffman code book */
int * InvBook; /* Inverse code book */
unsigned char * Code; /* input code buffer */
long Data_length; /* decoded data length (number of samples) */
int symbol_length; /* symbol length needed for extended code book */
{
  long k,c;
  unsigned char *codptr; /* reading code cursor */
  unsigned char bitoff;	/* bit offset cursor */	 

  codptr = Code;
  bitoff = 0;
  for(k=0;k<Data_length;k++) /* number of samples to decode */
  {
	c=codptr[0]|(codptr[1]<<8)|(codptr[2]<<16)|(codptr[3]<<24); /* Retreive a long */

	if(Book[InvBook[(c>>bitoff)&MASK]].bits == MAX_CL)
	{
	  Data[k]=(c>>(MAX_CL+bitoff))&((1<<symbol_length)-1);
	  bitoff += (MAX_CL + symbol_length);
	}
	else
	  bitoff += Book[Data[k]=InvBook[(c>>bitoff)&MASK]].bits;

	while(bitoff >= 8)
	{
	  bitoff -= 8;
	  codptr ++;
	}
  }
}
// //huffman//////ZZZZZZZZZZZZZZZZZZZ
// void huffman_Code()
// {
// 	int i,j;
// 	for (i=0;i<num_regions;i++)
// 	{	
// 		for(j=0; j<(1<<MAX_BITS); j++)//10
// 		{
//             Huff_ptr[i][0][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][0][j].bits=0;
// 
// 			Huff_ptr[i][1][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][1][j].bits=0;
// 
// 			Huff_ptr[i][2][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][2][j].bits=0;
// 
// 			Huff_ptr[i][3][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][3][j].bits=0;
// 
// 			Huff_ptr[i][4][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][4][j].bits=0;
// 		
// 			Huff_ptr[i][5][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][5][j].bits=0;
// 		
// 			Huff_ptr[i][6][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][6][j].bits=0;
// 		
// 			Huff_ptr[i][7][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][7][j].bits=0;
// 			
// 			Huff_ptr[i][8][j].code=0; /* i: region; k: parameters; j: bits */ 
//             Huff_ptr[i][8][j].bits=0;
// 
// 		}
// 		CreateHuffmanCodeBook(Huff_ptr[i][0], 1<<huff_search_range, xy_p[i]);//前面是算出来的。中间是16
// 		CreateHuffmanCodeBook(Huff_ptr[i][1], (int)((MAX_ALPHA-MIN_ALPHA)*100)/5+1, s_p[i]);//中间是128
// 		CreateHuffmanCodeBook(Huff_ptr[i][2], (MAX_BETA-MIN_BETA)/5+1, o_p[i]);//中间是64
//         //Y
// 		CreateHuffmanCodeBook(Huff_ptr[i][3], 10, redual_ple[i]);//N_ple
// 		CreateHuffmanCodeBook(Huff_ptr[i][4], 10, redual_pru[i]);//N_pru
// 
// 		//UV
// 		CreateHuffmanCodeBook(Huff_ptr[i][5], 10, redual_ple_u[i]);//N_ple_u
// 		CreateHuffmanCodeBook(Huff_ptr[i][6], 10, redual_pru_u[i]);//N_pru_u
// 		CreateHuffmanCodeBook(Huff_ptr[i][7], 10, redual_ple_v[i]);//N_ple_v
// 		CreateHuffmanCodeBook(Huff_ptr[i][8], 10, redual_pru_v[i]);//N_pru_v
// 
// 		HuffmanEncoder(code[0][i], trans_count[i], x_trans[i], Huff_ptr[i][0],&code_length[0][i], huff_search_range);
// 		HuffmanEncoder(code[1][i], trans_count[i], y_trans[i], Huff_ptr[i][0],&code_length[1][i], huff_search_range);
// 		HuffmanEncoder(code[2][i], trans_count[i], s_trans[i], Huff_ptr[i][1],&code_length[2][i], (int)ceil(log((int)((MAX_ALPHA-MIN_ALPHA)*100)/5+1)/log(2)));
// 		HuffmanEncoder(code[3][i], trans_count[i], o_trans[i], Huff_ptr[i][2],&code_length[3][i], (int)ceil(log((MAX_BETA-MIN_BETA)/5+1)/log(2)));
// 
// 		//Y
// 		HuffmanEncoder(code[4][i], trans_count_re[i], redual_transle[i], Huff_ptr[i][3],&code_length[4][i], 10);//N_ple
// 		HuffmanEncoder(code[5][i], trans_count_re[i], redual_transru[i], Huff_ptr[i][4],&code_length[5][i], 10);//N_pru
// 
// 
// 		//UV
// 		HuffmanEncoder(code[6][i], trans_count_re_u[i], redual_transle_u[i], Huff_ptr[i][5],&code_length[6][i], 10);//N_ple_u
// 		HuffmanEncoder(code[7][i], trans_count_re_u[i], redual_transru_u[i], Huff_ptr[i][6],&code_length[7][i], 10);//N_pru_u
// 		HuffmanEncoder(code[8][i], trans_count_re_v[i], redual_transle_v[i], Huff_ptr[i][7],&code_length[8][i], 10);//N_ple_v
// 		HuffmanEncoder(code[9][i], trans_count_re_v[i], redual_transru_v[i], Huff_ptr[i][8],&code_length[9][i], 10);//N_pru_v
// 
// 
// 	}
// 	
// }