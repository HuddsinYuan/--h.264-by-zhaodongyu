#include "image.h"
/* *************************************************************** */
/* pack(int option, long data, int bits, FILE * output)

   Basic Usage :  file_ptr = fopen("filename", "wb+");
                  pack(0,data,bits,file_ptr); //(no initialization)
                  pack(0,data,bits,file_ptr);
                  ...
                  pack(2,0,0,file_ptr); //to retreive file size
                  pack(1,0,0,file_ptr); //to finalize file-write,
                                        //this must happen before
                                        //switching to another
                                        //file or before closing
                                        //it with fclose (empty
                                        //buffer and reset
                                        //pointer)
                  fclose(file_ptr);
   
                  **For multiple tests, reset files with :
                  pack(4,0,0,file_ptr);
   Return value : (long)number of bits successfully writen (0,1),
                  the incremental file size in bytes (2) or the   
                  number of bits deleted from buffer (3)
/* *************************************************************** */

long pack(option, data, bits, output)
int option;  /* option number (0-3):                                  */
             /* 0 - Write to file stream (usual)                      */ 
             /* 1 - Finalize file-write before closing or             */ 
             /*     switching to another file (empty buffer           */    
             /*     and reset pointer)                                */
             /* 2 - Get actual size of file in (Bytes)                */
             /* 3 - Delete bits : Delete last few bits written        */
             /*     in the buffer, if the buffer has just been        */
             /*     written, then a Warning will appear saying        */
             /*     it's too late.                                    */
             /* 4 - Reset file (clear and reset pointers and all)     */
long data;   /* data to be written                                    */
int bits;    /* number of bits to be written to buffer. For option 3, */
		     /* 'bits' is the number of bits deleted                  */ 
FILE * output; /* File stream pointer */    
{
  int i;
  long value=0;
  static int buffer_ptr = 0; /* buffer position pointer */
  static unsigned char buffer = 0; /* 8-bit buffer */

  switch(option)
  {
    case 0 : /* Write */
		   for(i=0;i<bits;i++)
		   {
             buffer = buffer<<1;
             if(data & 1)
                buffer |= 1;
             data = data>>1;
             buffer_ptr++;
             if(buffer_ptr==8)
			 {
               fwrite(&buffer,sizeof(unsigned char),1,output);
               buffer_ptr=0;
               buffer=0;
			 }
             value++;
		   }             
           break;

    case 1 : /* Finalize */

           if (buffer_ptr != 0)
		   {
             value = 8-buffer_ptr;
             buffer = buffer<<value;
             fwrite(&buffer,sizeof(unsigned char),1,output);
             buffer_ptr=0;
             buffer=0;
		   }
           fflush(output);
           break;

    case 2 : /* Get incremental size (Bytes) */
    
 //          fflush(output);
           value = (long) (ftell(output) + ceil((double)buffer_ptr/8));
           break;

    case 3 : /* Delete bits (must be in update mode - wb+) */
      
           if(bits<buffer_ptr)
		   {
             buffer = buffer>>bits;
             buffer_ptr-=bits;
             value = bits;
		   }
		   else /* A venir : Vider le buffer, reculer le pointeur de 4 octets faire fread (mode wb+) */
		   {
             printf("Warning, Buffer does not contain %d bits, will only delete %d bits.\n", bits, buffer_ptr);
             value = buffer_ptr;
		   }      
           break;

    case 4 : /* Reset file */
      
           rewind(output);
           buffer_ptr = 0;
           buffer = 0;
           value = 0;
  }
  return value;
}


/* ************************************************************** */
/* unpack (FILE * output, int bits)
   
  Basic Usage :   file_ptr = fopen("filename","rb");
                  unpack(bits,file_ptr);//(no initialization)
                  unpack(bits,file_ptr);
                  unpack(0,file_ptr);//empty buffer and reset
                                     //before going to next file
                  ...
                  fclose(file_ptr);//It's that simple!
  
  Return value : (long) data read from the file and then the buffer
/* ************************************************************** */

long unpack(bits, output)
int bits; /* number of bits to be read from buffer */
FILE * output; /* File stream pointer */
{
  int i;
  long data = 0;
  static int buffer_ptr = 0; /* buffer position pointer */
  static unsigned char buffer = 0; /* unpacked bits */
  
  if(bits==0) /* Don't unpack anything just reset buffer and pointer */
  {
    buffer_ptr=0;
    buffer=0;
    return 0;
  }

  for(i=0;i<bits;i++)
  {
    if(buffer_ptr==0) /* make sure buffer is never empty */
	{
      fread(&buffer,sizeof(unsigned char),1,output);
      buffer_ptr=8;
	}
    
	if(buffer & (1<<7))
       data |= 1<<i;
    
	buffer = buffer<<1;
    buffer_ptr--; 
  }
  return data;
}
