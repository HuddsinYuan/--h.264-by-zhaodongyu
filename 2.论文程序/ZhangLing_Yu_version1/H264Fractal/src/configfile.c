
#define INCLUDED_BY_CONFIGFILE_C

#include "configfile.h"

#define MAX_ITEMS_TO_PARSE  10000
extern InputParameters *input;

#define strcasecmp strcmpi
void PatchInp ();

int TestEncoderParams(int bitdepth_qp_scale)
{
	int i = 0;
	
	while (Map[i].TokenName != NULL)
	{
		if (Map[i].param_limits == 1)
		{
			if (Map[i].Type == 0)
			{
				if ( * (int *) (Map[i].Place) < (int) Map[i].min_limit || * (int *) (Map[i].Place) > (int) Map[i].max_limit )
				{
					printf("Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, (int) Map[i].min_limit,(int)Map[i].max_limit );
					exit(0);
				}
				
			}
			else if (Map[i].Type == 2)
			{
				if ( * (double *) (Map[i].Place) < Map[i].min_limit || * (double *) (Map[i].Place) > Map[i].max_limit )
				{
					printf("Error in input parameter %s. Check configuration file. Value should be in [%.2f, %.2f] range.", Map[i].TokenName,Map[i].min_limit ,Map[i].max_limit );
					exit(0);
				}        
			}            
		}
		else if (Map[i].param_limits == 2)
		{
			if (Map[i].Type == 0)
			{
				if ( * (int *) (Map[i].Place) < (int) Map[i].min_limit )
				{
					printf("Error in input parameter %s. Check configuration file. Value should not be smaller than %d.", Map[i].TokenName, (int) Map[i].min_limit);
					exit(0);
				}
				
			}
			else if (Map[i].Type == 2)
			{
				if ( * (double *) (Map[i].Place) < Map[i].min_limit )
				{
					printf("Error in input parameter %s. Check configuration file. Value should not be smaller than %2.f.", Map[i].TokenName,Map[i].min_limit);
					exit(0);
				}        
			}
		}
		else if (Map[i].param_limits == 3) // Only used for QPs
		{
			if (Map[i].Type == 0)
			{
				if ( * (int *) (Map[i].Place) !=  (int ) (Map[i].min_limit) && * (int *) (Map[i].Place) !=  (int ) (Map[i].max_limit) )
				{
					printf("Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, (int) (Map[i].min_limit - bitdepth_qp_scale),(int)Map[i].max_limit );
					exit(0);
				}
				
			}
		}
		
		i++;
	}
	return -1;
}

int ParameterNameToMapIndex (char *s)
{
	int i = 0;
	
	while (Map[i].TokenName != NULL)
		if (0==strcasecmp (Map[i].TokenName, s))
			return i;
		else
			i++;
		return -1;
}

int DisplayEncoderParams()
{
	int i = 0;
	
	printf("******************************************************\n");
	printf("*               Encoder Parameters                   *\n");
	printf("******************************************************\n");
	while (Map[i].TokenName != NULL)
	{
		if (Map[i].Type == 0)
			printf("Parameter %s = %d\n",Map[i].TokenName,* (int *) (Map[i].Place));
		else if (Map[i].Type == 1)
			printf("Parameter %s = ""%s""\n",Map[i].TokenName,(char *)  (Map[i].Place));
		else if (Map[i].Type == 2)
			printf("Parameter %s = %.2f\n",Map[i].TokenName,* (double *) (Map[i].Place));
		i++;
	}
	printf("******************************************************\n");
	return -1;
}

int InitEncoderParams()
{
	int i = 0;
	
	while (Map[i].TokenName != NULL)
	{
		if (Map[i].Type == 0)
			* (int *) (Map[i].Place) = (int) Map[i].Default;
		else if (Map[i].Type == 2)
			* (double *) (Map[i].Place) = Map[i].Default;
		i++;
	}
	return -1;
}

char *GetConfigFileContent (char *Filename)
{
	long FileSize;
	FILE *f;
	char *buf;
	
	if (NULL == (f = fopen (Filename, "r")))
	{
		printf ("Cannot open configuration file %s.", Filename);
		return NULL;
	}
	
	if (0 != fseek (f, 0, SEEK_END))
	{
		printf ("Cannot fseek in configuration file %s.", Filename);
		return NULL;
	}
	
	FileSize = ftell (f);
	if (FileSize < 0 || FileSize > 60000)
	{
		printf ("Unreasonable Filesize %ld reported by ftell for configuration file %s.", FileSize, Filename);
		return NULL;
	}
	if (0 != fseek (f, 0, SEEK_SET))
	{
		printf ("Cannot fseek in configuration file %s.", Filename);
		return NULL;
	}
	
	if ((buf = malloc (FileSize + 1))==NULL) 
		printf("Could not allocate memory: %s",Filename);
	
	// Note that ftell() gives us the file size as the file system sees it.  The actual file size,
	// as reported by fread() below will be often smaller due to CR/LF to CR conversion and/or
	// control characters after the dos EOF marker in the file.
	
	FileSize = fread (buf, 1, FileSize, f);
	buf[FileSize] = '\0';
	
	
	fclose (f);
	return buf;
}

void ParseContent (char *buf, int bufsize)
{
	char *items[MAX_ITEMS_TO_PARSE];
	int MapIdx;
	int item = 0;
	int InString = 0, InItem = 0;
	char *p = buf;
	char *bufend = &buf[bufsize];
	int IntContent;
	double DoubleContent;
	int i;
	
	// Stage one: Generate an argc/argv-type list in items[], without comments and whitespace.
	// This is context insensitive and could be done most easily with lex(1).
	
	while (p < bufend)
	{
		switch (*p)
		{
		case 13:
			p++;
			break;
		case '#':                 // Found comment
			*p = '\0';              // Replace '#' with '\0' in case of comment immediately following integer or string
			while (*p != '\n' && p < bufend)  // Skip till EOL or EOF, whichever comes first
				p++;
			InString = 0;
			InItem = 0;
			break;
		case '\n':
			InItem = 0;
			InString = 0;
			*p++='\0';
			break;
		case ' ':
		case '\t':              // Skip whitespace, leave state unchanged
			if (InString)
				p++;
			else
			{                     // Terminate non-strings once whitespace is found
				*p++ = '\0';
				InItem = 0;
			}
			break;
			
		case '"':               // Begin/End of String
			*p++ = '\0';
			if (!InString)
			{
				items[item++] = p;
				InItem = ~InItem;
			}
			else
				InItem = 0;
			InString = ~InString; // Toggle
			break;
			
		default:
			if (!InItem)
			{
				items[item++] = p;
				InItem = ~InItem;
			}
			p++;
		}
	}
	
	item--;

    for (i=0; i<item; i+= 3)
	{
		if (0 > (MapIdx = ParameterNameToMapIndex (items[i])))
		{
			printf (" Parsing error in config file: Parameter Name '%s' not recognized.", items[i]);
			exit(0);
		}
		if (strcasecmp ("=", items[i+1]))
		{
			printf (" Parsing error in config file: '=' expected as the second token in each line.");
			exit(0);
		}

		// Now interpret the Value, context sensitive...
		
		switch (Map[MapIdx].Type)
		{
		case 0:           // Numerical
			if (1 != sscanf (items[i+2], "%d", &IntContent))
			{
				printf (" Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
				exit(0);
			}
			* (int *) (Map[MapIdx].Place) = IntContent;
			printf (".");
			break;
		case 1:
			strncpy ((char *) Map[MapIdx].Place, items [i+2], FILE_NAME_SIZE);
			printf (".");
			break;
		case 2:           // Numerical double
			if (1 != sscanf (items[i+2], "%lf", &DoubleContent))
			{
				printf (" Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
				exit(0);
			}
			* (double *) (Map[MapIdx].Place) = DoubleContent;
			printf (".");
			break;
		default:
			printf("Unknown value type in the map definition of configfile.h\n");
		}
	}
	memcpy (input, &configinput, sizeof (InputParameters));
}

void Configure(int ac, char *av[])
{
	char *content;    
	char *filename=DEFAULTCONFIGFILENAME;
	if (0 == strncmp (av[1], "-d", 2))
	{
		filename=av[2];		
	}
/////原来分形没有定义 _DEBUG, 但是却当作是有的，所以我在这里定义了。
#ifndef _DEBUG
#define _DEBUG
#endif
//////////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
   
#else
	
	int str1=strlen(*av);

    char *pdest;
    pdest=strstr(*av,"Release");

	int str2=strlen(strstr(*av,"Release"));
    int str12=str1-str2;
	strncpy(filename,*av,strlen(*av)-strlen(strstr(*av,"Release")));
	filename[strlen(*av)-strlen(strstr(*av,"Release"))] = '\0';
#endif
   // strcat(filename,"encoder.cfg");
    memset (&configinput, 0, sizeof (InputParameters));
//     printf ("Setting Default Parameters...\n");ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	InitEncoderParams();		

// 	printf ("Parsing Configfile %s\n", filename);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	content = GetConfigFileContent (filename);	
	if (NULL == content)
	{
		printf("Error: No Encoder.cfg file!\n");
		exit(0);
	}

    ParseContent (content, strlen(content));
	printf ("\n");
    free (content);
    printf ("\n");
	PatchInp();
	if (input->displayEncoderParams)
		DisplayEncoderParams();
}

/*!
************************************************************************
* \brief
*    calculate Ceil(Log2(uiVal))
************************************************************************
*/
unsigned CeilLog2( unsigned uiVal)
{
	unsigned uiTmp = uiVal-1;
	unsigned uiRet = 0;
	
	while( uiTmp != 0 )
	{
		uiTmp >>= 1;
		uiRet++;
	}
	return uiRet;
}


void PatchInp()
{
	int i;
	int bitdepth_qp_scale = 6*(input->bitDepthLuma - 8);
	TestEncoderParams(bitdepth_qp_scale);
	if (input->frame_rate == 0.0)
		input->frame_rate = INIT_FRAME_RATE;



///+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
// 
// 	// Open Files
// 	if ((p_in=fopen(input->infile,"rb"))==NULL)
// 	{
// 		snprintf(errortext, ET_SIZE, "Input file %s does not exist",input->infile);
// 		error (errortext, 500);
// 	}
	// Open Files
// 	if ((p_in=fopen(input->infile,"rb"))==NULL)
// 	{
// 		snprintf(errortext, ET_SIZE, "Input file %s does not exist",input->infile);
// 		error (errortext, 500);
// 	}
	
	if (strlen (input->ReconFile) > 0 && (p_dec=fopen(input->ReconFile, "wb"))==NULL)
	{
		printf(errortext, ET_SIZE, "Error open file %s", input->ReconFile);
		error (errortext, 500);
	}
	
	if (strlen (input->TraceFile) > 0 && (p_trace=fopen(input->TraceFile,"w"))==NULL)
	{
		printf(errortext, ET_SIZE, "Error open file %s", input->TraceFile);
		error (errortext, 500);
	}



///+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////

////////////////////////////////////////////////////////////////////////////////////////////
	// Open Files
	/////////////////////////////////////////////R encode file
	if ((fp_in_c=fopen(input->infile_c,"rb"))==NULL)
	{
		printf("Input file_C %s does not exist\n",input->infile_c);
		exit(0);
	}
	if ((fp_in_r=fopen(input->infile_r,"rb"))==NULL)
	{
		printf("Input file_R %s does not exist\n",input->infile_r);
		exit(0);
	}
	if ((fp_in_l=fopen(input->infile_l,"rb"))==NULL)
	{
		printf("Input file_L %s does not exist\n",input->infile_l);
		exit(0);
	}
	if ((fp_out_all_c_rec=fopen(input->outall_c,"wb"))==NULL)
	{
		printf("Output_all_c %s does not exist\n",input->outall_c);
		exit(0);
	}
	if ((fp_out_264=fopen(input->outall_264,"wb"))==NULL)
	{
		printf("Output_all_c %s does not exist\n",input->outall_264);
		exit(0);
	}

	if (1<input->num_regions)
	{
		if ((fp_plane_c=fopen(input->infile_c_plane,"rb"))==NULL)
		{
			printf("Input file_Plane_C %s does not exist\n",input->infile_c_plane);
			exit(0);
		}
		for(i=0;i<input->num_regions;i++)
		{
			if((fp_out_c_rec[i]=fopen(input->outfile_c_rec[i],"wb"))==NULL)
			{
				printf(" Outfile_c_rec %s does not exist\n",input->outfile_c_rec[i]);
			    exit(0);
			}
			if (strlen (input->outfile_c[i]) > 0 && (fp_out_c[i]=fopen(input->outfile_c[i],"wb"))==NULL)
			{
				printf("Error open file_C %s\n", input->outfile_c[i]);
				exit(0);
			}
		}
	}
	else
	{
		if ((fp_out_all_c=fopen(input->outfile_all_c,"wb"))==NULL)
		{
			printf("Error open file_C %s\n", input->outfile_all_c);
			exit(0);
		}
	}
    /////////////////////////////////////////////if encode L, open the L video
	if (1==input->right)
	{
		if ((fp_in_r=fopen(input->infile_r,"rb"))==NULL)
		{
			printf("Input file_R %s does not exist\n",input->infile_r);
			exit(0);
		}
		if ((fp_out_all_r_rec=fopen(input->outall_r,"wb"))==NULL)
		{
			printf("Output_all_r %s does not exist\n",input->outall_r);
			exit(0);
		}
		
		if (1<input->num_regions)
		{
			if ((fp_plane_r=fopen(input->infile_r_plane,"rb"))==NULL)
			{
				printf("Input file_Plane_R %s does not exist\n",input->infile_r_plane);
				exit(0);
			}
			for(i=0;i<input->num_regions;i++)
			{
				if((fp_out_r_rec[i]=fopen(input->outfile_r_rec[i],"wb"))==NULL)
				{
					printf(" Outfile_r_rec %s does not exist\n",input->outfile_r_rec[i]);
					exit(0);
				}
				if (strlen (input->outfile_r[i]) > 0 && (fp_out_r[i]=fopen(input->outfile_r[i],"wb"))==NULL)
				{
					printf("Error open file_R %s\n", input->outfile_r[i]);
					exit(0);
				}
			}
		}
		else
		{
			if ((fp_out_all_r=fopen(input->outfile_all_r,"wb"))==NULL)
			{
				printf("Error open file_R %s\n", input->outfile_all_r);
				exit(0);
			}
		}
	}

	if (1==input->left)
	{
		if ((fp_in_l=fopen(input->infile_l,"rb"))==NULL)
		{
			printf("Input file_L %s does not exist\n",input->infile_l);
			exit(0);
		}
		if ((fp_out_all_l_rec=fopen(input->outall_l,"wb"))==NULL)
		{
			printf("Output_all_l %s does not exist\n",input->outall_l);
			exit(0);
		}
		
		if (1<input->num_regions)
		{
			if ((fp_plane_l=fopen(input->infile_l_plane,"rb"))==NULL)
			{
				printf("Input file_Plane_L %s does not exist\n",input->infile_l_plane);
				exit(0);
			}
			for(i=0;i<input->num_regions;i++)
			{
				if((fp_out_l_rec[i]=fopen(input->outfile_l_rec[i],"wb"))==NULL)
				{
					printf(" Outfile_l_rec %s does not exist\n",input->outfile_l_rec[i]);
					exit(0);
				}
				if (strlen (input->outfile_l[i]) > 0 && (fp_out_l[i]=fopen(input->outfile_l[i],"wb"))==NULL)
				{
					printf("Error open file_L %s\n", input->outfile_l[i]);
					exit(0);
				}
			}
		}
		else
		{
			if ((fp_out_all_l=fopen(input->outfile_all_l,"wb"))==NULL)
			{
				printf("Error open file_L %s\n", input->outfile_all_l);
				exit(0);
			}
		}
	}

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
	// Set block sizes
	
	input->blc_size[0][0]=16;
	input->blc_size[0][1]=16;
	
	input->blc_size[1][0]=16;
	input->blc_size[1][1]=16;
	
	input->blc_size[2][0]=16;
	input->blc_size[2][1]= 8;
	
	input->blc_size[3][0]= 8;
	input->blc_size[3][1]=16;
	
	input->blc_size[4][0]= 8;
	input->blc_size[4][1]= 8;
	
	input->blc_size[5][0]= 8;
	input->blc_size[5][1]= 4;
	
	input->blc_size[6][0]= 4;
	input->blc_size[6][1]= 8;
	
	input->blc_size[7][0]= 4;
	input->blc_size[7][1]= 4;
	
	// set proper log2_max_frame_num_minus4.
	{
		int storedBplus1 = (input->StoredBPictures ) ? input->successive_Bframe + 1: 1;
		
		log2_max_frame_num_minus4 = max( (int)(CeilLog2(1+ input->no_frames_h264 *storedBplus1 ))-4, 0);
		
		//    log2_max_frame_num_minus4 = 0;
	}
	
	log2_max_pic_order_cnt_lsb_minus4 = max( (int)(CeilLog2(1+2*input->no_frames_h264 * (input->successive_Bframe + 1))) -4, 0);
	

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////


} 	
/////////////////////////////////////////////////////////////////////////////////////////////

#define MAX_ITEMS_TO_PARSE  10000


/*!
 ***********************************************************************
 * \brief
 *   print help message and exit
 ***********************************************************************
 */
void JMHelpExit()
{
  fprintf( stderr, "\n   lencod [-h] [-p defenc.cfg] {[-f curenc1.cfg]...[-f curencN.cfg]}"
    " {[-p EncParam1=EncValue1]..[-p EncParamM=EncValueM]}\n\n"    
    "## Parameters\n\n"

    "## Options\n"
    "   -h :  prints function usage\n"
    "   -d :  use <defenc.cfg> as default file for parameter initializations.\n"
    "         If not used then file defaults to encoder.cfg in local directory.\n"
    "   -f :  read <curencM.cfg> for reseting selected encoder parameters.\n"
    "         Multiple files could be used that set different parameters\n"
    "   -p :  Set parameter <EncParamM> to <EncValueM>.\n"
    "         See default encoder.cfg file for description of all parameters.\n\n"
    
    "## Supported video file formats\n"
    "   RAW:  .yuv -> YUV 4:2:0\n\n"
    
    "## Examples of usage:\n"
    "   lencod\n"
    "   lencod  -h\n"
    "   lencod  -d default.cfg\n"
    "   lencod  -f curenc1.cfg\n"
    "   lencod  -f curenc1.cfg -p InputFile=\"e:\\data\\container_qcif_30.yuv\" -p SourceWidth=176 -p SourceHeight=144\n"  
    "   lencod  -f curenc1.cfg -p FramesToBeEncoded=30 -p QPFirstFrame=28 -p QPRemainingFrame=28 -p QPBPicture=30\n");

  exit(-1);
}

/*!
 ***********************************************************************
 * \brief
 *    Parse the command line parameters and read the config files.
 * \param ac
 *    number of command line parameters
 * \param av
 *    command line parameters
 ***********************************************************************
 */
void Configure_h264(int ac, char *av[])
{
  char *content;
  int CLcount, ContentLen, NumberParams;
  char *filename=DEFAULTCONFIGFILENAME;

  memset (&configinput, 0, sizeof (InputParameters));
  //Set some initial parameters.
  configinput.LevelIDC   = LEVEL_IDC;
  configinput.ProfileIDC = PROFILE_IDC;
  // Process default config file
  CLcount = 1;

  if (ac==2)
  {
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMHelpExit();
    }
  }

  if (ac>=3)
  {
    if (0 == strncmp (av[1], "-d", 2))
    {
      filename=av[2];
      CLcount = 3;
    }
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMHelpExit();
    }
  }
//   printf ("Parsing Configfile %s", filename);ZZZZZZZZZZZZZZZZZZZZZZZZZ
  content = GetConfigFileContent (filename);
  ParseContent (content, strlen(content));
//   printf ("\n");ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
  free (content);

  // Parse the command line

  while (CLcount < ac)
  {
    if (0 == strncmp (av[CLcount], "-h", 2))
    {
      JMHelpExit();
    }
    
    if (0 == strncmp (av[CLcount], "-f", 2))  // A file parameter?
    {
      content = GetConfigFileContent (av[CLcount+1]);
      printf ("Parsing Configfile %s", av[CLcount+1]);
      ParseContent (content, strlen (content));
      printf ("\n");
      free (content);
      CLcount += 2;
    } else
    {
      if (0 == strncmp (av[CLcount], "-p", 2))  // A config change?
      {
        // Collect all data until next parameter (starting with -<x> (x is any character)),
        // put it into content, and parse content.

        CLcount++;
        ContentLen = 0;
        NumberParams = CLcount;

        // determine the necessary size for content
        while (NumberParams < ac && av[NumberParams][0] != '-')
          ContentLen += strlen (av[NumberParams++]);        // Space for all the strings
        ContentLen += 1000;                     // Additional 1000 bytes for spaces and \0s


        if ((content = malloc (ContentLen))==NULL) no_mem_exit("Configure: content");;
        content[0] = '\0';

        // concatenate all parameters identified before

        while (CLcount < NumberParams)
        {
          char *source = &av[CLcount][0];
          char *destin = &content[strlen (content)];

          while (*source != '\0')
          {
            if (*source == '=')  // The Parser expects whitespace before and after '='
            {
              *destin++=' '; *destin++='='; *destin++=' ';  // Hence make sure we add it
            } else
              *destin++=*source;
            source++;
          }
          *destin = '\0';
          CLcount++;
        }
//         printf ("Parsing command line string '%s'", content);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
        ParseContent (content, strlen(content));
        free (content);
//         printf ("\n");ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
      }
      else
      {
        printf (errortext, ET_SIZE, "Error in command line, ac %d, around string '%s', missing -f or -p parameters?", CLcount, av[CLcount]);
        error (errortext, 300);
      }
    }
  }
  printf ("\n");
  PatchInp();//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZz因为和分形的重复了。暂时删除掉，想试试test_rec有关系不？
}
void PatchInputNoFrames()
{
	// Tian Dong: May 31, 2002
	// If the frames are grouped into two layers, "FramesToBeEncoded" in the config file
	// will give the number of frames which are in the base layer. Here we let input->no_frames_h264
	// be the total frame numbers.
	input->no_frames_h264 = 1+ (input->no_frames_h264-1) * (input->NumFramesInELSubSeq+1);
	if ( input->NumFrameIn2ndIGOP )
		input->NumFrameIn2ndIGOP = 1+(input->NumFrameIn2ndIGOP-1) * (input->NumFramesInELSubSeq+1);
	FirstFrameIn2ndIGOP = input->no_frames_h264;
}
