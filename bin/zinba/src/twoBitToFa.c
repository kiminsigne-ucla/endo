/* twoBitToFa - Convert all or part of twoBit file to fasta. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dnaseq.h"
#include "fa.h"
#include "twoBit.h"
#include <R.h>
#include "bPlusTree.h"

static char const rcsid[] = "$Id: twoBitToFa.c,v 1.12 2009/01/09 10:14:33 kent Exp $";
/*
void usage()
/* Explain usage and exit. 
{
errAbort(
  "twoBitToFa - Convert all or part of .2bit file to fasta\n"
  "usage:\n"
  "   twoBitToFa input.2bit output.fa\n"
  "options:\n"
  "   -seq=name - restrict this to just one sequence\n"
  "   -start=X  - start at given position in sequence (zero-based)\n"
  "   -end=X - end at given position in sequence (non-inclusive)\n"
  "   -seqList=file - file containing list of the desired sequence names \n"
  "                    in the format seqSpec[:start-end], e.g. chr1 or chr1:0-189\n"
  "                    where coordinates are half-open zero-based, i.e. [start,end)\n"
  "   -noMask - convert sequence to all upper case\n"
  "   -bpt=index.bpt - use bpt index instead of built in one\n"
  "\n"
  "Sequence and range may also be specified as part of the input\n"
  "file name using the syntax:\n"
  "      /path/input.2bit:name\n"
  "   or\n"
  "      /path/input.2bit:name\n"
  "   or\n"
  "      /path/input.2bit:name:start-end\n"
  );
}
*/


//char *clSeq = NULL;	/* Command line sequence. */
//int clStart = 0;	/* Start from command line. */
//int clEnd = 0;		/* End from command line. */
char *clSeqList = NULL; /* file containing list of seq names */
bool noMask = FALSE;  /* convert seq to upper case */
//char *clBpt = NULL;	/* External index file. */


/*static struct optionSpec options[] = {
   {"seq", OPTION_STRING},
   {"seqList", OPTION_STRING},
   {"start", OPTION_INT},
   {"end", OPTION_INT},
   {"noMask", OPTION_BOOLEAN},
   {"bpt", OPTION_STRING},
   {NULL, 0},
};*/

void outputOne(struct twoBitFile *tbf, char *seqSpec, FILE *f, 
	int start, int end)
/* Output sequence. */
{
struct dnaSeq *seq = twoBitReadSeqFrag(tbf, seqSpec, start, end);
if (noMask)
    toUpperN(seq->dna, seq->size);
faWriteNext(f, seq->name, seq->dna, seq->size);
dnaSeqFree(&seq);
}

static void processAllSeqs(struct twoBitFile *tbf, FILE *outFile)
/* get all sequences in a file */
{
struct twoBitIndex *index;
for (index = tbf->indexList; index != NULL; index = index->next)
    outputOne(tbf, index->name, outFile, 0, 0);
}

static void processSeqSpecs(struct twoBitFile *tbf, struct twoBitSeqSpec *tbss,
                            FILE *outFile)
/* process list of twoBitSeqSpec objects */
{
struct twoBitSeqSpec *s;
for (s = tbss; s != NULL; s = s->next)
    outputOne(tbf, s->name, outFile, s->start, s->end);
}

void twoBitToFa(char **RinSeq,int *Rstart,int *Rstop,char **Rtwobitfile,char **RgcOut)
/* twoBitToFa - Convert all or part of twoBit file to fasta. */
{

	char *clSeq = RinSeq[0];	/* Command line sequence. */
	int clStart = *Rstart;	/* Start from command line. */
	int clEnd = *Rstop;		/* End from command line. */
	//char *clSeqList = Rtwobitfile[0]; /* file containing list of seq names */
	noMask = TRUE;  /* convert seq to upper case */
	char *clBpt = NULL;/* External index file. */
	char *inName = Rtwobitfile[0];	
	char *outName = RgcOut[0];
	
//char* inName=RinName[0];
//char* outName=RoutName[0] ;
struct twoBitFile *tbf;
FILE *outFile = mustOpen(outName, "w");
struct twoBitSpec *tbs;

if (clSeq != NULL)
    {
    char seqSpec[2*PATH_LEN];
    if (clEnd > clStart)
        safef(seqSpec, sizeof(seqSpec), "%s:%s:%d-%d", inName, clSeq, clStart, clEnd);
    else
        safef(seqSpec, sizeof(seqSpec), "%s:%s", inName, clSeq);
    tbs = twoBitSpecNew(seqSpec);
    }
else if (clSeqList != NULL)
    tbs = twoBitSpecNewFile(inName, clSeqList);
else
    tbs = twoBitSpecNew(inName);
if (tbs->seqs != NULL && clBpt != NULL)
    tbf = twoBitOpenExternalBptIndex(tbs->fileName, clBpt);
else
    tbf = twoBitOpen(tbs->fileName);
if (tbs->seqs == NULL)
    processAllSeqs(tbf, outFile);
else
    processSeqSpecs(tbf, tbs->seqs, outFile);
twoBitSpecFree(&tbs);
carefulClose(&outFile);
twoBitClose(&tbf);
}


void twoBitToFa2(char *clSeq,int clStart,int clEnd,char *inName,char *outName)
/* twoBitToFa - Convert all or part of twoBit file to fasta. */
{

	//char *clSeq = RinSeq[0];	/* Command line sequence. */
	//int clStart = Rstart;	/* Start from command line. */
	//int clEnd = Rstop;		/* End from command line. */
	//char *clSeqList = Rtwobitfile[0]; /* file containing list of seq names */
	noMask = TRUE;  /* convert seq to upper case */
	char *clBpt = NULL;/* External index file. */
	//char *inName = Rtwobitfile[0];	
	//char *outName = RgcOut[0];
	
//char* inName=RinName[0];
//char* outName=RoutName[0] ;
struct twoBitFile *tbf;
FILE *outFile = mustOpen(outName, "w");
struct twoBitSpec *tbs;

if (clSeq != NULL)
    {
    char seqSpec[2*PATH_LEN];
    if (clEnd > clStart)
        safef(seqSpec, sizeof(seqSpec), "%s:%s:%d-%d", inName, clSeq, clStart, clEnd);
    else
        safef(seqSpec, sizeof(seqSpec), "%s:%s", inName, clSeq);
    tbs = twoBitSpecNew(seqSpec);
    }
else if (clSeqList != NULL)
    tbs = twoBitSpecNewFile(inName, clSeqList);
else
    tbs = twoBitSpecNew(inName);
if (tbs->seqs != NULL && clBpt != NULL)
    tbf = twoBitOpenExternalBptIndex(tbs->fileName, clBpt);
else
    tbf = twoBitOpen(tbs->fileName);
if (tbs->seqs == NULL)
    processAllSeqs(tbf, outFile);
else
    processSeqSpecs(tbf, tbs->seqs, outFile);
twoBitSpecFree(&tbs);
carefulClose(&outFile);
twoBitClose(&tbf);
}

/*
int main(int argc, char *argv[])
/* Process command line. 
{
optionInit(&argc, argv, options);
if (argc != 3)
    usage();
clSeq = optionVal("seq", clSeq);
clStart = optionInt("start", clStart);
clEnd = optionInt("end", clEnd);
clSeqList = optionVal("seqList", clSeqList);
clBpt = optionVal("bpt", clBpt);
if ((clStart > clEnd) && (clSeq == NULL))
    errAbort("must specify -seq with -start and -end");
if ((clSeq != NULL) && (clSeqList != NULL))
    errAbort("can't specify both -seq and -seqList");
noMask = optionExists("noMask");
dnaUtilOpen();
twoBitToFa(argv[1], argv[2]);
return 0;
}*/
