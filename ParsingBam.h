#ifndef PARSINGBAM_H
#define PARSINGBAM_H

#include "Util.h"
#include "PhasingProcess.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#include <zlib.h>

enum VariantType {
    UNDEFINED = -1,
    SNP_HET = -2, // heterozygous Single Nucleotide Polymorphism
    SNP_HOM = -3, // homozygous Single Nucleotide Polymorphism
    INDEL_HET = -4, // heterozygous insertion/deletion
    INDEL_HOM = -5, // homozygous insertion/deletion
    SV_HET = -6, // heterozygous structure variation
    SV_HOM = -7, // homozygous structure variation
    MOD_HET_FORWARD_STRAND = -8, // heterozygous forward strand modification
    MOD_HET_REVERSE_STRAND = -9, // heterozygous reverse strand modification
    MOD_HET = -10, // heterozygous modification
    MOD_HOM = -11, // homozygous modification
};

struct RefAlt{
    std::string Ref;
    std::string Alt;
    bool is_reverse;
    bool is_modify;
};

class FastaParser{
    private:
        // file name
        std::string fastaFile ;
        std::vector<std::string> chrName;
        std::vector<int> last_pos;
    public:
        FastaParser(std::string fastaFile, std::vector<std::string> chrName, std::vector<int> last_pos, int numThreads);
        ~FastaParser();
        
        // chrName, chr string
        std::map<std::string, std::string > chrString;
    
};

enum Fields {
    CHROM = 0,
    POS = 1,
    ID = 2,
    REF = 3,
    ALT = 4,
    QUAL = 5,
    FILTER = 6,
    INFO = 7,
    FORMAT = 8,
    SAMPLE = 9
};

class FormatSample {

    private:
        std::vector<std::string> &fields;

        int findFlagColon(const int flagStart);
        int findValueStart(const int flagColon);

        void eraseColon(const Fields field, const int modifyStart);
        void resetGTValue(int modifyStart);

        void addValues(const std::string& flag, const std::string& value);
        void setGTValue(int modifyStart, const std::string& value);

    public:
        FormatSample(std::vector<std::string>& fields);
        void eraseFormatSample(const std::string& flag);
        void addFlagAndValue(const std::string& flagBase, const int value, const std::string& flagAdd = "");
        void setGTFlagAndValue(const std::string& flagBase, const std::string& value, const std::string& flagAdd = "");
};

struct FormatInfo {
    std::string id;
    std::string number;
    std::string type;
    std::string description;
    bool is_present;
};

using FormatDefs = std::vector<FormatInfo>;

class BaseVairantParser{

    protected:
        virtual bool checkType(const std::string& chr, int pos) const = 0;
        FormatDefs formatDefs = {
            {"GT", "1", "String", "Genotype", false},
            {"PS", "1", "Integer", "Phase set identifier", false},
            // {"SH", "1", "Integer", "Source haplotype", false},
            // {"GT2", "1", "String", "Sub genotype", false},
            // {"PS2", "1", "Integer", "Sub phase set identifier", false},
            // {"SH2", "1", "Integer", "Sub source haplotype", false},
        };
        bool commandLine;
        PhasingParameters *params;

    public:
        BaseVairantParser();
        virtual ~BaseVairantParser();
        // input parser
        void compressParser(std::string &variantFile);
        void unCompressParser(std::string &variantFile);
        virtual void parserProcess(std::string &input)=0;
        // output parser
        void compressInput(std::string variantFile, std::string resultFile, ChrPhasingResult &chrPhasingResult);
        void unCompressInput(std::string variantFile, std::string resultFile, ChrPhasingResult &chrPhasingResult);
        virtual void writeLine(std::string &input, std::ofstream &resultVcf, ChrPhasingResult &chrPhasingResult);
};

class SnpParser : public BaseVairantParser{
    
    private:
        // chr, variant position (0-base), allele haplotype
        std::map<std::string, std::map<int, RefAlt> > *chrVariant;
        // id and idx
        std::vector<std::string> chrName;
        // chr, variant position (0-base)
        std::map<std::string, std::map<int, bool> > chrVariantHomopolymer;
        
        // override input parser
        void parserProcess(std::string &input);

    public:

        SnpParser(PhasingParameters &in_params);
        ~SnpParser();
            
        std::map<int, RefAlt> getVariants(std::string chrName);  

        std::vector<std::string> getChrVec();
        
        bool findChromosome(std::string chrName);
        
        int getLastSNP(std::string chrName);
        
        void writeResult(ChrPhasingResult &chrPhasingResult);

        bool findSNP(std::string chr, int posistion);
        
        void filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec, std::string &chr_reference);

        bool checkType(const std::string& chr, int pos) const override;
};

class SVParser : public BaseVairantParser{
    
    private:
        SnpParser *snpFile;

        // chr , variant position (0-base), read
        std::map<std::string, std::map<int, std::map<std::string ,bool> > > *chrVariant;
        // chr, variant position (0-base)
        std::map<std::string, std::map<int, bool> > posDuplicate;
        
        // override input parser
        void parserProcess(std::string &input);
        
    public:
    
        SVParser(PhasingParameters &params, SnpParser &snpFile);
        ~SVParser();
            
        std::map<int, std::map<std::string ,bool> > getVariants(std::string chrName);  

        void writeResult(ChrPhasingResult &chrPhasingResult);

        bool findSV(std::string chr, int posistion);

        bool checkType(const std::string& chr, int pos) const override;
};

class METHParser : public BaseVairantParser{
    
    private:
        SnpParser *snpFile;
        SVParser *svFile;
        
        int representativePos;
        int upMethPos;
        
        // chr , variant position (0-base), read, (is reverse strand)
        std::map<std::string, std::map<int, std::map<std::string ,RefAlt> > > *chrVariant;
        // In a series of consecutive methylation positions, 
        // the first methylation position will be used as the representative after merging.
        // This map is used to find the coordinates that originally represented itself
        std::map<int, int > *representativeMap;

        // override input parser
        void parserProcess(std::string &input);
        
    public:
        
        std::map<int, std::map<std::string ,RefAlt> > getVariants(std::string chrName);  
        
        METHParser(PhasingParameters &params, SnpParser &snpFile, SVParser &svFile);
        ~METHParser();
		
		void writeResult(ChrPhasingResult &chrPhasingResult);

        bool checkType(const std::string& chr, int pos) const override;
};

struct Alignment{
    std::string chr;
    std::string qname;
    int refStart;
    int qlen;
    char *qseq;
    int cigar_len;
    int *op;
    int *ol;
    char *quality;
    bool is_reverse;
};

class BamParser{
    
    private:
        std::string chrName;
        std::vector<std::string> BamFileVec;
        // SNP map and iter
        std::map<int, RefAlt> *currentVariants;
        std::map<int, RefAlt>::iterator firstVariantIter;
        // SV map and iter
        std::map<int, std::map<std::string ,bool> > *currentSV;
        std::map<int, std::map<std::string ,bool> >::iterator firstSVIter;
        // mod map and iter
        std::map<int, std::map<std::string ,RefAlt> > *currentMod;
        std::map<int, std::map<std::string ,RefAlt> >::iterator firstModIter;
        void get_snp(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::vector<ReadVariant> &readVariantVec, const std::string &ref_string, bool isONT);
   
    public:
        BamParser(std::string chrName, std::vector<std::string> inputBamFileVec, SnpParser &snpMap, SVParser &svFile, METHParser &modFile);
        ~BamParser();
        
        void direct_detect_alleles(int lastSNPPos, htsThreadPool &threadPool, PhasingParameters params, std::vector<ReadVariant> &readVariantVec , const std::string &ref_string);

};


#endif