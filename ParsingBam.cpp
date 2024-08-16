#include "ParsingBam.h"

#include <string.h>
#include <sstream>


// vcf parser modify from 
// http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
// bam parser modify from 
// https://github.com/whatshap/whatshap/blob/9882248c722b1020321fea6e9491c1cc5b75354b/whatshap/variants.py


// FASTA
FastaParser::FastaParser(std::string fastaFile,  std::vector<std::string> chrName, std::vector<int> last_pos, int numThreads):
    fastaFile(fastaFile),
    chrName(chrName),
    last_pos(last_pos)
{
    // init map
    for(std::vector<std::string>::iterator iter = chrName.begin() ; iter != chrName.end() ; iter++)
        chrString.insert(std::make_pair( (*iter) , ""));

    // iterating all chr
    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for(std::vector<std::string>::iterator iter = chrName.begin() ; iter != chrName.end() ; iter++){
        
        faidx_t *fai = NULL;
        fai = fai_load(fastaFile.c_str());
        int index = iter - chrName.begin();

        // ref_len is a return value that is length of retrun string
        int ref_len = 0;
        // read file
        std::string chr_info(faidx_fetch_seq(fai , (*iter).c_str() , 0 ,last_pos.at(index)+5 , &ref_len));
        if(ref_len == 0){
            std::cout<<"nothing in reference file \n";
        }
        // update map
        chrString[(*iter)] = chr_info;

    }
}

FastaParser::~FastaParser(){
 
}

FormatSample::FormatSample(std::vector<std::string> &fields) : fields(fields) {}

int FormatSample::findFlagColon(const int flagStart) {
    int flagColon = 0;
    for (int i = 0; i < flagStart; ++i) {
        if (fields[FORMAT][i] == ':') flagColon++;
    }
    return flagColon;
}

int FormatSample::findValueStart(const int flagStart) {
    int flagColon = findFlagColon(flagStart);
    int currentColon = 0;
    int valueStart = 0;
    for (unsigned int i = 0; i < fields[SAMPLE].length(); ++i) {
        if (currentColon >= flagColon) break;
        if (fields[SAMPLE][i] == ':') currentColon++;
        valueStart++;
    }
    return valueStart;
}

void FormatSample::eraseColon(const Fields field, const int modifyStart) {
    auto endPos = fields[field].find(":", modifyStart + 1);
    if (endPos != std::string::npos) {
        fields[field].erase(modifyStart, endPos - modifyStart + 1);
    } else {
        fields[field].erase(modifyStart - 1, fields[field].length() - modifyStart + 1);
    }
}

void FormatSample::resetGTValue(int modifyStart) {
    if (fields[SAMPLE][modifyStart + 1] == '|') {
        char &firstChar = fields[SAMPLE][modifyStart];
        char &secondChar = fields[SAMPLE][modifyStart + 2];

        if ((firstChar == '0' && secondChar == '.') || (firstChar == '.' && secondChar == '0')) {
            firstChar = '0';
            secondChar = '1';
        } else if ((firstChar == '1' && secondChar == '.') || (firstChar == '.' && secondChar == '1')) {
            firstChar = '1';
            secondChar = '1';
        } else if (firstChar > secondChar) {
            std::swap(firstChar, secondChar);
        }
        fields[SAMPLE][modifyStart + 1] = '/';
    }
}

void FormatSample::addValues(const std::string &flag, const std::string &value) {
    fields[FORMAT] += ":" + flag;
    fields[SAMPLE] += ":" + value;
}

void FormatSample::setGTValue(int modifyStart, const std::string &value) {
    fields[SAMPLE][modifyStart] = value[0];
    fields[SAMPLE][modifyStart+1] = value[1];
    fields[SAMPLE][modifyStart+2] = value[2];
}

void FormatSample::eraseFormatSample(const std::string &flag) {
    auto flagStart = fields[FORMAT].find(flag);
    if (flagStart != std::string::npos) {
        int valueStart = findValueStart(flagStart);
        if (flag != "GT") {
            eraseColon(FORMAT, flagStart);
            eraseColon(SAMPLE, valueStart);
        } else {
            resetGTValue(valueStart);
        }
    }
}

void FormatSample::addFlagAndValue(const std::string &flagBase, const int inValue, const std::string &flagAdd) {
    std::string value = inValue != 0 ? std::to_string(inValue) : ".";
    addValues(flagBase + flagAdd, value);
}

void FormatSample::setGTFlagAndValue(const std::string &flagBase, const std::string &value, const std::string &flagAdd) {
    const std::string mergeFlag = flagBase + flagAdd;
    auto flagStart = fields[FORMAT].find(mergeFlag);
    if (flagStart != std::string::npos) {
        int valueStart = findValueStart(flagStart);
        if (value != ".") {
            setGTValue(valueStart, value);
        } else {
            resetGTValue(valueStart);
        }
    } else {
        if (value != ".") {
            addValues(mergeFlag, value);
        } else {
            addValues(mergeFlag, "./.");
        }
    }
}

BaseVairantParser::BaseVairantParser() : commandLine(false) {}

BaseVairantParser::~BaseVairantParser(){}

void BaseVairantParser::compressParser(std::string &variantFile){
    gzFile file = gzopen(variantFile.c_str(), "rb");
    if(variantFile=="")
        return;
    if(!file){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{  
        int buffer_size = 1048576; // 1M
        char* buffer = (char*) malloc(buffer_size);
        if(!buffer){
            std::cerr<<"Failed to allocate buffer\n";
            exit(EXIT_FAILURE);
        }
        char* offset = buffer;
            
        while(true) {
            int len = buffer_size - (offset - buffer);
            if (len == 0){
                buffer_size *= 2; // Double the buffer size
                char* new_buffer = (char*) realloc(buffer, buffer_size);
                if(!new_buffer){
                    std::cerr<<"Failed to allocate buffer\n";
                    free(buffer);
                    exit(EXIT_FAILURE);
                }
                buffer = new_buffer;
                offset = buffer + buffer_size / 2; // Update the offset pointer to the end of the old buffer
                len = buffer_size - (offset - buffer);
            }

            len = gzread(file, offset, len);
            if (len == 0) break;    
            if (len <  0){ 
                int err;
                fprintf (stderr, "Error: %s.\n", gzerror(file, &err));
                exit(EXIT_FAILURE);
            }

            char* cur = buffer;
            char* end = offset+len;
            for (char* eol; (cur<end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
            {
                std::string input = std::string(cur, eol);
                parserProcess(input);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
        free(buffer);
    }    
}

void BaseVairantParser::unCompressParser(std::string &variantFile){
    std::ifstream originVcf(variantFile);
    if(variantFile=="")
        return;
    if(!originVcf.is_open()){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
        exit(1);
    }
    else{
        std::string input;
        while(! originVcf.eof() ){
            std::getline(originVcf, input);
            parserProcess(input);
        }
    }
}

void BaseVairantParser::compressInput(std::string variantFile, std::string resultFile, ChrPhasingResult &chrPhasingResult){
    
    gzFile file = gzopen(variantFile.c_str(), "rb");
    std::ofstream resultVcf(resultFile);
    
    if(!resultVcf.is_open()){
        std::cout<< "Fail to open write file: " << resultFile << "\n";
    }
    else if(!file){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{
        int buffer_size = 1048576; // 1M
        char* buffer = (char*) malloc(buffer_size);
        if(!buffer){
            std::cerr<<"Failed to allocate buffer\n";
            exit(EXIT_FAILURE);
        }
        char* offset = buffer;
            
        while(true) {
            int len = buffer_size - (offset - buffer);
            if (len == 0){
                buffer_size *= 2; // Double the buffer size
                char* new_buffer = (char*) realloc(buffer, buffer_size);
                if(!new_buffer){
                    std::cerr<<"Failed to allocate buffer\n";
                    free(buffer);
                    exit(EXIT_FAILURE);
                }
                buffer = new_buffer;
                offset = buffer + buffer_size / 2; // Update the offset pointer to the end of the old buffer
                len = buffer_size - (offset - buffer);
            }

            len = gzread(file, offset, len);
            if (len == 0) break;    
            if (len <  0){ 
                int err;
                fprintf (stderr, "Error: %s.\n", gzerror(file, &err));
                exit(EXIT_FAILURE);
            }

            char* cur = buffer;
            char* end = offset+len;
            for (char* eol; (cur<end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
            {
                std::string input = std::string(cur, eol);
                //parserProcess(input);
                writeLine(input, resultVcf, chrPhasingResult);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
        free(buffer);
    }
}

void BaseVairantParser::unCompressInput(std::string variantFile, std::string resultFile, ChrPhasingResult &chrPhasingResult){
    std::ifstream originVcf(variantFile);
    std::ofstream resultVcf(resultFile);
    
    if(!resultVcf.is_open()){
        std::cout<< "Fail to open write file: " << resultFile << "\n";
    }
    else if(!originVcf.is_open()){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{
        std::string input;
        while(! originVcf.eof() ){
            std::getline(originVcf, input);
            if( input != "" ){
                writeLine(input, resultVcf, chrPhasingResult);
            }
        }
    }
}

void BaseVairantParser::writeLine(std::string &input, std::ofstream &resultVcf, ChrPhasingResult &chrPhasingResult){
    // header
    if( input.substr(0, 2) == "##" ){
        // avoid double definition
        for (auto formatDef = formatDefs.begin(); formatDef != formatDefs.end(); ++formatDef){
            if (formatDef->is_present) {
                continue;
            }
            std::string format_string = "##FORMAT=<ID=" + formatDef->id + ",";
            if( input.substr(0, 16) == format_string ){
                formatDef->is_present = true;
            }
        }
        resultVcf << input << "\n";
    }
    else if ( input.substr(0, 6) == "#CHROM" || input.substr(0, 6) == "#chrom" ){
        // format line 
        if( commandLine == false ){
            for (auto formatDef = formatDefs.begin(); formatDef != formatDefs.end(); ++formatDef){
                if (formatDef->is_present == false){
                    resultVcf<< "##FORMAT=<ID=" + formatDef->id + ",Number=" + formatDef->number + ",Type=" + formatDef->type + ",Description=\"" + formatDef->description + "\">\n";
                }
            }
            resultVcf << "##longphaseVersion=" << params->version << "\n";
            resultVcf << "##commandline=\"" << params->command << "\"\n";
            commandLine = true;
        }
        resultVcf << input << "\n";
    }
    else{
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
        if( fields.size() == 0 )
            return;
        FormatSample formatSample(fields);

        // if flag already exist, reset flag info
        for (auto formatDef = formatDefs.begin(); formatDef != formatDefs.end(); ++formatDef){
            if (formatDef->is_present == true){
                formatSample.eraseFormatSample(formatDef->id);
            }
        }

        int pos = std::stoi(fields[POS]) - 1;
        const auto& posPhasingResult = chrPhasingResult[fields[CHROM]];
        const auto& posPhasingResultIter = posPhasingResult.find(pos);
        // Check if the variant is extracted from this VCF
        bool isType = checkType(fields[CHROM], pos);
        std::string flagAdd = "";

        // this pos is phase
        if( posPhasingResultIter != posPhasingResult.end() && isType ){
            const auto& phasingResult = posPhasingResultIter->second;
            auto refHaplotype = phasingResult.refHaplotype;
            auto altHaplotype = refHaplotype == 0 ? 1 :0;
            std::string gtValue = std::to_string(refHaplotype) + "|" + std::to_string(altHaplotype);
            formatSample.setGTFlagAndValue("GT", gtValue, flagAdd);
            formatSample.addFlagAndValue("PS", phasingResult.phaseSet + 1, flagAdd);
        }
        // this pos has not been phased
        else{
            formatSample.addFlagAndValue("PS", 0, flagAdd);
        }

        for(std::vector<std::string>::iterator fieldIter = fields.begin(); fieldIter != fields.end(); ++fieldIter){
            if( fieldIter != fields.begin() )
                resultVcf<< "\t";
            resultVcf << (*fieldIter);
        }
        resultVcf << "\n";
    }
}

bool SnpParser::checkType(const std::string& chr, int pos) const {
    return (*chrVariant)[chr].find(pos)!= (*chrVariant)[chr].end();
}
bool SVParser::checkType(const std::string& chr, int pos) const {
    return (*chrVariant)[chr].find(pos)!= (*chrVariant)[chr].end();
}
bool METHParser::checkType(const std::string& chr, int pos) const {
    return (*chrVariant)[chr].find(pos)!= (*chrVariant)[chr].end();
}

// SNP
SnpParser::SnpParser(PhasingParameters &in_params){

    chrVariant = new std::map<std::string, std::map<int, RefAlt> >;
    
    params = &in_params;
    // open vcf file
    htsFile * inf = bcf_open(params->snpFile.c_str(), "r");
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    // counters
    int nseq = 0;
    // report names of all the sequences in the VCF file
    const char **seqnames = NULL;
    // chromosome idx and name
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    // store chromosome
    for (int i = 0; i < nseq; i++) {
        // bcf_hdr_id2name is another way to get the name of a sequence
        chrName.push_back(seqnames[i]);
    }
    // set all sample string
    std::string allSmples = "-";
    // limit the VCF data to the sample name passed in
    int is_file = bcf_hdr_set_samples(hdr, allSmples.c_str(), 0);
    if( is_file != 0 ){
        std::cout << "error or a positive integer if the list contains samples not present in the VCF header\n";
    }

    // struct for storing each record
    bcf1_t *rec = bcf_init();
    VariantType parserType = UNDEFINED;

    while (bcf_read(inf, hdr, rec) == 0) {
        // snp or indel
        if (bcf_is_snp(rec) || params->phaseIndel) {
            parserType = confirmRequiredGT(hdr, rec, "GT", rec->pos);
            if ( parserType != UNDEFINED ) {
                // get chromosome string
                std::string chr = seqnames[rec->rid];
                recordVariant(chr, rec, parserType, chrVariant);
            }
        }
    }
}

SnpParser::~SnpParser(){
    delete chrVariant;
}

std::map<int, RefAlt> SnpParser::getVariants(std::string chrName){
    std::map<int, RefAlt> targetVariants;
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant->find(chrName);
    
    if( chrIter != chrVariant->end() )
        targetVariants = (*chrIter).second;
    
    return targetVariants;
}

std::vector<std::string> SnpParser::getChrVec(){
    return chrName;
}

/*bool SnpParser::findChromosome(std::string chrName){
    std::map<std::string, std::map< int, RefAlt > >::iterator chrVariantIter = chrVariant->find(chrName);
    // this chromosome not exist in this file. 
    if(chrVariantIter == chrVariant->end())
        return false;
    return true;
}*/

int SnpParser::getLastSNP(std::string chrName){
    std::map<std::string, std::map< int, RefAlt > >::iterator chrVariantIter = chrVariant->find(chrName);
    // this chromosome not exist in this file. 
    if(chrVariantIter == chrVariant->end())
        return -1;
    // get last SNP
    std::map< int, RefAlt >::reverse_iterator lastVariantIter = (*chrVariantIter).second.rbegin();
    // there are no SNPs on this chromosome
    if(lastVariantIter == (*chrVariantIter).second.rend())
        return -1;
    return (*lastVariantIter).first;
}

void SnpParser::writeResult(ChrPhasingResult &chrPhasingResult){

    if( params->snpFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressInput(params->snpFile, params->resultPrefix+".vcf", chrPhasingResult);
    }
    else if( params->snpFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressInput(params->snpFile, params->resultPrefix+".vcf", chrPhasingResult);
    }
    return;
}

void SnpParser::parserProcess(std::string &input){
}

void SnpParser::fetchAndValidateTag(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int **dst, int *ndst, hts_pos_t pos){
    int checkTag = bcf_get_format_int32(hdr, line, tag, dst, ndst);

    if(checkTag<0){
        std::cerr<< "pos " << pos << " missing " << tag <<" value" << "\n";
        exit(1);
    }
}

VariantType SnpParser::confirmRequiredGT(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, hts_pos_t pos){
    
    int ngt_arr = 0;
    int *gt = NULL;
    fetchAndValidateTag(hdr, line, tag, &gt, &ngt_arr, pos);
    
    // hetero SNP
    if ((gt[0] == 2 && gt[1] == 4) || (gt[0] == 4 && gt[1] == 2) || // 0/1, 1/0
        (gt[0] == 2 && gt[1] == 5) || (gt[0] == 4 && gt[1] == 3) || // 0|1, 1|0
        (gt[0] == 0 && gt[1] == 3) || (gt[0] == 2 && gt[1] == 1)) { // .|0, 0|.
        return SNP_HET;
    }
    // homo SNP
    else if ((gt[0] == 4 && gt[1] == 4) || // 1/1
             (gt[0] == 4 && gt[1] == 5) || // 1|1
             (gt[0] == 0 && gt[1] == 5) || (gt[0] == 4 && gt[1] == 1)) { // .|1, 1|.
        return SNP_HOM;
    }
    else {
        return UNDEFINED;
    }
}

void SnpParser::recordVariant(std::string &chr, bcf1_t *rec, VariantType parserType, std::map<std::string, std::map<int, RefAlt> > *chrVariant) {
    // position is 0-base
    int variantPos = rec->pos;
    // get r alleles
    RefAlt tmp;
    tmp.Ref = rec->d.allele[0];
    tmp.Alt = rec->d.allele[1];
    tmp.is_homozygous = false;
    if (parserType == SNP_HOM) {
        tmp.is_homozygous = true;
    }
    
    //prevent the MAVs calling error which makes the GT=0/1
    if ( rec->d.allele[1][tmp.Alt.size()+1] != '\0' ){
        return;
    }
    
    // record 
    (*chrVariant)[chr][variantPos] = tmp;
}

bool SnpParser::findSNP(std::string chr, int position){
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant->find(chr);
    // empty chromosome
    if( chrIter == chrVariant->end() )
        return false;
    
    std::map<int, RefAlt>::iterator posIter = (*chrVariant)[chr].find(position);
    // empty position
    if( posIter == (*chrVariant)[chr].end() )
        return false;

    // If the variant is homozygous (is_homozygous is true), remove the SNP at this position from chrVariant 
    // and prioritize other variants.
    if (posIter->second.is_homozygous){
        (*chrVariant)[chr].erase(posIter);
        return false;
    }
        
    return true;
}

void SnpParser::filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec, std::string &chr_reference){
    
    // pos, <allele, <strand, True>
    std::map<int, std::map<int, std::map<int, bool> > > posAlleleStrand;
    std::map< int, bool > methylation;
    
    /*
    // iter all variant, record the strand contained in each SNP
    for( auto readSNPVecIter : readVariantVec ){
        // tag allele on forward or reverse strand
        for(auto variantIter : readSNPVecIter.variantVec ){
            posAlleleStrand[variantIter.position][variantIter.allele][readSNPVecIter.is_reverse] = true;
        }
    }
    
    // iter all SNP, both alleles that require SNP need to appear in the two strand
    for(auto pos: posAlleleStrand){
        // this position contain two allele, REF allele appear in the two strand, ALT allele appear in the two strand
        if(pos.second.size() == 2 && pos.second[0].size() == 2 && pos.second[1].size() == 2 ){
            // high confident SNP
        }
        else{
            //methylation[pos.first] = true;
        }
    }
    */

    // Filter SNPs that are not easy to phasing due to homopolymer
    // get variant list
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter =  chrVariant->find(chr);
    std::map< int, bool > errorProneSNP;
    
    if( chrIter != chrVariant->end() ){
        std::map<int, int> consecutiveAllele;
        // iter all SNP and tag homopolymer
        for(std::map<int, RefAlt>::iterator posIter = (*chrIter).second.begin(); posIter != (*chrIter).second.end(); posIter++ ){
            consecutiveAllele[(*posIter).first] = homopolymerLength((*posIter).first, chr_reference);
        }
        std::map<int, RefAlt>::iterator currSNPIter = (*chrIter).second.begin();
        std::map<int, RefAlt>::iterator nextSNPIter = std::next(currSNPIter,1);
        // check whether each SNP pair falls in an area that is not easy to phasing
        while( currSNPIter != (*chrIter).second.end() && nextSNPIter != (*chrIter).second.end() ){
            if (!currSNPIter->second.is_homozygous && !nextSNPIter->second.is_homozygous){
                int currPos = (*currSNPIter).first;
                int nextPos = (*nextSNPIter).first;
                // filter one of SNP if this SNP pair falls in homopolymer and distance<=2
                if( consecutiveAllele[currPos] >= 3 && consecutiveAllele[nextPos] >= 3 && std::abs(currPos-nextPos)<=2 ){
                    errorProneSNP[nextPos]=true;
                nextSNPIter = (*chrIter).second.erase(nextSNPIter);
                    continue;
                }
            }
            currSNPIter++;
            nextSNPIter++;
        }
        
    }
    
    // iter all reads
    for( std::vector<ReadVariant>::iterator readSNPVecIter = readVariantVec.begin() ; readSNPVecIter != readVariantVec.end() ; readSNPVecIter++ ){
        // iter all SNPs in this read
        for(std::vector<Variant>::iterator variantIter = (*readSNPVecIter).variantVec.begin() ; variantIter != (*readSNPVecIter).variantVec.end() ; ){
            std::map< int, bool >::iterator delSNPIter = methylation.find((*variantIter).position);
            std::map< int, bool >::iterator homoIter = errorProneSNP.find((*variantIter).position);
            
            if( delSNPIter != methylation.end() ){
                variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
            }
            else if( homoIter != errorProneSNP.end() ){
                variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
            }
            else{
                variantIter++;
            }
        }
    }
}

// SV
SVParser::SVParser(PhasingParameters &in_params, SnpParser &in_snpFile){
    params = &in_params;
    snpFile = &in_snpFile;
    
    chrVariant = new std::map<std::string, std::map<int, std::map<std::string ,bool> > >;
    
    if( params->svFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(params->svFile);
    }
    else if( params->svFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(params->svFile);
    }

    // erase SV pos if this pos appear two or more times
    for(std::map<std::string, std::map<int, bool> >::iterator chrIter = posDuplicate.begin(); chrIter != posDuplicate.end() ; chrIter++){
        for(std::map<int, bool>::iterator posIter = (*chrIter).second.begin() ; posIter != (*chrIter).second.end() ; posIter++ ){
            if( (*posIter).second == true){
                std::map<int, std::map<std::string ,bool> >::iterator erasePosIter = (*chrVariant)[(*chrIter).first].find((*posIter).first);
                if(erasePosIter != (*chrVariant)[(*chrIter).first].end()){
                    (*chrVariant)[(*chrIter).first].erase(erasePosIter);
                }
            }
        }
    }
}

SVParser::~SVParser(){ 
    delete chrVariant;
}

void SVParser::parserProcess(std::string &input){
    if( input.substr(0, 2) == "##" ){

    }
    else if ( input.substr(0, 1) == "#" ){
        
    }
    else{
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;
        // trans to 0-base
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];
                
        // find GT flag
        int colon_pos = 0;
        int gt_pos = fields[8].find("GT");
        for(int i =0 ; i< gt_pos ; i++){
            if(fields[8][i]==':')
                colon_pos++;
        }
        // find GT value start
        int current_colon = 0;
        int modify_start = 0;
        for(unsigned int i =0; i < fields[9].length() ; i++){
            if( current_colon >= colon_pos )
                break;
            if(fields[9][i]==':')
                current_colon++;  
            modify_start++;
        }
        bool filter = false;
                
        // homo GT
        if(fields[9][modify_start]==fields[9][modify_start+2]){
            filter = true;
        }
        // conflict pos with SNP
        if( (*snpFile).findSNP(chr,pos) ){
            filter = true;
        }
                
        std::map<int, bool>::iterator posIter = posDuplicate[chr].find(pos);
        // conflict pos with SV
        if( posIter == posDuplicate[chr].end() )
            posDuplicate[chr][pos] = false;
        else{
            posDuplicate[chr][pos] = true;
            filter = true;
        }
                
        if(filter){
            return;
        }
        
        // get read INFO
        int read_pos = fields[7].find("RNAMES=");
        // detected RNAMES included in INFO
        if( read_pos != -1 ){
            // Extract the position of "=" in "RNAMES="
            read_pos = fields[7].find("=",read_pos);
            read_pos++;
            // Capture the position of the ";" symbol at the end of the information in "RNAMES="
            int next_field = fields[7].find(";",read_pos);
            // Capture the range of read IDs included in the entire RNAMES
            std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
            std::stringstream totalReadStream(totalRead);
            // Extract each read ID individually
            std::string read;
            while(std::getline(totalReadStream, read, ','))
            {
               (*chrVariant)[chr][pos][read]= true;
            }
        }
    }
}

std::map<int, std::map<std::string ,bool> > SVParser::getVariants(std::string chrName){
    std::map<int, std::map<std::string ,bool> > targetVariants;
    std::map<std::string, std::map<int, std::map<std::string ,bool> > >::iterator chrIter = chrVariant->find(chrName);
    
    if( chrIter != chrVariant->end() )
        targetVariants = (*chrIter).second;
        
    return targetVariants;
}

void SVParser::writeResult(ChrPhasingResult &chrPhasingResult){
    
    if( params->svFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressInput(params->svFile, params->resultPrefix+"_SV.vcf", chrPhasingResult);
    }
    else if( params->svFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressInput(params->svFile, params->resultPrefix+"_SV.vcf", chrPhasingResult);
    }
    return;
}

bool SVParser::findSV(std::string chr, int position){
    std::map<std::string, std::map<int, std::map<std::string ,bool> > >::iterator chrIter = chrVariant->find(chr);
    // empty chromosome
    if( chrIter == chrVariant->end() )
        return false;
    
    std::map<int, std::map<std::string ,bool>>::iterator posIter = (*chrVariant)[chr].find(position);
    // empty position
    if( posIter == (*chrVariant)[chr].end() )
        return false;
        
    return true;
}

BamParser::BamParser(std::string inputChrName, std::vector<std::string> inputBamFileVec, SnpParser &snpMap, SVParser &svFile, METHParser &modFile):chrName(inputChrName),BamFileVec(inputBamFileVec){
    
    currentVariants = new std::map<int, RefAlt>;
    currentSV = new std::map<int, std::map<std::string ,bool> >;
    currentMod = new std::map<int, std::map<std::string ,RefAlt> >;
    
    // use chromosome to find recorded snp map
    (*currentVariants) = snpMap.getVariants(chrName);
    // set skip variant start iterator
    firstVariantIter = currentVariants->begin();
    if( firstVariantIter == currentVariants->end() ){
        std::cerr<< "error chromosome name or empty map.\n";
        exit(1);
    }
    // set current chromosome SV map
    (*currentSV) = svFile.getVariants(chrName);
    firstSVIter = currentSV->begin();
    // set current chromosome MOD map
    (*currentMod) = modFile.getVariants(chrName);
    firstModIter = currentMod->begin();
    
}

BamParser::~BamParser(){
    delete currentVariants;
    delete currentSV;
    delete currentMod;
}

void BamParser::direct_detect_alleles(int lastSNPPos, htsThreadPool &threadPool, PhasingParameters params, std::vector<ReadVariant> &readVariantVec, ClipCount &clipCount, const std::string &ref_string){
    
    // record SNP start iter
    std::map<int, RefAlt>::iterator tmpFirstVariantIter = firstVariantIter;
    // record SV start iter
    std::map<int, std::map<std::string ,bool> >::iterator tmpFirstSVIter = firstSVIter;
    // record MOD start iter
    std::map<int, std::map<std::string ,RefAlt> >::iterator tmpFirstModIter = firstModIter;
    
    for( auto bamFile: BamFileVec ){
        
        firstVariantIter = tmpFirstVariantIter;
        firstSVIter = tmpFirstSVIter;
        firstModIter = tmpFirstModIter;

        // open bam file
        samFile *fp_in = hts_open(bamFile.c_str(),"r"); 
        // load reference file
        hts_set_fai_filename(fp_in, params.fastaFile.c_str() );
        // read header
        bam_hdr_t *bamHdr = sam_hdr_read(fp_in); 
        // initialize an alignment
        bam1_t *aln = bam_init1(); 
        hts_idx_t *idx = NULL;
        
        if ((idx = sam_index_load(fp_in, bamFile.c_str())) == 0) {
            std::cout<<"ERROR: Cannot open index for bam file\n";
            exit(1);
        }
        
        std::string range = chrName + ":1-" + std::to_string(lastSNPPos);
        hts_itr_t* iter = sam_itr_querys(idx, bamHdr, range.c_str());

        
        hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &threadPool);
        int result;
        while ((result = sam_itr_multi_next(fp_in, iter, aln)) >= 0) { 
            int flag = aln->core.flag;

            if (    aln->core.qual < params.mappingQuality  // mapping quality
                 || (flag & 0x4)   != 0  // read unmapped
                 || (flag & 0x100) != 0  // secondary alignment. repeat. 
                                         // A secondary alignment occurs when a given read could align reasonably well to more than one place.
                 || (flag & 0x400) != 0  // duplicate 
                 // (flag & 0x800) != 0  // supplementary alignment
                                         // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
                 ){
                continue;
            }

            get_snp(*bamHdr, *aln, readVariantVec, clipCount, ref_string, params.isONT);
        }
        hts_idx_destroy(idx);
        bam_hdr_destroy(bamHdr);
        bam_destroy1(aln);
        sam_close(fp_in);
    }
    
}

void BamParser::get_snp(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::vector<ReadVariant> &readVariantVec, ClipCount &clipCount, const std::string &ref_string, bool isONT){

    ReadVariant *tmpReadResult = new ReadVariant();
    (*tmpReadResult).read_name = bam_get_qname(&aln);
    (*tmpReadResult).source_id = bamHdr.target_name[aln.core.tid];
    (*tmpReadResult).reference_start = aln.core.pos;
    (*tmpReadResult).is_reverse = bam_is_rev(&aln);
    
    // position relative to reference
    int ref_pos = aln.core.pos;

    // position relative to read
    int query_pos = 0;
    
    // Skip variants that are to the left of this read
    while( firstVariantIter != currentVariants->end() && (*firstVariantIter).first < ref_pos )
        firstVariantIter++;
    
    // Skip structure variants that are to the left of this read
    while( firstSVIter != currentSV->end() && (*firstSVIter).first < ref_pos )
        firstSVIter++;
    
    // Skip modify that are to the left of this read
    while( firstModIter != currentMod->end() && (*firstModIter).first < ref_pos )
        firstModIter++;

    // set variant start for current alignment
    std::map<int, RefAlt>::iterator currentVariantIter = firstVariantIter;
    
    // set structure variant start for current alignment
    std::map<int, std::map<std::string ,bool> >::iterator currentSVIter = firstSVIter;
     
    // set modify start for current alignment
    std::map<int, std::map<std::string ,RefAlt> >::iterator currentModIter = firstModIter;

    // set cigar pointer and number of cigar
    uint32_t *cigar = bam_get_cigar(&aln);
    int aln_core_n_cigar = aln.core.n_cigar;
    
    // reading cigar to detect varaint on this read
    for(int i = 0; i < aln_core_n_cigar ; i++ ){
        
        // get current cigar type and cigar length
        int cigar_op = bam_cigar_op(cigar[i]);
        int cigar_oplen = bam_cigar_oplen(cigar[i]);
        
        // get the starting position of each variant currently.
        int modPos = (*currentModIter).first;
        int svPos = (*currentSVIter).first;
        int variantPos = (*currentVariantIter).first;
        
        // get the first variant detected by the alignment.
        while( currentVariantIter != currentVariants->end() && variantPos < ref_pos ){
            currentVariantIter++;
            variantPos = (*currentVariantIter).first;
        }
        
        // Processing the region covered by the current CIGAR operator
        // Determine if any variant is included in the current CIGAR operator
        while( ( currentModIter     != currentMod->end()      && modPos     < ref_pos + cigar_oplen ) || 
               ( currentSVIter      != currentSV->end()       && svPos      < ref_pos + cigar_oplen ) || 
               ( currentVariantIter != currentVariants->end() && variantPos < ref_pos + cigar_oplen )){
            
            // modification's position is minimal
            if( ( currentVariantIter == currentVariants->end() || modPos < variantPos ) &&
                ( currentSVIter      == currentSV->end()       || modPos < svPos ) &&
                  currentModIter     != currentMod->end() ){
                
                // check this read contain modification
                std::map<std::string ,RefAlt>::iterator readIter = (*currentModIter).second.find(bam_get_qname(&aln));
                if( readIter != (*currentModIter).second.end() && modPos < variantPos ){
                    
                    // check varaint strand in vcf file is same as bam file
                    if( (*readIter).second.is_reverse == bam_is_rev(&aln) ){
                        int allele = ((*readIter).second.is_modify ? 0 : 1) ;
                        // using this quality to identify modification forward/reverse
                        int quality = (bam_is_rev(&aln) ? MOD_HET_REVERSE_STRAND : MOD_HET_FORWARD_STRAND);
                        Variant *tmpVariant = new Variant(modPos, allele, quality );
                        // push mod into result vector
                        (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                        delete tmpVariant;
                    }
                }
                currentModIter++;
                modPos = (*currentModIter).first;
            }
            // SV's position is minimal
            else if( ( currentVariantIter == currentVariants->end() || svPos < variantPos ) &&
                     ( currentModIter     == currentMod->end()      || svPos < modPos ) &&
                       currentSVIter      != currentSV->end()){
                        
                std::map<std::string ,bool>::iterator readIter = (*currentSV)[svPos].find(bam_get_qname(&aln));
                // If this read not contain SV, it means this read is the same as reference genome.
                // default this read the same as ref
                int allele = 0; 
                // this read contain SV.
                if( readIter != (*currentSV)[svPos].end() ){
                    allele = 1;
                }
                // use quality SV_HET to identify SVs
                // push this SV to vector
                Variant *tmpVariant = new Variant(svPos, allele, SV_HET );
                (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                delete tmpVariant;
                // next SV iter 
                currentSVIter++;
                svPos = (*currentSVIter).first;
            }

            // SNP's position is minimal
            else if( ( currentSVIter      == currentSV->end()  || variantPos < svPos ) &&
                     ( currentModIter     == currentMod->end() || variantPos < modPos ) &&
                       currentVariantIter != currentVariants->end() ){
                
                // CIGAR operators: MIDNSHP=X correspond 012345678
                // 0: alignment match (can be a sequence match or mismatch)
                // 7: sequence match
                // 8: sequence mismatch                
                if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
                    int refAlleleLen = (*currentVariantIter).second.Ref.length();
                    int altAlleleLen = (*currentVariantIter).second.Alt.length();
                    int offset = variantPos - ref_pos;
                    int base_q = UNDEFINED;
                    int allele = -1;
                    
                    // The position of the variant exceeds the length of the read.
                    if( query_pos + offset + 1 > int(aln.core.l_qseq) ){
                        return;
                    }

                    // SNP
                    if( refAlleleLen == 1 && altAlleleLen == 1){
                        char base = seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos + offset)];
                        if(base == (*currentVariantIter).second.Ref[0])
                            allele = 0;
                        else if(base == (*currentVariantIter).second.Alt[0])
                            allele = 1;

                        // using this quality to identify snp homozygous/heterozygous
                        base_q = (currentVariantIter->second.is_homozygous ? SNP_HOM : bam_get_qual(&aln)[query_pos + offset]);
                    } 
            
                    // insertion
                    //if( refAlleleLen == 1 && altAlleleLen != 1 && align.op[i+1] == 1 && i+1 < align.cigar_len){
                    if( refAlleleLen == 1 && altAlleleLen != 1 && i+1 < aln_core_n_cigar){
                
                        // currently, qseq conversion is not performed. Below is the old method for obtaining insertion sequence.
                
                        // uint8_t *qstring = bam_get_seq(aln); 
                        // qseq[i] = seq_nt16_str[bam_seqi(qstring,i)]; 
                        // std::string prevIns = ( align.op[i-1] == 1 ? qseq.substr(prev_query_pos, align.ol[i-1]) : "" );

                        if ( ref_pos + cigar_oplen - 1 == variantPos && bam_cigar_op(cigar[i+1]) == 1 ) {
                            allele = 1 ;
                        }
                        else {
                            allele = 0 ;
                        }
                        // using this quality to identify indel
                        base_q = (currentVariantIter->second.is_homozygous ? INDEL_HOM : INDEL_HET);
                    } 
            
                    // deletion
                    //if( refAlleleLen != 1 && altAlleleLen == 1 && align.op[i+1] == 2 && i+1 < align.cigar_len){
                    if( refAlleleLen != 1 && altAlleleLen == 1 && i+1 < aln_core_n_cigar) {

                        if ( ref_pos + cigar_oplen - 1 == variantPos && bam_cigar_op(cigar[i+1]) == 2 ) {
                            allele = 1 ;
                        }
                        else {
                            allele = 0 ;
                        }
                        // using this quality to identify indel
                        base_q = (currentVariantIter->second.is_homozygous ? INDEL_HOM : INDEL_HET);
                    } 
            
                    if( allele != -1 ){
                        // record snp result
                        Variant *tmpVariant = new Variant(variantPos, allele, base_q);
                        (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                        delete tmpVariant;                        
                    }
                    currentVariantIter++;
                    variantPos = (*currentVariantIter).first;
                }
                else break;
            }
        }
        
        // Preparing to process the next CIGAR operator.
        
        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
            query_pos += cigar_oplen;
            ref_pos += cigar_oplen;
        }
        // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += cigar_oplen;
        }
        else if( cigar_op == 2 ){
            
            // If a reference is given
            // it will determine whether the SNP falls in the homopolymer
            // and start the processing of the SNP fall in the alignment GAP
            if(ref_string != ""){
                int del_len = cigar_oplen;
                if ( ref_pos + del_len + 1 == (*currentVariantIter).first ){
                    //if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        // special case
                    //}
                }
                else if( (*currentVariantIter).first >= ref_pos  && (*currentVariantIter).first < ref_pos + del_len ){
                    // check snp in homopolymer
                    if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        
                        int refAlleleLen = (*currentVariantIter).second.Ref.length();
                        int altAlleleLen = (*currentVariantIter).second.Alt.length();
                        int base_q = UNDEFINED;
                        
                        if( query_pos + 1 > aln.core.l_qseq ){
                            return;
                        }

                        int allele = -1;
                        // SNP
                        if( refAlleleLen == 1 && altAlleleLen == 1){
                            // get the next match
                            char base = seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos)];
                            if( base == (*currentVariantIter).second.Ref[0] ){
                                allele = 0;
                            }
                            else if( base == (*currentVariantIter).second.Alt[0] ){
                                allele = 1;
                            }
                            base_q = (currentVariantIter->second.is_homozygous ? SNP_HOM : bam_get_qual(&aln)[query_pos]);
                        }
                        // the read deletion contain VCF's deletion
                        else if( refAlleleLen != 1 && altAlleleLen == 1 ){

                            if( refAlleleLen != 1 && altAlleleLen == 1){
                                //std::string delSeq = ref_string.substr(ref_pos - 1, align.ol[i] + 1);
                                //std::string refSeq = (*currentVariantIter).second.Ref;
                                //std::string altSeq = (*currentVariantIter).second.Alt;

                                allele = 1;
                                // using this quality to identify indel
                                base_q = (currentVariantIter->second.is_homozygous ? INDEL_HOM : INDEL_HET);
                            }
                            else if ( allele == -1 ) {
                                allele = 0;
                                // using this quality to identify indel
                                base_q = (currentVariantIter->second.is_homozygous ? INDEL_HOM : INDEL_HET);
                            }
                            
                        }
                        
                        if(allele != -1){
                            Variant *tmpVariant = new Variant((*currentVariantIter).first, allele, base_q);
                            (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                            currentVariantIter++;
                            delete tmpVariant;
                        }

                    }
                }
            }
            ref_pos += cigar_oplen;
        }
        // 3: skipped region from the reference
        else if( cigar_op == 3 ){
            ref_pos += cigar_oplen;
        }
        // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            query_pos += cigar_oplen;
            getClip(ref_pos, i, cigar_oplen, clipCount);
        }
        // 5: hard clipping (clipped sequences NOT present in SEQ)
        else if( cigar_op == 5 ){
            getClip(ref_pos, i, cigar_oplen, clipCount);
        }
        // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 6 ){
            // do nothing
        }
        else{
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
            exit(1);
        }
    }
    if( (*tmpReadResult).variantVec.size() > 0 )
        readVariantVec.push_back((*tmpReadResult));
    
    delete tmpReadResult;
}

void BamParser::getClip(int pos, int clipFrontBack, int len, ClipCount &clipCount){
    if (len > 5){
        if (clipFrontBack == FRONT){
            clipCount[pos][FRONT]++;
        }
        else {
            clipCount[pos][BACK]++;
        }
    }
}

METHParser::METHParser(PhasingParameters &in_params, SnpParser &in_snpFile, SVParser &in_svFile){
	params = &in_params;
    snpFile = &in_snpFile;
    svFile = &in_svFile;
    representativePos=-1;
    upMethPos = -1;
    
    chrVariant = new std::map<std::string, std::map<int, std::map<std::string ,RefAlt> > >;
    representativeMap = new std::map<int, int >;
    
    if( params->modFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(params->modFile);
    }
    else if( params->modFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(params->modFile);
    }
}

void METHParser::writeResult(ChrPhasingResult &chrPhasingResult){
    
    if( params->modFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressInput(params->modFile, params->resultPrefix+"_mod.vcf", chrPhasingResult);
    }
    else if( params->modFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressInput(params->modFile, params->resultPrefix+"_mod.vcf", chrPhasingResult);
    }
    return;
}

METHParser::~METHParser(){
	delete chrVariant;
    delete representativeMap;
}

void METHParser::parserProcess(std::string &input){
    // header
    if( input.substr(0, 1) != "#" ){
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 ){
            return;
        }
        
        // trans 1-base position to 0-base
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];
                
        // find GT flag
        int colon_pos = 0;
        int gt_pos = fields[8].find("GT");
        for(int i =0 ; i< gt_pos ; i++){
            if(fields[8][i]==':')
                colon_pos++;
        }
        
        // In a series of consecutive methylation positions, 
        // the first methylation position will be used as the representative after merging.
        if( upMethPos + 1 != pos ){
            representativePos = pos;
        }
        
        // find GT value start
        int current_colon = 0;
        int modify_start = 0;
        for(unsigned int i =0; i < fields[9].length() ; i++){
            if( current_colon >= colon_pos )
                break;
            if(fields[9][i]==':')
                current_colon++;  
            modify_start++;
        }
       
        // homo GT
        if(fields[9][modify_start]==fields[9][modify_start+2]){
            return;
        }
        
        // conflict pos with SNP and SV
        if( (*snpFile).findSNP(chr,pos) || (*svFile).findSV(chr,pos) ){
            return;
        }
        
        bool is_reverse;
        // get strand
        if(fields[7].find("RS=P") != std::string::npos){
            is_reverse = false;
        }
        else if(fields[7].find("RS=N") != std::string::npos){
            is_reverse = true;
        }
        else{
            return;
        }
        
        // parse MR and NR reads
        // get modification reads (MR) INFO
        int read_pos = fields[7].find("MR=");
        read_pos = fields[7].find("=",read_pos);
        read_pos++;
                
        int next_field = fields[7].find(";",read_pos);
        std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
        std::stringstream totalMonReadStream(totalRead);
        
        // Extract the reads in the MR string. These reads contain methylation.
        std::string read;
        while(std::getline(totalMonReadStream, read, ',')){
            RefAlt tmp;
            tmp.is_reverse = is_reverse;
            tmp.is_modify = true;
            (*chrVariant)[chr][representativePos][read]= tmp;
        }
        
        // get nonmodification reads (NR) INFO
        read_pos = fields[7].find("NR=");
        read_pos = fields[7].find("=",read_pos);
        read_pos++;
        
        next_field = fields[7].find(";",read_pos);
        totalRead = fields[7].substr(read_pos,next_field-read_pos);
        std::stringstream totalnonModReadStream(totalRead);
        
        // Extract the reads in the NR string. These reads not contain methylation.
        while(std::getline(totalnonModReadStream, read, ',')){
            RefAlt tmp;
            tmp.is_reverse = is_reverse;
            tmp.is_modify = false;
            (*chrVariant)[chr][representativePos][read]= tmp;
        }
        
        // Record the positions corresponding to the current representative position
        (*representativeMap)[pos]= representativePos;  
        upMethPos = pos;
    }
}

std::map<int, std::map<std::string ,RefAlt> > METHParser::getVariants(std::string chrName){
    std::map<int, std::map<std::string ,RefAlt> > targetVariants;
    std::map<std::string, std::map<int, std::map<std::string ,RefAlt> > >::iterator chrIter = chrVariant->find(chrName);
    
    if( chrIter != chrVariant->end() )
        targetVariants =  (*chrIter).second;
        
    return targetVariants;
}

