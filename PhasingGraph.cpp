#include "PhasingGraph.h"

//SubEdge

SubEdge::SubEdge():readCount(0){ 
    refRead = new std::map<int, std::vector<std::string> >;
    altRead = new std::map<int, std::vector<std::string> >;
    refQuality = new std::map<int, int>;
    altQuality = new std::map<int, int>;
    refReadCount = new std::map<int, float>;
    altReadCount = new std::map<int, float>;
}

SubEdge::~SubEdge(){ 
}

void SubEdge::destroy(){
    delete refRead;
    delete altRead;
    delete refQuality;
    delete altQuality;
    delete refReadCount;
    delete altReadCount;
}

void SubEdge::addSubEdge(Variant &currentNode, Variant &connectNode, std::string readName, int conditionQuality, double lowQualityWeight){
    double edgeWeight = 1;
    if( (currentNode.quality > UNDEFINED && currentNode.quality < conditionQuality ) || (connectNode.quality > UNDEFINED && connectNode.quality < conditionQuality ) ){
        edgeWeight = lowQualityWeight;
    }

    // target noded is REF allele
    if(connectNode.allele == 0 ){
        // debug, this parameter will record the names of all reads between two podoubles
        //(*refRead)[connectNode.position].push_back(readName);
        // quality sum
        // (*refQuality)[connectNode.position] += currentQuality + connectNode.quality;
        (*refReadCount)[connectNode.position] += edgeWeight;
    }
    // target noded is ALT allele
    else if(connectNode.allele == 1 ){
        // debug, this parameter will record the names of all reads between two points
        // (*altRead)[connectNode.position].push_back(readName);
        // quality sum
        // (*refQuality)[connectNode.position] += currentQuality + connectNode.quality;
        (*altReadCount)[connectNode.position] += edgeWeight;
    }
    // readCount++;
}

std::pair<float,float> SubEdge::BestPair(int targetPos){
    return std::make_pair( getRefReadCount(targetPos), getAltReadCount(targetPos) );
}

float SubEdge::getRefReadCount(int targetPos){
    std::map<int, float>::iterator posIter = refReadCount->find(targetPos);
    if( posIter != refReadCount->end() ){
        return (*refReadCount)[targetPos];
    }
    return 0;
}

float SubEdge::getAltReadCount(int targetPos){
    std::map<int, float>::iterator posIter = altReadCount->find(targetPos);
    if( posIter != altReadCount->end() ){
        return (*altReadCount)[targetPos];
    }
    return 0;
}

std::vector<std::string> SubEdge::showEdge(std::string message){
    std::vector<std::string> result;
    for(std::map<int, float >::iterator edgeIter = refReadCount->begin() ; edgeIter != refReadCount->end() ; edgeIter++ ){
        result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    for(std::map<int, float >::iterator edgeIter = altReadCount->begin() ; edgeIter != altReadCount->end() ; edgeIter++ ){
        result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    return result;
}

std::vector<std::pair<int,int>> SubEdge::getConnectPos(){
    std::vector<std::pair<int,int>> result;
    for(std::map<int, float >::iterator edgeIter = refReadCount->begin() ; edgeIter != refReadCount->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 0 ) );
    }
    for(std::map<int, float >::iterator edgeIter = altReadCount->begin() ; edgeIter != altReadCount->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 1 ) );
    }
    return result;
}

int SubEdge::getQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality->find(targetPos.first);
        if( qIter == refQuality->end() )
            return 0;
        else
            return (*refQuality)[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality->find(targetPos.first);
        if( qIter == altQuality->end() )
            return 0;
        else
            return (*altQuality)[targetPos.first];
    }
    return 0;
}

int SubEdge::getAvgQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality->find(targetPos.first);
        if( qIter == refQuality->end() )
            return 0;
        else
            return (*refQuality)[targetPos.first]/(*refReadCount)[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality->find(targetPos.first);
        if( qIter == altQuality->end() )
            return 0;
        else
            return (*altQuality)[targetPos.first]/(*altReadCount)[targetPos.first];
    }
    return 0;
}

VariantEdge::VariantEdge(int inCurrPos){
    currPos = inCurrPos;
    alt = new SubEdge();
    ref = new SubEdge();
}

//VariantEdge
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(int targetPos, bool isONT, double edgeThreshold, bool debug){
    std::pair<float,float> refBestPair  = ref->BestPair(targetPos);
    std::pair<float,float> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    float rr = refBestPair.first;
    float ra = refBestPair.second;
    float ar = altBestPair.first;
    float aa = altBestPair.second;
    // initialize the edge connection
    // -1 : not connect
    int refAllele = -1;
    int altAllele = -1;
    
    double edgeSimilarRatio = (double)std::min((rr+aa),(ar+ra)) / (double)std::max((rr+aa),(ar+ra));

    if( rr + aa > ra + ar ){
        // RR conect
        refAllele = 1;
        altAllele = 2;
    }
    else if( rr + aa < ra + ar ){
        // RA connect
        refAllele = 2;
        altAllele = 1;
    }
    else if( rr + aa == ra + ar ){
        // no connect 
        // not sure which is better
    }
    
    if( edgeSimilarRatio > edgeThreshold ){
        refAllele = -1;
        altAllele = -1;
    }
    
    if(debug){
        std::cout<< currPos << "\t->\t" << targetPos << "\t|rr aa | ra ar\t" << "\t" << rr << "\t" << aa << "\t" << ra << "\t" << ar  << "\n";
    }
    
    // create edge pairs
    PosAllele refEdge = std::make_pair( targetPos, refAllele );
    PosAllele altEdge = std::make_pair( targetPos, altAllele );
    // return edge pair
    return std::make_pair( refEdge, altEdge );
}

std::pair<int,int> VariantEdge::findNumberOfRead(int targetPos){
    std::pair<int,int> refBestPair  = ref->BestPair(targetPos);
    std::pair<int,int> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    int rr = refBestPair.first;
    int ra = refBestPair.second;
    int ar = altBestPair.first;
    int aa = altBestPair.second;
    return std::make_pair( rr + aa , ra +ar );
}

//BlockRead
void BlockRead::recordRead(std::string readName){
    std::map<std::string,int>::iterator readIter = readVec.find(readName);
    if( readIter == readVec.end() )
        readVec[readName] = 1;
    else
        readVec[readName]++;
}

//VairiantGraph
void VairiantGraph::edgeConnectResult(){
    // current snp, haplotype (1 or 2), support snp
    std::map<int, std::map<int,std::vector<int> > > *hpCountMap = new std::map<int, std::map<int,std::vector<int> > >;
    
    int blockStart = -1;
    int currPos = -1;
    int nextPos = -1;
    int lastConnectPos = -1;
    
    // Visit all position and assign SNPs to haplotype.
    // Avoid recording duplicate information,
    // only one of the two alleles needs to be used for each SNP
    for(auto variantIter = variantPosType->begin() ; variantIter != variantPosType->end() ; variantIter++ ){
        // check next position
        auto nextNodeIter = std::next(variantIter, 1);
        if( nextNodeIter == variantPosType->end() ){
             break;
        }
        
        currPos = variantIter->first;
        nextPos = nextNodeIter->first;
        
        // There should not be a large distance between any two variants, 
        // with the default being a distance of 300000bp, equivalent to one centromere length.
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        // get the number of HP1 and HP2 supported reference allele
        int h1 = (*hpCountMap)[currPos][HAPLOTYPE1].size();
        int h2 = (*hpCountMap)[currPos][HAPLOTYPE2].size();

        // new block, set this position as block start 
        if( h1 == h2 ){
            // No new blocks should be created if the next SNP has already been picked up
            if( currPos < lastConnectPos ){
                continue;
            }
            // check block size to remove one node island
            if(!posPhasingResult->empty() && posPhasingResult->rbegin()->first == blockStart){
                posPhasingResult->erase(blockStart);
            }

            blockStart = currPos;
            posPhasingResult->emplace(currPos, PhasingResult(HAPLOTYPE1, blockStart));
        }
        else{
            if( h1 > h2 || h1 < h2 ){
                int currHP = ( h1 > h2 ? HAPLOTYPE1 : HAPLOTYPE2 );
                posPhasingResult->emplace(currPos, PhasingResult(currHP, blockStart));
            }
        }

        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( currPos );
        if( edgeIter==edgeList->end() ){
            continue;
        }

        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent ; i++ ){
            if( nextNodeIter->second == SNP_HOM ){
                i--;
            }
            else{
                // consider reads from the currnt SNP and the next (i+1)'s SNP
                std::pair<PosAllele,PosAllele> tmp = edgeIter->second->findBestEdgePair(nextNodeIter->first, params->isONT, params->edgeThreshold, false);
                // -1 : no connect  
                //  1 : the haplotype of next (i+1)'s SNP are same as previous
                //  2 : the haplotype of next (i+1)'s SNP are different as previous
                if( tmp.first.second != -1 ){
                    // record the haplotype resut of next (i+1)'s SNP
                    if( (*posPhasingResult)[currPos].refHaplotype == HAPLOTYPE1 ){
                        if( tmp.first.second == 1 ){
                            (*hpCountMap)[nextNodeIter->first][HAPLOTYPE1].push_back(currPos);
                        }
                        if( tmp.first.second == 2 ){
                            (*hpCountMap)[nextNodeIter->first][HAPLOTYPE2].push_back(currPos);
                        }
                    }
                    if( (*posPhasingResult)[currPos].refHaplotype == HAPLOTYPE2 ){
                        if( tmp.first.second == 1 ){
                            (*hpCountMap)[nextNodeIter->first][HAPLOTYPE2].push_back(currPos);
                        }
                        if( tmp.first.second == 2 ){
                            (*hpCountMap)[nextNodeIter->first][HAPLOTYPE1].push_back(currPos);
                        }
                    }
                    if( params->generateDot ){
                        std::string e1 = std::to_string(currPos+1) + ".1\t->\t" + std::to_string(tmp.first.first+1) + "." + std::to_string(tmp.first.second);
                        std::string e2 = std::to_string(currPos+1) + ".2\t->\t" + std::to_string(tmp.second.first+1) + "." + std::to_string(tmp.second.second);

                        dotResult.push_back(e1);
                        dotResult.push_back(e2);
                    }
                    
                    lastConnectPos = nextNodeIter->first;
                }
            }
            nextNodeIter++;
            if( nextNodeIter == variantPosType->end() ){
                break;
            }
        }
    }

    delete hpCountMap;
}

VairiantGraph::VairiantGraph(std::string &in_ref, PhasingParameters &in_params){
    params=&in_params;
    ref=&in_ref;
    
    variantPosType = new std::map<int,int>;
    edgeList = new std::map<int,VariantEdge*>;
    readHpMap = new std::map<std::string,int>;
}

VairiantGraph::~VairiantGraph(){
}

void VairiantGraph::destroy(){
    dotResult.clear();
    dotResult.shrink_to_fit();

    for( auto edgeIter = edgeList->begin() ; edgeIter != edgeList->end() ; edgeIter++ ){
        edgeIter->second->ref->destroy();
        edgeIter->second->alt->destroy();
        delete edgeIter->second->ref;
        delete edgeIter->second->alt;
    }
    
    delete variantPosType;
    delete edgeList;
    delete readHpMap;
}

void VairiantGraph::addEdge(std::vector<ReadVariant> &in_readVariant){
    readVariant = &in_readVariant;
    std::map<std::string,ReadVariant> mergeReadMap;

    // each read will record fist and list variant posistion
    std::map<std::string, std::pair<int,int>> alignRange;
    // record an iterator for all alignments of a read.
    std::map<std::string, std::vector<int>> readIdxVec;
    // record need del read index
    std::vector<int> delReadIdx;

    // Check for overlaps among different alignments of a read and filter out the shorter overlapping alignments.
    for(int readIter = 0 ; readIter < (int)in_readVariant.size() ; readIter++ ){
        std::string readName = in_readVariant[readIter].read_name;
        int firstVariantPos = in_readVariant[readIter].variantVec[0].position;
        int lastVariantPos  = in_readVariant[readIter].variantVec[in_readVariant[readIter].variantVec.size()-1].position;
        
        auto rangeIter = alignRange.find(readName);
        // this read name appears for the first time
        if( rangeIter == alignRange.end() ){
            alignRange[readName]=std::make_pair(firstVariantPos,lastVariantPos);
        }
        // the read appears more than once, check if the alignments overlap
        else{
            // overlap
            if( alignRange[readName].first <= firstVariantPos && firstVariantPos <= alignRange[readName].second ){
                double alignStart   = std::min(alignRange[readName].first, firstVariantPos);
                double alignEnd     = std::max(alignRange[readName].second, lastVariantPos);
                double alignSpan    = alignEnd - alignStart + 1;
                double overlapStart = std::max(alignRange[readName].first, firstVariantPos);
                double overlapEnd   = std::min(alignRange[readName].second, lastVariantPos);
                double overlapLen   = overlapEnd - overlapStart + 1;
                double overlapRatio = overlapLen / alignSpan;
                
                //filtering highly overlapping alignments.
                if( overlapRatio >= params->overlapThreshold ){
                    int alignLen1 = alignRange[readName].second - alignRange[readName].first + 1;
                    int alignLen2 = lastVariantPos - firstVariantPos + 1;
                    
                    // filter shorter alignment
                    // current alignment is shorter
                    if( alignLen2 <= alignLen1 ){
                        delReadIdx.push_back(readIter);
                    }
                    // previous alignment is shorter
                    else{
                        // iterate all previous alignments
                        for(int iter = 0 ; iter < (int)readIdxVec[readName].size() ; iter++ ){
                            delReadIdx.push_back(readIdxVec[readName][iter]);
                        }

                        // update range
                        alignRange[readName].first  = firstVariantPos;
                        alignRange[readName].second = lastVariantPos;
                        readIdxVec[readName].clear();
                        readIdxVec[readName].push_back(readIter);
                    }
                    continue;
                }
            }
            // update range
            alignRange[readName].second = lastVariantPos;
        }
        readIdxVec[readName].push_back(readIter);
    }

    // sort read index
    std::sort(delReadIdx.begin(), delReadIdx.end());
    // remove overlap alignment
    for( int idx = delReadIdx.size() -1 ; idx > 0 ; idx-- ){
        in_readVariant.erase( in_readVariant.begin() + delReadIdx[idx] );
    }

    int readCount=0;
    // merge alignment
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
        // Creating a pseudo read which allows filtering out variants that should not be phased
        //ReadVariant tmpRead;
        // Visiting all the variants on the read
        for( auto variant : (*readIter).variantVec ){
            readCount++;
            if( variant.quality <= UNDEFINED ){
                if (variant.quality != SNP_HOM){
                    (*variantPosType)[variant.position] = variant.quality;
                }
            }
            else{
                (*variantPosType)[variant.position] = SNP_HET;
            }

            mergeReadMap[(*readIter).read_name].variantVec.push_back(variant);
        }
    }    
    
    for(std::map<std::string,ReadVariant>::iterator readIter = mergeReadMap.begin() ; readIter != mergeReadMap.end() ; readIter++){
        (*readIter).second.sort();
        
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = readIter->second.variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        
        while(variant1Iter != readIter->second.variantVec.end() && variant2Iter != readIter->second.variantVec.end() ){
            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = edgeList->find(variant1Iter->position);
            if( posIter == edgeList->end() )
                (*edgeList)[variant1Iter->position] = new VariantEdge(variant1Iter->position);

            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent; nextNode++){
                if (variant2Iter->quality == SNP_HOM && variant1Iter->quality != SNP_HOM){
                    nextNode--;
                }
                // this allele support ref
                if( variant1Iter->allele == 0 )
                    (*edgeList)[variant1Iter->position]->ref->addSubEdge((*variant1Iter), (*variant2Iter),readIter->first,params->baseQuality,params->edgeWeight);
                // this allele support alt
                if( (*variant1Iter).allele == 1 )
                    (*edgeList)[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter), (*variant2Iter),readIter->first,params->baseQuality,params->edgeWeight);
                
                // next snp
                variant2Iter++;
                if( variant2Iter == readIter->second.variantVec.end() ){
                    break;
                }
            }

            variant1Iter++;
            variant2Iter = std::next(variant1Iter,1);
        }
    }

} 

void VairiantGraph::readCorrection(){
    
    
    std::map<std::string,std::map<int,std::map<int,int>>> readBlockHP;
    //
    std::map<std::string,std::map<int,std::map<int,int>>> readBlockHPcount;
    
    
    // haplotype, <position <allele, base count>>
    std::map<int,std::map<int,std::map<int,int>>> *hpAlleleCountMap = new std::map<int,std::map<int,std::map<int,int>>>;

    // phasing result and variant allele mapping
    const int variantHaplotype[2][2] = {
        {HAPLOTYPE1, HAPLOTYPE2},
        {HAPLOTYPE2, HAPLOTYPE1}
    };

    // iter all read, determine the haplotype of the read
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        double haplotype1Count = 0;
        double haplotype2Count = 0;
        //int block;
        
        // loop all variant 
        for( auto variant : (*readIter).variantVec ){
            auto posPhasingResultIter = posPhasingResult->find(variant.position);
            if( posPhasingResultIter != posPhasingResult->end() ){
                const PhasingResult& phasingResult = posPhasingResultIter->second;
                if(variantHaplotype[phasingResult.refHaplotype][variant.allele] == HAPLOTYPE1)haplotype1Count++;
                else haplotype2Count++;
            }
        }
        
        // tag high confident reads
        if( std::max(haplotype1Count,haplotype2Count)/(haplotype1Count+haplotype2Count) > params->readConfidence && (haplotype1Count + haplotype2Count) > 1 ){
            // tag read with the corresponding haplotype
            int belongHP = ( haplotype1Count > haplotype2Count ? 0 : 1 );
            (*readHpMap)[(*readIter).read_name] = belongHP;
            
            //readBlockHP[(*readIter).read_name][(*readIter).reference_start][block]=belongHP;
            //readBlockHPcount[(*readIter).read_name][block][belongHP]++;
            
            for(auto variantIter = (*readIter).variantVec.begin() ; variantIter != (*readIter).variantVec.end() ; variantIter++ ){
                if( (*variantIter).allele == 0 || (*variantIter).allele == 1){
                    (*hpAlleleCountMap)[belongHP][(*variantIter).position][(*variantIter).allele]++;
                }
            }
        }
        else{
            (*readHpMap)[(*readIter).read_name] = -1;
        }
    }
    
    /*
    for(auto readIter = readBlockHP.begin() ; readIter != readBlockHP.end() ; readIter++ ){
        
        int max = 0;
        for(auto blockIter = readBlockHPcount[readIter->first].begin() ; blockIter != readBlockHPcount[readIter->first].end() ; blockIter++ ){
            if( (*blockIter).second.size() > max ){
                max = (*blockIter).second.size();
            }
        }

        std::cout<< readIter->first << "\t" << max;
        
        for(auto startIter = readIter->second.begin() ; startIter != readIter->second.end() ; startIter++ ){
            std::cout<< "\t" << startIter->first << ",";
            for(auto blockIter = startIter->second.begin() ; blockIter != startIter->second.end() ; blockIter++ ){
                std::cout<< blockIter->first << "," <<  blockIter->second;
            }
        }
        std::cout<< "\n";
    }
    */
    
    double snpConfidenceThreshold = params->snpConfidence;
    std::map<int,std::map<int,int>> hpAllele;
    // reassign allele result
    for(auto variantIter = variantPosType->begin() ; variantIter != variantPosType->end() ; variantIter++ ){
        int position = variantIter->first;
        auto posPhasingResultIter = posPhasingResult->find(position);
        if (posPhasingResultIter == posPhasingResult->end()) {
            continue;
        }
        
        double hp1Ref = (*hpAlleleCountMap)[0][position][0];
        double hp1Alt = (*hpAlleleCountMap)[0][position][1];
        double hp2Ref = (*hpAlleleCountMap)[1][position][0];
        double hp2Alt = (*hpAlleleCountMap)[1][position][1];
        double result1reads = hp1Ref + hp2Alt;
        double result2reads = hp2Ref + hp1Alt;
        double resultConfidence = std::max(result1reads, result2reads) / (result1reads + result2reads);
        
        int hp1Result = -1;
        int hp2Result = -1;
        
        if( resultConfidence > snpConfidenceThreshold ){
            if( result1reads > result2reads ){
                hp1Result = HAPLOTYPE1;
                hp2Result = HAPLOTYPE2;
            }
            else if( result1reads < result2reads ){
                hp1Result = HAPLOTYPE2;
                hp2Result = HAPLOTYPE1;
            }
        }

        if( hp1Result != -1 && hp2Result != -1 ){
            posPhasingResultIter->second.refHaplotype = hp1Result;
        }
        else{
            posPhasingResult->erase(posPhasingResultIter);
        }
    }

    delete hpAlleleCountMap;
}

void VairiantGraph::writingDotFile(std::string dotPrefix){
    
    std::ofstream resultVcf(dotPrefix+".dot");

    if(!resultVcf.is_open()){
        std::cerr<< "Fail to open write file: " << dotPrefix+".vcf" << "\n";
    }
    else{
        resultVcf << "digraph G {\n";

        for(auto edge : dotResult){
            resultVcf << edge << "\n";
        }
        resultVcf << "}\n";
    }
    return;
}

std::map<std::string,int>* VairiantGraph::getReadHP(){
    return readHpMap;
}

int VairiantGraph::totalNode(){
    return variantPosType->size();
}

void VairiantGraph::phasingProcess(PosPhasingResult &inPosPhasingResult){
    posPhasingResult = &inPosPhasingResult;

    // This step involves converting all reads into a graph structure, which will be stored as an edge list
    // in a two-layer map. The first layer of the map uses the starting coordinate as the key and contains
    // a second layer map as the value. The second layer map uses the destination coordinate as the key and
    // stores the number of support read as values. (There is another map used for debugging purposes that
    // treats the read name vector as a value.) The method begins by visiting the coordinates covered by each 
    // read and recording this information in 'variantPosType.' Subsequently, it connects the coordinates contained 
    // in each read on the graph. Specifically, each coordinate is connected to the next N coordinates in a 
    // linear fashion.
    this->edgeConnectResult();

    // This step will utilize the results of graph phasing to attempt to separate all the reads into two 
    // haplotypes and then identify high-confidence SNPs using reads from the two distinct haplotypes.
    this->readCorrection();  
}

Clip::Clip(std::string &chr, ClipCount &clipCount){
    this->chr = chr;
    getCNVInterval(clipCount);
}

Clip::~Clip(){
}

void Clip::getCNVInterval(ClipCount &clipCount){
    bool push=0;
    bool up =1;
    bool convert =0;

    int push_pos=-2;
    int firstClipOldPos=-10000;
    int secondClipOldPos=-10000;

    int upCount = 0;
    int downCount = 0;
    int passDept=5;
    int intervalSecondMax=0;
    int tmp=0;

    for(auto posIter = clipCount.begin(); posIter != clipCount.end() ; posIter++ ){
        upCount = posIter->second[FRONT];
        downCount = posIter->second[BACK];
        if (upCount>=5 || downCount>=5){
            if (((posIter->first-secondClipOldPos)>10000)){
                if (push_pos!=-1 && push_pos!=-2){
                    firstClipOldPos=push_pos;
                    cnvInterval.push_back(push_pos);
                    push_pos=-2;
                }
                else if (push_pos==-1){
                    push_pos=-2;
                }
                push=0;
                convert=0;
                passDept=5;
                intervalSecondMax=0;
                if((up && upCount<downCount)||(!up && upCount>downCount)){
                    up=!up;
                }
            }
            if (upCount>=passDept || downCount>=passDept){
                if(!convert&&((up && upCount<downCount)||(!up && upCount>downCount))){
                    up=!up;
                    convert=1;
                }
                if (abs(upCount - downCount) > 12){
                    if(((up && upCount<downCount)||(!up && upCount>downCount))){
                        up=!up;
                    }
                    if (((posIter->first-firstClipOldPos)>10000) && !push){
                        push_pos=posIter->first;
                        passDept=(upCount-downCount>10?upCount:downCount)/4;
                        passDept=(passDept>5?passDept:5);
                        push=1;
                    }
                    else if (((push&&convert)||(!push&&!convert))){
                        tmp= up ? upCount : downCount;
                        if (intervalSecondMax < tmp) {
                            intervalSecondMax = tmp;
                        }
                    }
                }
                else if (((push&&convert)||(!push&&!convert))){
                    tmp= up ? upCount : downCount;
                    if (intervalSecondMax < tmp) {
                        intervalSecondMax = tmp;
                    }
                }
                if ((push && convert && passDept<=intervalSecondMax)){
                    push_pos=-1;
                }
                secondClipOldPos=posIter->first;
            }
        }
    }
}
