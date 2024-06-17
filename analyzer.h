#ifndef STK_ANALYZER_H
#define STK_ANALYZER_H
#include "param.h"
using namespace slib;
using namespace slib::sbio;
using namespace slib::sapp;
namespace stk {
	/**
	*
	*/
	class Analyzer {
	public:
		Param par;
	public:
		Analyzer();
		Analyzer(SDictionary& pref);
		~Analyzer();
		Response analyze();
		void setParam(SDictionary& pref);
		void reset();
	};
	/**
	*
	*/
	class ReadCounter {
		int bin;
		Array<svecf>* depths;
	public:
		ReadCounter();
		ReadCounter(NGSData *summary, Param *par);
		~ReadCounter();
		void set(NGSData* summary, Param* par);
		void count(sbpos *pos, CigarArray* cigars);
		void count(sbam::ReadInfo* read);
	};	
	//
	class BamReader;
	/**
	*
	*/
	class SRAnalysis {
		friend BamReader;

		String _buffer;
	public:
		NGSData* summary;
		Param* par;
		
		Status* status;
		SWork* threads;
		SLogger* logger;
		
		DNASeqTrie2 trie;
		SeqSearch search;

		RecycleArray<SVar> svars;

	public:
		SRAnalysis(Analyzer* an = nullptr);
		SRAnalysis(NGSData* summary, stk::Param *p);
		~SRAnalysis();
		void setReader(stk::BamReader* br);
		void splitreads(NGSData* summary, IOStream& stream);

		void addClipRead(sbam::ReadInfo* ri, DIRECTION clip);
		void addClipReadPair(sbam::ReadInfo* ri1, sbam::ReadInfo* ri2, bool inv = false);
		void realign();

		//void detect(BamFile* bam, NGSData* summary, Array<uintegerarray>* clips);
		
		void setData(NGSData* d);
		void setParam(Param* p);
		void reset();
	};
	/*

	*/
	class BamReader {
		friend SRAnalysis;
	public:
		Param* par;
		Status* status;
		SLogger* logger;
		SWork* threads;		
	public:
		BamReader(Analyzer *an);
		~BamReader();
		void readinfo(BamFile* bam, IOStream& stream);
		void summary1(int i, BamFile* bam, NGSData* summary, SRAnalysis *svd, vchunk range, bool async);
		void summary2_(int i, BamFile* bam, NGSData* summary, Array<Map<String, sbam::ReadInfo>>* pairs, SRAnalysis* svd, vchunk range, bool async);
		void summary2(BamFile* bam, NGSData* summary, Array<Map<String, sbam::ReadInfo>> *pairs, SRAnalysis* svd, vchunk range);
		void summarize(BamFile *bam, NGSData * summary);
		void totalize(NGSData* summery);
	};
	/*
	*
	*/
	class CNAnalysis {
		Param* par;
		Status* status;
		SLogger* logger;
		NGSData* smpl, * ctrl;
	public:
		Array<svecf> normalized[2], copies;
	public:
		CNAnalysis();
		CNAnalysis(Analyzer *an);
		~CNAnalysis();
		
		void setParam(stk::Param* p);
		void setData(NGSData* data, NGSData *control = nullptr);

		void depth(smath::Vector<svecf>& values);
		void copynum(smath::Vector<svecf>& values);
		
		//float count(const sbpos& pos, NGSData* data = nullptr);
		float count(int idx, const srange& range, bool control = false);
		float ncount(int idx, const srange& range, bool control = false);
		float copy(int idx, const srange& range);
		//void copies(svecf& vec, int idx, const srange& range);
		//subyte check();
		void analyze(CNData& cnd, int r, srange range);
		void analyze();
		void reset();

	};

	extern void makeVar1Index(const size_t length, Array<ArrayIterator<SVar>> *index, Array<SVar> *variants);
	extern void makeVar2Index(const size_t length, Array<ArrayIterator<SVar>>* findex, Array<ArrayIterator<SVar>>* rindex, Array<SVar>* variants);
	extern void makeVar3Index(const sveci* length, Array<ArrayIterator<SVar>>* tindex, Array<SVar>* variants);
	extern void makeVar4Index(const sveci* length, Array<ArrayIterator<SVar>>* tindex, Array<SVar>* variants);

	extern void selectDelCandidate(int t, Array<SVar*>* candidate, Array<SVar>* variants, stk::Param* par);
	extern void selectDupCandidate(int t, Array<SVar*>* candidate, Array<SVar>* variants, stk::Param* par);
	extern void selectInsCandidate(int t, Array<SVar*>* candidate, Array<SVar>* variants, stk::Param* par);
	extern void selectComplexCandidate(int t, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants1, Array<SVar>* variants2, stk::Param* vpar);
	extern void selectInvCandidate(int t, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants, stk::Param* vpar);
	extern void selectTrsCandidate(int t, int i1, int i2, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants1, Array<SVar>* variants2, stk::Param* vpar);
	extern void selectTrInvCandidate(int t, int i1, int i2, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants, stk::Param* vpar);

	extern void getDelCpy(Array<SPointer<sbio::Variant>> *primary, Array<SVar*>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getDupCpy(Array<SPointer<sbio::Variant>>* primary, Array<SVar*>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getInsCpy(Array<SPointer<sbio::Variant>>* primary, Array<SVar*>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getComplexCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getInvCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getInvInsCpy1(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getInvInsCpy2(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getTrsCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getTrsInsCpy1(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);	
	extern void getTrsInsCpy2(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getTrsInvCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getTrsInvInsCpy1(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);
	extern void getTrsInvInsCpy2(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par);

	class VarSearch {
		Param* par;
		SWork* threads;
		Status* status;
		SLogger* logger;

		stk::CNAnalysis cna;
		Array<Array<SVar*>> candidates1;
		Array<Array<Pair<SVar*, SVar*>>> candidates2;
		Array<Array<SPointer<Variant>>> primary;

		VarFilter filter;
	public:
		VarList variants;
	private:
		void selectCandidate(NGSData* data);
		void copyCheck(NGSData* data);
		void conflictCheck();
		void annotate();
	public:
		VarSearch(Analyzer *a);
		~VarSearch();

		void detect(NGSData* data);
		void reset();
	};
	
}
/*








#define SHOW_READS 0x01
#define BAM_SUMMERIZE 0x02
#define SPLIT_LIST 0x04
#define CPNUM_ARRAY 0x08
#define SEARCH_VAR 0x10

#define INSUFFICIENT_ERROR 0x0001
#define FILE_NOT_FOUND_ERROR 0x0002
#define SEQENCE_SUMMARY_FILE_ERROR 0x0002
#define VARIANT_FILE_ERROR 0x0004
#define DATABASE_ERROR 0x0008
#define TARGET_FILE_ERROR 0x0010

#include "sbioinfo.h"
#include "param.h"
#include "util.h"

using namespace slib;
using namespace slib::sbio;
//using namespace slib::sapp;

namespace stk {

#define ANALYSIS_COMPLETED 0xffff	
	
	class Analyzer {
	public:
		Param par;
		Status status;
		//SWork threads;
		//SBam bam;
		//SNGSData data;
		//SVarList vlist;
		//SLogger logger;

	public:
		Analyzer(SDictionary &pref);
		~Analyzer();

		void summerize();

		void analyze();
		void setParam(SDictionary& pref);
		void init();
	};
	*/
	/*
	extern inline bool isClipped(SCigar& cigar, sint& min) {
		return cigar.option == SCigar::SCLIP && min < cigar.length;
	}
	extern inline void countRead(sbpos *pos, SCigarArray *cigars, SNGSData* data) {
		auto& bin = data->summary.bin;
		sint rpos = pos->begin;
		auto depth = data->depth.count[pos->idx].ptr(rpos / bin);
		sforeach(*cigars) {
			sint beg = rpos / bin, end = (rpos + E_.length - 1) / bin;
			if (E_.option == SCigar::MATCH || E_.option == SCigar::PMATCH || E_.option == SCigar::MMATCH) {
				sforin(p, beg, end) depth[p - beg] += 1.0f;
				depth[0] -= (float)(rpos % bin) / bin;
				depth[end - beg] += (float)(rpos + E_.length - end * bin) / bin;
				rpos += E_.length;
				depth += (rpos / bin) - beg;
			}
			else if (E_.option == SCigar::DELETION || E_.option == SCigar::SKIP) {
				depth[0] += (float)(rpos % bin) / bin;
				depth[end - beg] += 1.0 - (float)(rpos + E_.length - end * bin) / bin;
				rpos += E_.length;
				depth += (rpos / bin) - beg;
			}
		}
	}
	*/
	

	


	
	/*
	extern inline double depthSum(sint r, srange range, SNGSData* data) {
		double val = 0.0;
		if (range.begin < 0) range.begin = 0;
		if (data->summary.reflen[r] <= range.end) range.end = data->summary.reflen[r] - 1;
		sint init = range.begin / data->summary.bin, term = range.end / data->summary.bin;
		float* dp = data->depth.count[r].ptr(init);
		if (init == term) return (double)(*dp) * range.length(true);
		else {
			val = (*dp) * (data->summary.bin - (range.begin % data->summary.bin));
			++init; ++dp;
			while (init < term) {
				val += (double)(*dp) * data->summary.bin; ++init; ++dp;
			}
			val += (double)(*dp) * ((range.end % data->summary.bin) + 1);
			return val;
		}
	}
	extern inline double depthValue(sint r, srange range, SNGSData* data) {
		double val = 0.0;
		if (range.begin < 0) range.begin = 0;
		if (data->summary.reflen[r] <= range.end) range.end = data->summary.reflen[r] - 1;
		sint init = range.begin / data->summary.bin, term = range.end / data->summary.bin;
		float* dp = data->depth.count[r].ptr(init);
		if (init == term) return (double)(*dp);
		else {
			val = (*dp) * (data->summary.bin - (range.begin % data->summary.bin));
			++init; ++dp;
			while (init < term) {
				val += (double)(*dp) * data->summary.bin; ++init; ++dp;
			}
			val += (double)(*dp) * ((range.end % data->summary.bin) + 1);
			return val / range.length(true);
		}
	}
	extern inline void copyValues(sint r, srange range, CNAnalysis* cna, double* values) {
		if (range.begin < 0) range.begin = 0;
		if (cna->data->summary.reflen[r] <= range.end) range.end = cna->data->summary.reflen[r] - 1;
		memset(values, 0, 2 * sizeof(double));
		auto len = range.length(true);
		if (len < MAX_CPCOUNT_SIZE) {
			values[0] = depthSum(r, range, cna->data);
			values[1] = cna->control ? depthSum(r, range, cna->control) : len;
		}
		else {
			sforeach(cna->cnvs[r]) {
				if (E_.pos.include(range)) {
					if (2 * len < E_.pos.length(true)) {
						values[0] = depthSum(r, range, cna->data);
						values[1] = cna->control ? depthSum(r, range, cna->control) : len;
					}
					else {
						values[0] = E_.copy[0];
						values[1] = E_.copy[1];
						if (E_.pos.begin < range.begin) {
							srange r1(E_.pos.begin, range.begin - 1);
							values[0] -= depthSum(r, r1, cna->data);
							values[1] -= cna->control ? depthSum(r, r1, cna->control) : r1.length(true);
						}
						if (range.end < E_.pos.end) {
							srange r2(range.end, E_.pos.end - 1);
							values[0] -= depthSum(r, r2, cna->data);
							values[1] -= cna->control ? depthSum(r, r2, cna->control) : r2.length(true);
						}
					}
				}
				else if (range.include(E_.pos)) {
					values[0] += E_.copy[0];
					values[1] += E_.copy[1];
				}
				else if (E_.pos.overlap(range)) {
					if (E_.pos.begin < range.begin) {
						if (range.begin - E_.pos.begin < E_.pos.end - range.begin) {
							srange r1(E_.pos.begin, range.begin - 1);
							values[0] += E_.copy[0] - depthSum(r, r1, cna->data);
							values[1] += E_.copy[1] - (cna->control ? depthSum(r, r1, cna->control) : r1.length(true));
						}
						else {
							srange r2(range.begin, E_.pos.end);
							values[0] += depthSum(r, r2, cna->data);
							values[1] += (cna->control ? depthSum(r, r2, cna->control) : r2.length(true));
						}
					}
					else {
						if (range.end - E_.pos.begin < E_.pos.end - range.end) {
							srange r1(E_.pos.begin, range.end);
							values[0] += depthSum(r, r1, cna->data);
							values[1] += (cna->control ? depthSum(r, r1, cna->control) : r1.length(true));
						}
						else {
							srange r2(range.end + 1, E_.pos.end);
							values[0] += E_.copy[0] - depthSum(r, r2, cna->data);
							values[1] += E_.copy[1] - (cna->control ? depthSum(r, r2, cna->control) : r2.length(true));
						}
					}
				}
			}
		}
	}
	
	
	
}
*/
#endif
