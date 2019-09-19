#ifndef __ANALYZER_H__
#define __ANALYZER_H__

#define READ_INFO_VIEW 0x0001
#define MAKE_SUMMARY 0x0002
#define COPY_NUM_VIEW 0x0004
#define SPLIT_READ_VIEW 0x0008
#define SUMMARY_INTEGRATE 0x0010
#define VAR_SEARCH 0x0020
#define MERGE_VLIST 0x0040
#define SUBTRACT_VLIST 0x0080
#define COMMON_VLIST 0x0100
#define TEMPLATE_EXPORT 0x1000

#define CNV_DETECTION 0x1000
#define SRV_DETECTION 0x2000
#define CNSRV_DETECTION 0x4000

#define INSUFFICIENT_ERROR 0x0001
#define FILE_NOT_FOUND_ERROR 0x0002
#define SEQENCE_SUMMARY_FILE_ERROR 0x0002
#define VARIANT_FILE_ERROR 0x0004
#define DATABASE_ERROR 0x0008
#define TARGET_FILE_ERROR 0x0010

#include "slib/sexception.hpp"
#include "slib/sutil.hpp"
#include "sio/sjson.hpp"
#include "sio/sxml.hpp"

#include "sbioinfo/sbioexception.hpp"
#include "sngs/sbam.hpp"
#include "sbioinfo/sbioseq.hpp"

#include "param.hpp"

using namespace slib;
using namespace slib::sio;
using namespace slib::smath;
using namespace slib::sbio;
using namespace slib::sngs;

namespace stk {
    struct status_t {
        bool analysis;
        int err;
        sint total_task;
        vecNd<sint> current_task;
        sstring src, msg;
        std::mutex mtx;
        
        status_t();
        ~status_t();
        
        void abort(sexception *ex);
        double progress();
        void init();
    };
    
    class analyzer;
    
    class bamreader;
    class realigner;
    class sranalysis;
    class cnanalysis;
    class varsearch;
    
    extern void mapping(sbamsummary *sum, salign_t *al, bool sync = false);
    extern void addSR(sbamsummary *sum, variant_t *var, bool sync = false);
    extern int srvtype(variant_t *var);
    
    class analyzer {
    private:
        bamreader *reader;
        sranalysis *sran;
        cnanalysis *cnan;
        varsearch *search;
        
    public:
        stk_param_t *par;
        slist<sbamfile *> bam;
        sbamsummary *summary;
        varlist *variants;
        stpool threads;
        status_t *status;
        
    public:
        analyzer();
        analyzer(stk_param_t *p);
        ~analyzer();
        
        void readview(const sfile &infile, sfile &outfile);
        void summarize(const sfile &infile, sfile &outfile);
        void sumintegrate(const slist<sfile> &infiles, sfile &outfile);
        void copynum(const sfile &infile, sfile &outfile);
        void splitread(const sfile &infile, sfile &outfile);
        void vsearch(int method, const sfile &file, sfile &outfile);
        void analyze(int method, const sfile &file, sfile &outfile);
        void vlmerge(const slist<sfile> &infiles, sfile &outfile);
        void subtract(const sfile &infile, sfile &outfile);
        void common(const slist<sfile> &infiles, sfile &outfile);
        
        void setParam(stk_param_t *p);
        void init();
    };
    
    #define clip_vec slist<uint64_t>
    
    class bamreader {
    private:
        stk_param_t *par;
        status_t *status;
        sdnafile *ref;
        slist<sbamfile *> *bam;
        sbamsummary *summary;
        stpool *threads;
        
    public:
        slist<clip_vec> clips;
        
    public:
        bamreader(analyzer *a);
        ~bamreader();
        
        void read(slist<sbreadinfo> *list);
        void summarize(bool paired);
        void init();
    };
    
    struct split_read_t {
        salign_t *align;
        bool head_clip, read_dir;
        sbioseq seq;
        
        split_read_t();
        split_read_t(const split_read_t &sr);
        ~split_read_t();
        
        split_read_t &operator=(const split_read_t &sr);
        void setHead(salign_t *a);
        void setTail(salign_t *a);
        inline void init();
    };
    
#define ri_vec std::vector<sbreadinfo>
#define sr_vec std::vector<split_read_t>
    
    class sranalysis {
    private:
        stk_param_t *par;
        sdnafile *ref;
        slist<sbamfile *> *bam;
        sbamsummary *summary;
        status_t *status;
        stpool *threads;
        std::vector<std::mutex> mtx;
        std::vector<ri_vec> rivs;
        std::vector<sr_vec> srvs;
        slist<sbquery *> ques;
        slist<slocalsearch *> searches;
        slist<slocalextend *> extends;
        
    private:
        void sranalyze(int r, clip_vec *cv);
        
    public:
        sranalysis(analyzer *a);
        ~sranalysis();
        
        void analyze();
        void analyze(slist<clip_vec> &cv);
        void init();
    };
    
    struct cn_var_t {
        int8_t type;
        int32_t pos, len;
        double sread, sread2, bread, bread2;
        
        cn_var_t();
        cn_var_t(const cn_var_t &cnv);
        ~cn_var_t();
        
        bool operator<(const cn_var_t &cnv) const;
        bool operator==(const cn_var_t &cnv) const;
    };
    
#define cn_vec slist<cn_var_t>
    
    class cnanalysis {
    private:
        variant_param_t *vpar;
        cnvar_param_t *par;
        sdnafile *ref;
        sbamsummary *summary, *background;
        float ratio;
        stpool *threads;
        std::mutex mtx;
        
    public:
        slist<cn_vec> variants;
        slist<smath::vecNd<int32_t>> cnvidx;
        
    private:
        void searchCNV(int r);
        
    public:
        cnanalysis(analyzer *a);
        ~cnanalysis();
        
        void search();
        void makelist(varlist *vl);
        float depthIn(int32_t ref, int32_t pos, int32_t len, bool bg = false);
        void init();
    };
    
    class varsearch {
    private:
        stk_param_t *par;
        variant_param_t *vpar;
        cnvar_param_t *cvpar;
        srvar_param_t *svpar;
        sdnafile *ref;
        sbamsummary *summary, *background;
        varlist *vlist;
        cnanalysis *cnan;
        stpool *threads;
        
        slist<slist<int>> sr_dels, sr_dups;
        slist<slist<std::pair<int, int>>> sr_indels, sr_invs;
        slist<slist<slist<std::pair<int, int>>>> sr_trs, sr_trinvs;
        
    private:
        void makeDelList(int r);
        void makeDelInfoList(int r, int v);
        void makeDelCPInfoList(int r, int v);
        void makeDupList(int r);
        void makeDupInfoList(int r, int v);
        void makeDupCPInfoList(int r, int v);
        void makeInDelList1(int r);
        void makeInDelList2(int r);
        void makeInDelInfoList(int r, int v, int i);
        void makeInDelCPInfoList(int r, int v, int i);
        void makeInvList(int r);
        void makeInvInfoList(int r, int v);
        void makeInvCPInfoList(int r, int v);
        void makeInvInsList1(int r);
        void makeInvInsList2(int r);
        void makeInvInsInfoList(int r, int v, int i);
        void makeInvInsCPInfoList(int r, int v, int i);
        void makeTrsList(int i, int j);
        void makeTrsInfoList(int i, int j, int v);
        void makeTrsCPInfoList(int i, int j, int v);
        void makeTrsInsList1(int i, int j);
        void makeTrsInsList2(int i, int j);
        void makeTrsInsInfoList(int i, int j, int v, int k);
        void makeTrsInsCPInfoList(int i, int j, int v, int k);
        void makeTrsInvList(int i, int j);
        void makeTrsInvInfoList(int i, int j, int v);
        void makeTrsInvCPInfoList(int i, int j, int v);
        void makeTrsInvInsList1(int i, int j);
        void makeTrsInvInsList2(int i, int j);
        void makeTrsInvInsInfoList(int i, int j, int v, int k);
        void makeTrsInvInsCPInfoList(int i, int j, int v, int k);
        
        void searchSRVar();
        void makeSRVInfo();
        void makeCNSRVInfo();
        
    public:
        varsearch(analyzer *a);
        ~varsearch();
        
        void search(int mode);
        void init();
    };
}

#endif
