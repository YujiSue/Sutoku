#ifndef STK_PARAM_H
#define STK_PARAM_H
#include "sobj.h"
#include "sbioinfo.h"
#include "sapp/sapp.h"
#include "util.h"
using namespace slib;
using namespace slib::sbio;
namespace stk {
    class Param {
    protected:
        String _logpath;
    public:
        // Target files
        stringarray inputs;
        String command, inputdir, outdir, output, oformat;

        // Reference
        SeqList reference;
        AnnotDB annotdb;

        // Target region
        Array<sregion> target;

        // Control
        NGSData control;
        VarList vcontrol;

        // Parameter set from slib
        SeqSearchParam seqp;
        VarParam varp;

        // Unique parameters
        // Sequence type : SINGLE / PAIRED
        sngs::SEQ_TYPE seqtype;
        // Read BAM in async. mode
        bool async_load;
        // Detect SV 
        bool detect_sv;
        // Ignore reads with duplicated flag (marked by Picard)
        bool ignore_dp;
        // Select only read with clipped region
        bool clipped;

        // Bin size for depth count
        int depth_bin;
        // Min. clip size to analyze
        size_t min_clip;
        // Buffer size to store clip information
        size_t clip_buffer;
        // Minimum error probability (<=> Maximum phred score)
        double min_err_prob;

        // Variant types to detect
        bool detect[10]; // DEL,DUP,INS,COMPLEX,INV,INV+INS,TRS,TRS+INS,TRS+INV,TRS+INV+INS

        // Use machine learning model to check copy number
        bool mlcheck;

        // Gene annotaion flag
        bool annotation;
        sushort var_site;

        // Multi thread proc.
        size_t max_thread;
        SWork threads;

        // Logger
        String logpath;
        SLogger logger;
        // Status
        Status status;        
        // Locker for sync. proc.
        SLock slock;
        Array<SLock> mlock;

    public:
        Param();
        Param(SDictionary& pref);
        ~Param();

        void load(const char* path);
        void save(const char *path);
        void set(SDictionary& pref);
        sobj toObj();
    };
}
#endif