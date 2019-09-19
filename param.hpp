#ifndef __PARAM_H__
#define __PARAM_H__

#include "slib/sutil.hpp"
#include "sbioinfo/sbioannot.hpp"
#include "sbioinfo/slsearch.hpp"
#include "sngs/sbam.hpp"

using namespace slib;
using namespace slib::sio;
using namespace slib::sbio;
using namespace slib::sngs;

#define DEFAULT_THREAD 16
#define DEFAULT_CP_BIN 1000
#define DEFAULT_CLIP_LENGTH 20
#define DEFAULT_CLIP_COUNT 200

namespace stk {
    struct analyze_param_t {
        int32_t thread_count, cp_bin, max_clip_count, min_clip_length;
        
        analyze_param_t();
        ~analyze_param_t();
        
        void set(sobj_ptr &obj);
        sobj_ptr toObj();
    };
    
    struct annot_param_t {
        bool var_filter[5];
        
        annot_param_t();
        ~annot_param_t();
        
        void set(sobj_ptr &obj);
        sobj_ptr toObj();
    };
    
    struct stk_param_t {
        sdnafile *ref;
        sbannotdb *db;
        sbamsummary *bg;
        varlist *ctrl;
        slist<sregion> target, ntarget;
        sstring outformat;
        
        analyze_param_t an_par;
        sbam_param_t bam_par;
        slsearch_param_t ls_par;
        variant_param_t var_par;
        annot_param_t ant_par;
        
        stk_param_t();
        ~stk_param_t();
        
        void set(sobj_ptr &obj);
        sobj_ptr toObj();
        
        void load(const char *path);
        
        void saveA(const char *path);
        void saveV(const char *path);
        void setRefPath(const char *path);
        void setDBPath(const char *path);
        void setBGPath(const char *path);
        void setCtrlPath(const char *path);
        void setTgtPath(const char *path);
        void setNtgtPath(const char *path);
    };
}




#endif
