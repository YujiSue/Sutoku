#include "analyzer.hpp"

using namespace stk;

const float default_bg = 1.0f;

cn_var_t::cn_var_t() {
    type = 0;
    pos = -1;
    len = 0;
    sread = 0.0;
    sread2 = 0.0;
    bread = 0.0;
    bread2 = 0.0;
}
cn_var_t::cn_var_t(const cn_var_t &cnv) {
    type = cnv.type; pos = cnv.pos; len = cnv.len;
    sread = cnv.sread; sread2 = cnv.sread2; bread = cnv.bread; bread2 = cnv.bread2;
}
cn_var_t::~cn_var_t() {}

bool cn_var_t::operator<(const cn_var_t &cnv) const {
    return pos < cnv.pos;
}
bool cn_var_t::operator==(const cn_var_t &cnv) const {
    return type == cnv.type && pos == cnv.pos && len == cnv.len;
}

cnanalysis::cnanalysis(analyzer *a) {
    vpar = &a->par->var_par;
    par = &a->par->var_par.cnv_par;
    ref = a->par->ref;
    summary = a->summary;
    background = a->par->bg;
    threads = &a->threads;
    variants.resize(ref->size());
    cnvidx.resize(ref->size());
    sforeachi(i, *ref) cnvidx[i].resize(((ref->length(i)-1)>>14)+1, -1);
}
cnanalysis::~cnanalysis() {}
/*
inline bool sequential(cn_var_t &v1, int i, cn_vec &variants, double &ratio, cnvar_param_t *par) {
    cn_var_t &v2 = variants[i], &v3 = variants[i+1];
    if (!v2.type && v1.type == v3.type) {
        if (v1.type == DELETION &&
            (v1.sread+v2.sread+v3.sread)/(v1.bread+v2.bread+v3.bread)*ratio <= par->del_cp)
            return true;
        else if(v1.type == DUPLICATION &&
                par->dup_cp <= (v1.sread+v2.sread+v3.sread)/(v1.bread+v2.bread+v3.bread)*ratio)
            return true;
    }
    return false;
}
 */
inline float cpyRatio(float *sdp, float *bdp, float &ratio) {
    if (bdp) {
        if (*bdp < SMATH_FEPS) return INFINITY;
        return (*sdp)/(*bdp)*ratio;
    }
    else return (*sdp)*ratio;
}
inline bool isLCVInit(float &cp, float &edge, cnvar_param_t *par) {
    return cp <= par->del_cp && edge <= -par->cp_diff;
}
inline bool isLCV(float &cp, cnvar_param_t *par) { return cp <= par->del_cp; }

inline bool isLCVEnd(float &cp, float &edge, cnvar_param_t *par) {
    return par->del_cp < cp || par->cp_diff < smath::abs(edge);
}
inline bool isHCVInit(float &cp, float &edge, cnvar_param_t *par) {
    return par->dup_cp <= cp && par->cp_diff <= edge;
}
inline bool isHCV(float &cp, cnvar_param_t *par) { return par->dup_cp <= cp; }
inline bool isHCVEnd(float &cp, float &edge, cnvar_param_t *par) {
    return cp < par->dup_cp || par->cp_diff <= edge;
}
inline void makeCV(int type, int &idx, float *sdp, float *bdp, int &bin, cn_var_t &cnv) {
    cnv.type = type;
    cnv.pos = idx*bin;
    cnv.len = bin;
    cnv.sread = (*sdp)*bin;
    cnv.sread2 = cnv.sread*cnv.sread;
    cnv.bread = (bdp?(*bdp):1.0f)*bin;
    cnv.bread2 = cnv.bread*cnv.bread;
}
inline void addCP(cn_var_t &cnv, float *sdp, float *bdp, int &bin) {
    cnv.len += bin;
    cnv.sread += *sdp;
    cnv.sread2 += (*sdp)*(*sdp);
    cnv.bread += bdp?*bdp:1.0f;
    cnv.bread2 += bdp?(*bdp)*(*bdp):1.0f;
}
inline void extendCV(cn_var_t &cnv, float *cp, float &edge, float *sdp, float *bdp, float &ratio,
                     int &bin, int &idx, int &size, cn_vec &vec, cnvar_param_t *par) {
    switch (cnv.type) {
        case SB_DELETION:
        {
            while (idx < size) {
                cp[1] = cpyRatio(sdp, bdp, ratio);
                edge = cp[1]-cp[0];
                if (isLCVEnd(cp[1], edge, par)) {
                    vec.add(cnv);
                    if (isLCV(cp[1], par)) {
                        makeCV(SB_DELETION, idx, sdp, bdp, bin, cnv);
                        extendCV(cnv, cp, edge, sdp, bdp, ratio, bin, idx, size, vec, par);
                    }
                    break;
                }
                else addCP(cnv, sdp, bdp, bin);
                ++idx; ++sdp; if(bdp) ++bdp; cp[0] = cp[1];
            }
            break;
        }
        case SB_DUPLICATION:
        {
            while (idx < size) {
                cp[1] = cpyRatio(sdp, bdp, ratio);
                edge = cp[1]-cp[0];
                if (isHCVEnd(cp[1], edge, par)) {
                    vec.add(cnv);
                    if (isHCV(cp[1], par)) {
                        
                        std::cout<<"HCV:"<<cp[0]<<"-"<<cp[1]<<std::endl;
                        
                        makeCV(SB_DUPLICATION, idx, sdp, bdp, bin, cnv);
                        extendCV(cnv, cp, edge, sdp, bdp, ratio, bin, idx, size, vec, par);
                    }
                    break;
                }
                else addCP(cnv, sdp, bdp, bin);
                ++idx; ++sdp; if(bdp) ++bdp; cp[0] = cp[1];
            }
            break;
        }
        default:
        {
            while (idx < size) {
                cp[1] = cpyRatio(sdp, bdp, ratio);
                edge = cp[1]-cp[0];
                if (isLCVInit(cp[1], edge, par) || isHCVInit(cp[1], edge, par)) {
                    
                    std::cout<<"NCV:"<<cp[0]<<"-"<<cp[1]<<std::endl;
                    
                    
                    vec.add(cnv);
                    break;
                }
                else addCP(cnv, sdp, bdp, bin);
                ++idx; ++sdp; if(bdp) ++bdp; cp[0] = cp[1];
            }
            break;
        }
    }
}

void cnanalysis::searchCNV(int r) {
    if (!summary->depth_size[r]) return;
    int size = summary->depth_size[r], idx = 0, bin = summary->bin;
    auto sdp = summary->depth[r].data(),
    bdp = background?background->depth[r].data():nullptr;
    float cp[2], edge;
    cn_var_t cnv;
    cp[0] = 1.0;
    cp[1] = cpyRatio(sdp, bdp, ratio);
    edge = cp[1]-cp[0];
    while (idx < size) {
        if (isLCVInit(cp[1], edge, par))
            makeCV(SB_DELETION, idx, sdp, bdp, bin, cnv);
        else if (isHCVInit(cp[1], edge, par))
            makeCV(SB_DUPLICATION, idx, sdp, bdp, bin, cnv);
        else makeCV(0, idx, sdp, bdp, bin, cnv);
        ++idx; ++sdp; if(bdp) ++bdp; cp[0] = cp[1];
        extendCV(cnv, cp, edge, sdp, bdp, ratio, bin, idx, size, variants[r], par);
    }
}

void cnanalysis::search() {
    if (SMATH_FEPS < summary->average_depth)
        ratio = (background?background->average_depth:default_bg)/summary->average_depth;
    else return;
    
    
    sforeachi(i, *ref) searchCNV(i);
    /*
    sforeachi(i, *ref) threads->addTask(&cnanalysis::searchCNV, this, i);
    threads->complete();
     */
}

void cnanalysis::makelist(varlist *vl) {
    sforeachi(r, *ref) {
        cn_vec &cnvs = variants[r];
        for (int i = 0; i < cnvs.size()-2; ++i) {
            if (!cnvs[i].type) continue;
            else {
                cn_var_t v = cnvs[i];
                /*
                if(sequential(v, i+1, cnvs, ratio, par)) {
                    while (i < cnvs.size()-2 && sequential(v, i+1, cnvs, ratio, par)) {
                        v += cnvs[++i];
                        v += cnvs[++i];
                    }
                }
                 */
                if ((v.type == SB_DELETION && vpar->min_length[0] <= v.len) ||
                    (v.type == SB_DUPLICATION && vpar->min_length[1] <= v.len)) {
                    varinfo_t vi;
                    vi.var_method = COPY_NUM_VARIANT;
                    vi.variant.type = v.type;
                    vi.variant.pos1.idx = r;
                    vi.variant.pos1.pos = v.pos*summary->bin+1;
                    if(v.pos+v.len == summary->depth_size[r])
                        vi.variant.pos1.len = ref->length(r)-vi.variant.pos1.pos+1;
                    else
                        vi.variant.pos1.len = v.len*summary->bin;
                    vi.cp.copy[0] = (double)(v.sread*summary->bin)/vi.variant.pos1.len;
                    vi.cp.bgcopy[0] = background?(double)(v.bread*summary->bin)/vi.variant.pos1.len:default_bg;
                    vi.cp.copy[1] = vi.cp.copy[0]/summary->average_depth;
                    vi.cp.bgcopy[1] = (background&&0.0<background->average_depth)?
                    vi.cp.bgcopy[0]/background->average_depth:default_bg;
                    vi.cp.copy[2] = vi.cp.copy[1]/vi.cp.bgcopy[1];
                    if (vi.variant.type == SB_DELETION &&
                        vi.cp.copy[2] <= par->homo_del_cp) vi.cp.homo = true;
                    else if (vi.variant.type == SB_INSERTION &&
                             par->mul_cp <= vi.cp.copy[2]) vi.variant.type = SB_MULTIPLICATION;
                    vl->add(vi);
                }
            }
        }
    }
}

float cnanalysis::depthIn(int32_t ref, int32_t pos, int32_t len, bool bg) {
    srange vrange(pos, pos+len-1), cnvrange(-1, -1);
    double val = 0.0f;
    auto it = variants[ref].begin()+cnvidx[ref][pos>>14];
    while (it < variants[ref].end()) {
        if (vrange.end < it->pos) break;
        if (vrange.overlap(srange(it->pos, it->pos+it->len-1))) {
            if (cnvrange.begin < 0) cnvrange.begin = it->pos;
            cnvrange.end = it->pos+it->len-1;
            if (bg) val += it->bread;
            else val += it->sread;
        }
        ++it;
    }
    if (cnvrange.begin < vrange.begin)
        val -= summary->totalDepthIn(ref, cnvrange.begin, cnvrange.begin-cnvrange.begin);
    else if (vrange.begin < cnvrange.begin)
        val += summary->totalDepthIn(ref, pos, cnvrange.begin-cnvrange.begin);
    if (cnvrange.end < vrange.end)
        val += summary->totalDepthIn(ref, cnvrange.end, vrange.end-cnvrange.end);
    else if (vrange.end < cnvrange.end)
        val -= summary->totalDepthIn(ref, vrange.end, cnvrange.end-vrange.end);
    return val/len;
}
void cnanalysis::init() {
    sforeachi(i, *ref) {
        variants[i].clear();
        cnvidx[i].reset(-1);
    }
}
