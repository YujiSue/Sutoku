#include "analyzer.hpp"

using namespace stk;

const int MAX_DIRECT_COUNT = 1<<14;

varsearch::varsearch(analyzer *a) {
    par = a->par;
    vpar = &par->var_par;
    cvpar = &vpar->cnv_par;
    svpar = &vpar->srv_par;
    ref = a->par->ref;
    summary = a->summary;
    background = a->par->bg;
    threads = &a->threads;
    
    sr_dels.resize(ref->size());
    sr_dups.resize(ref->size());
    sr_indels.resize(2*ref->size());
    
    sr_invs.resize(3*ref->size());
    sr_trs.resize(ref->size());
    sr_trinvs.resize(ref->size());
    for (int i = 0; i < ref->size(); ++i) {
        sr_trs[i].resize(3*ref->size());
        sr_trinvs[i].resize(3*ref->size());
    }
    cnan = new cnanalysis(a);
    vlist = a->variants;
}
varsearch::~varsearch() {
    if (cnan) delete cnan;
}

inline void combinate(varinfo_t *vi, variant_t *v) {
    vi->variant += *v;
    vi->variant.alt += "/"+v->alt;
}
inline void oneStart(varinfo_t *vi) {
    vi->variant.pos1.pos++;
    if (-1 < vi->variant.pos2.idx) vi->variant.pos2.pos++;
}

inline bool srVarCheck1(int i, variant_t &var, variant_param_t *vp) {
    return vp->srv_par.min_sr[i] <= var.total() &&
    readBias(var.pread, var.nread) <= vp->srv_par.max_fr_bias;/* &&
    vp->min_qual <= phredVal(var.qual);*/
}
inline bool srVarCheck2(int i, variant_t &v1, variant_t &v2, variant_param_t *vp) {
    return vp->srv_par.min_comp_sr[i] <= v1.total()+v2.total() &&
    readBias(v1.total(), v2.total()) <= vp->srv_par.max_comp_bias;/* &&
    vp->min_qual <= phredVal(1.0-(1.0-v1.qual)*(1.0-v2.qual));*/
}
inline void calcCP(sbpos_t *pos, double *copy, double *bgcopy, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan) {
    /*
    if (pos->len < MAX_DIRECT_COUNT) {
        copy[0] = cnan->depthIn(pos->idx, pos->pos-1, pos->len);
        if (bg) bgcopy[0] = cnan->depthIn(pos->idx,  pos->pos-1,  pos->len, true);
        else bgcopy[0] = 1.0f;
    }
    else {
     */
        copy[0] = sum->depthIn(pos->idx, pos->pos-1, pos->len);
        if (bg) bgcopy[0] = bg->depthIn(pos->idx, pos->pos-1, pos->len);
        else bgcopy[0] = 1.0f;
    //}
    if (0.0 < sum->average_depth) copy[1]=copy[0]/sum->average_depth;
    if (bg) {
        if (0.0 < bg->average_depth) bgcopy[1]=bgcopy[0]/bg->average_depth;
        else bgcopy[1] = bgcopy[0];
    }
    else bgcopy[1] = 1.0f;
    if (0.0 < bgcopy[1]) copy[2] = copy[1]/bgcopy[1];
    else copy[2] = INFINITY;
}
inline void calcCP2(sbpos_t *pos1, sbpos_t *pos2, double *copy, double *bgcopy, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan) {
    /*
    if (pos1->len < MAX_DIRECT_COUNT) {
        copy[0] = cnan->depthIn(pos1->idx, pos1->pos-1, pos1->len);
        if (bg) bgcopy[0] = cnan->depthIn(pos1->idx,  pos1->pos-1,  pos1->len, true);
        else bgcopy[0] = 1.0f;
    }
    else {
     */
        copy[0] = sum->depthIn(pos1->idx, pos1->pos-1, pos1->len);
        if (bg) bgcopy[0] = bg->depthIn(pos1->idx, pos1->pos-1, pos1->len);
        else bgcopy[0] = 1.0f;
    //}
    /*
    if (pos2->len < MAX_DIRECT_COUNT) {
        copy[0] += cnan->depthIn(pos2->idx, pos2->pos-1, pos2->len);
        if (bg) bgcopy[0] += cnan->depthIn(pos2->idx,  pos2->pos-1,  pos2->len, true);
        else bgcopy[0] += 1.0f;
    }
    else {
     */
        copy[0] += sum->depthIn(pos2->idx, pos2->pos-1, pos2->len);
        if (bg) bgcopy[0] += bg->depthIn(pos2->idx, pos2->pos-1, pos2->len);
        else bgcopy[0] += 1.0f;
    //}
    copy[0]/=2.0; bgcopy[0]/=2.0;
    if (0.0 < sum->average_depth) copy[1]=copy[0]/sum->average_depth;
    if (bg) {
        if (0.0 < bg->average_depth) bgcopy[1] = bgcopy[0]/bg->average_depth;
        else bgcopy[1] = bgcopy[0];
    }
    else bgcopy[1] = 1.0f;
    if (0.0 < bgcopy[1]) copy[2] = copy[1]/bgcopy[1];
    else copy[2] = INFINITY;
}

void varsearch::makeDelList(int r) {
    auto &variants = summary->variants[5*r];
    auto &list = sr_dels[r];
    if (!vpar->srv_par.detect_var[0] || variants.empty()) return;
    sforeach(variants) {
        if(srVarCheck1(0, selement, vpar) && vpar->min_length[0] <= selement.pos1.len) list.add(it-variants.begin());
    }
}
inline void makeDelInfo(variant_t *var, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var;
    vi->variant.pos1.pos++;
    vi->variant.pos2.init();
    vi->variant.pos2.len = vi->variant.alt.length();
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    oneStart(vi);
}
inline void calcDelCp(varinfo_t *vi, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan, variant_param_t *par) {
    calcCP(&vi->variant.pos1, vi->cp.copy, vi->cp.bgcopy, sum, bg, cnan);
    if (vi->cp.copy[2] <= par->cnv_par.homo_del_cp) vi->cp.homo = true;
    if (par->cnv_par.del_cp < vi->cp.copy[2] || (bg && vi->cp.bgcopy[0] < par->cnv_par.bg_depth)) vi->variant.type = 0;
    else vi->var_method |= COPY_NUM_VARIANT;
}
void varsearch::makeDelInfoList(int r, int v) {
    auto &variants = summary->variants[5*r];
    auto &list = sr_dels[r];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeachi(i, list) {
        makeDelInfo(&variants[list[i]], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeDelCPInfoList(int r, int v) {
    auto &variants = summary->variants[5*r];
    auto &list = sr_dels[r];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeachi(i, list)  {
        makeDelInfo(&variants[list[i]], vptr, summary, ref);
        calcDelCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}
void varsearch::makeDupList(int r) {
    auto &variants = summary->variants[5*r+1];
    auto &list = sr_dups[r];
    if (!vpar->srv_par.detect_var[1] || variants.empty()) return;
    sforeach(variants) {
        if(srVarCheck1(1, selement, vpar) && vpar->min_length[1] <= selement.pos1.len) list.add(it-variants.begin());
    }
}
inline void makeDupInfo(variant_t *var, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var;
    vi->variant.pos1.pos = var->pos2.pos;
    vi->variant.pos2.init();
    vi->variant.pos2.len = vi->variant.alt.length();
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][vi->variant.pos1.pos/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len-1)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    oneStart(vi);
}
inline void calcDupCp(varinfo_t *vi, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan, variant_param_t *par) {
    calcCP(&vi->variant.pos1, vi->cp.copy, vi->cp.bgcopy, sum, bg, cnan);
    if (par->cnv_par.homo_freq <= vi->cp.frequency) vi->cp.homo = true;
    if (vi->cp.copy[2] < par->cnv_par.dup_cp || (bg && vi->cp.bgcopy[0] < par->cnv_par.bg_depth)) vi->variant.type = 0;
    else {
        if (par->cnv_par.mul_cp <= vi->cp.copy[2]) vi->variant.type = SB_MULTIPLICATION;
        vi->var_method |= COPY_NUM_VARIANT;
    }
}
void varsearch::makeDupInfoList(int r, int v) {
    auto &variants = summary->variants[5*r+1];
    auto &list = sr_dups[r];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeachi(i, list) {
        makeDupInfo(&variants[list[i]], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeDupCPInfoList(int r, int v) {
    auto &variants = summary->variants[5*r+1];
    auto &list = sr_dups[r];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeachi(i, list) {
        makeDupInfo(&variants[list[i]], vptr, summary, ref);
        calcDupCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//InDel////////////
//  1 2     2'1' //
// -> ->   -> -> //
//  2 1     1'2' //
// -> ->   -> -> //
///////////////////
void varsearch::makeInDelList1(int r) {
    auto &del = summary->variants[5*r], &ins = summary->variants[5*r+1];
    auto &list = sr_indels[2*r];
    auto &vidx = summary->insidx[r];
    if (!vpar->srv_par.detect_comp_var[0] || del.empty() || ins.empty()) return;
    for (auto dit = del.begin(); dit < del.end(); ++dit) {
        auto start = ins.begin()+vidx[dit->pos2.pos>>14];
        for (auto iit = start; iit < ins.end(); ++iit) {
            if (srVarCheck1(0, *dit, vpar) &&
                srVarCheck1(1, *iit, vpar) &&
                srVarCheck2(0, *dit, *iit, vpar) &&
                iit->pos2.pos < dit->pos2.pos &&
                dit->pos2.pos < iit->pos1.pos &&
                dit->pos1.pos-vpar->max_dist <= iit->pos2.pos &&
                (vpar->min_length[0] <= iit->pos2.pos-dit->pos1.pos ||
                 vpar->min_length[1] <= iit->pos1.pos-dit->pos2.pos+1))
                list.add(std::make_pair(dit-del.begin(), iit-ins.begin()));
        }
    }
}
void varsearch::makeInDelList2(int r) {
    auto &del = summary->variants[5*r], &ins = summary->variants[5*r+1];
    auto &list = sr_indels[2*r+1];
    auto &vidx = summary->delidx[r];
    if (!vpar->srv_par.detect_comp_var[0] || del.empty() || ins.empty()) return;
    for (auto iit = ins.begin(); iit < ins.end(); ++iit) {
        auto start = del.begin()+vidx[iit->pos2.pos>>14];
        for (auto dit = start; dit < del.end(); ++dit) {
            if (srVarCheck1(0, *dit, vpar) &&
                srVarCheck1(1, *iit, vpar) &&
                srVarCheck2(0, *dit, *iit, vpar) &&
                iit->pos2.pos < dit->pos1.pos &&
                dit->pos1.pos < iit->pos1.pos &&
                iit->pos1.pos-vpar->max_dist <= dit->pos2.pos &&
                (vpar->min_length[0] <= dit->pos2.pos-iit->pos1.pos ||
                 vpar->min_length[1] <= dit->pos1.pos-iit->pos2.pos+1))
                list.add(std::make_pair(iit-ins.begin(), dit-del.begin()));
        }
    }
}
inline void makeInDelInfo(variant_t *var1, variant_t *var2, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var1;
    vi->variant.type = SB_INSERTION;
    vi->variant.pos1.len = var2->pos2.pos-var1->pos1.pos-1;
    if(0 < vi->variant.pos1.len) {
        vi->variant.pos1.pos++;
        vi->variant.type |= SB_DELETION;
    }
    else {
        vi->variant.pos1.pos = var2->pos2.pos;
        vi->variant.pos1.len = 1;
    }
    vi->variant.pos2.len = var2->pos1.pos-var1->pos2.pos+1;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, var2);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][vi->variant.pos2.pos/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len-1)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
}
inline void calcInDelCp(varinfo_t *vi, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan, variant_param_t *par) {
    calcCP(&vi->variant.pos2, &vi->cp.copy[3], &vi->cp.bgcopy[2], sum, bg, cnan);
    if(vi->variant.type&SB_DELETION) {
        calcCP(&vi->variant.pos1, vi->cp.copy, vi->cp.bgcopy, sum, bg, cnan);
        if (vi->cp.copy[2] <= par->cnv_par.homo_del_cp) vi->cp.homo = true;
        if (par->cnv_par.del_cp < vi->cp.copy[2] ||
            (bg && vi->cp.bgcopy[0] < par->cnv_par.bg_depth)) vi->variant.type = 0;
    }
    else if(par->cnv_par.homo_freq <= vi->cp.frequency) vi->cp.homo = true;
    if (vi->variant.type) {
        if (par->cnv_par.dup_cp <= vi->cp.copy[5]) vi->variant.type |= SB_DUPLICATION;
        if (par->cnv_par.mul_cp <= vi->cp.copy[5]) vi->variant.type |= SB_MULTIPLICATION;
        vi->var_method |= COPY_NUM_VARIANT;
    }
}
void varsearch::makeInDelInfoList(int r, int v, int i) {
    auto &del = summary->variants[5*r], &ins = summary->variants[5*r+1];
    auto &list = sr_indels[2*r+i];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeInDelInfo((i?(&ins[selement.first]):(&del[selement.first])),
                      (i?(&del[selement.second]):(&ins[selement.second])), vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeInDelCPInfoList(int r, int v, int i) {
    auto &del = summary->variants[5*r], &ins = summary->variants[5*r+1];
    auto &list = sr_indels[2*r+i];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeInDelInfo((i?(&ins[selement.first]):(&del[selement.first])),
                      (i?(&del[selement.second]):(&ins[selement.second])), vptr, summary, ref);
        calcInDelCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//Inv////////////////
//  1 1''   1'1''' //
// -> <-   <- ->   //
/////////////////////
void varsearch::makeInvList(int r) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r];
    auto &vidx = summary->invidx[2*r+1];
    if (!vpar->srv_par.detect_comp_var[1] || inv.empty()) return;
    auto beg = inv.begin()+vidx[0];
    if (beg->pos2.dir) return;
    for (auto iit = inv.begin(); iit < inv.end(); ++iit) {
        if (iit->pos1.dir) break;
        for (auto iit_ = beg; iit_ < inv.end(); ++iit_) {
            if (iit->pos2.pos - par->an_par.min_clip_length < iit_->pos1.pos) break;
            if (srVarCheck1(2, *iit, vpar) &&
                srVarCheck1(2, *iit_, vpar) &&
                srVarCheck2(1, *iit, *iit_, vpar) &&
                iit->pos1.pos-vpar->max_dist <= iit_->pos1.pos &&
                iit->pos2.pos-vpar->max_dist <= iit_->pos2.pos &&
                (vpar->min_length[2] <= iit->pos2.pos-iit_->pos1.pos+1))
                list.add(std::make_pair(iit-inv.begin(), iit_-inv.begin()));
        }
    }
}
inline void makeInvInfo(variant_t *var1, variant_t *var2,
                        varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var1;
    vi->variant.pos2 = var2->pos1;
    vi->variant.pos2.len = var1->pos2.pos-var2->pos1.pos+1;
    if (var1->pos1.pos < var2->pos1.pos) vi->variant.pos1.pos++;
    else vi->variant.pos1.pos = var2->pos1.pos;
    if (var1->pos2.pos < var2->pos2.pos) vi->variant.pos1.len = var2->pos2.pos-vi->variant.pos1.pos;
    else vi->variant.pos1.len = var1->pos2.pos-vi->variant.pos1.pos+1;
    if(vi->variant.pos2.len < vi->variant.pos1.len)
        vi->variant.type |= SB_DELETION;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, var2);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][vi->variant.pos2.pos/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len-1)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
}
inline void calcInvCp(varinfo_t *vi, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan, variant_param_t *par) {
    calcCP(&vi->variant.pos2, &vi->cp.copy[3], &vi->cp.bgcopy[2], sum, bg, cnan);
    if(vi->variant.type&SB_DELETION) {
        sbpos_t p1 = vi->variant.pos1, p2 = vi->variant.pos2;
        p1.len = vi->variant.pos2.pos-vi->variant.pos1.pos;
        p2.pos += p2.len;
        p2.len = p1.pos+p1.len-p2.pos-p2.len;
        calcCP2(&p1, &p2, vi->cp.copy, vi->cp.bgcopy, sum, bg, cnan);
        if (vi->cp.copy[2] <= par->cnv_par.homo_del_cp) vi->cp.homo = true;
        if (par->cnv_par.del_cp < vi->cp.copy[2] ||
            (bg && vi->cp.bgcopy[0] < par->cnv_par.bg_depth)) vi->variant.type = 0;
    }
    else if(par->cnv_par.homo_freq <= vi->cp.frequency) vi->cp.homo = true;
    if (vi->variant.type) {
        if (par->cnv_par.dup_cp <= vi->cp.copy[5]) vi->variant.type |= SB_DUPLICATION;
        if (par->cnv_par.mul_cp <= vi->cp.copy[5]) vi->variant.type |= SB_MULTIPLICATION;
        vi->var_method |= COPY_NUM_VARIANT;
    }
}
void varsearch::makeInvInfoList(int r, int v) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeInvInfo(&inv[selement.first], &inv[selement.second], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeInvCPInfoList(int r, int v) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeInvInfo(&inv[selement.first], &inv[selement.second], vptr, summary, ref);
        calcInvCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//InvIns/////////////////////
//  1 2'    2 1'      1'2  //
// -> <-   <- ->  =  <- -> //
//  1'2     2 1'      1 2' //
// -> <- = -> <-     <- -> //
/////////////////////////////
void varsearch::makeInvInsList1(int r) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r+1];
    auto &vidx = summary->invidx[2*r+1];
    if (!vpar->srv_par.detect_comp_var[0] || inv.empty()) return;
    auto beg = inv.begin()+vidx[0];
    if (beg->pos2.dir) return;
    for (auto iit = inv.begin(); iit < inv.end(); ++iit) {
        if (iit->pos1.dir) break;
        beg = inv.begin()+vidx[(iit->pos1.pos-vpar->max_dist)>>14];
        for (auto iit_ = beg; iit_ < inv.end(); ++iit_) {
            if (srVarCheck1(2, *iit, vpar) &&
                srVarCheck1(2, *iit_, vpar) &&
                srVarCheck2(0, *iit, *iit_, vpar) &&
                iit_->pos2.pos < iit->pos2.pos &&
                iit->pos1.pos-vpar->max_dist <= iit_->pos1.pos &&
                (vpar->min_length[0] <= iit_->pos1.pos-iit->pos1.pos ||
                 vpar->min_length[1] <= iit->pos2.pos-iit_->pos2.pos+1))
                list.add(std::make_pair(iit-inv.begin(), iit_-inv.begin()));
        }
    }
}
void varsearch::makeInvInsList2(int r) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r+1];
    auto &vidx = summary->invidx[2*r+1];
    if (!vpar->srv_par.detect_comp_var[0] || inv.empty()) return;
    auto beg = inv.begin()+vidx[0];
    if (beg->pos2.dir) return;
    for (auto iit = inv.begin(); iit < inv.end(); ++iit) {
        if (iit->pos1.dir) break;
        for (auto iit_ = beg; iit_ < inv.end(); ++iit_) {
            if (iit->pos1.pos < iit_->pos1.pos) break;
            if (srVarCheck1(2, *iit, vpar) &&
                srVarCheck1(2, *iit_, vpar) &&
                srVarCheck2(0, *iit, *iit_, vpar) &&
                iit->pos2.pos-vpar->max_dist <= iit_->pos2.pos &&
                (vpar->min_length[0] <= iit_->pos2.pos-iit->pos2.pos ||
                 vpar->min_length[1] <= iit->pos1.pos-iit_->pos1.pos+1))
                list.add(std::make_pair(iit_-inv.begin(), iit-inv.begin()));
        }
    }
}
inline void makeInvInsInfo(variant_t *var1, variant_t *var2, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    variant_t var_;
    if (var1->pos1.dir) {
        var_ = *var1;
        vi->variant = *var2;
        vi->variant.comp();
    }
    else {
        vi->variant = *var1;
        var_ = *var2;
        var_.comp();
    }
    vi->variant.type |= SB_INSERTION;
    vi->variant.pos1.len = var_.pos2.pos-vi->variant.pos1.pos-1;
    if (0 < vi->variant.pos1.len) {
        vi->variant.pos1.pos++;
        vi->variant.type |= SB_DELETION;
    }
    else {
        vi->variant.pos1.pos = var_.pos2.pos;
        vi->variant.pos1.len = 1;
    }
    vi->variant.pos2.len = vi->variant.pos2.pos-var_.pos1.pos+1;
    vi->variant.pos2.pos = var_.pos1.pos;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, &var_);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][vi->variant.pos2.pos/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len-1)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
}
void varsearch::makeInvInsInfoList(int r, int v, int i) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r+i];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeInvInfo(&inv[selement.first], &inv[selement.second], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeInvInsCPInfoList(int r, int v, int i) {
    auto &inv = summary->variants[5*r+2];
    auto &list = sr_invs[3*r+i];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeInvInfo(&inv[selement.first], &inv[selement.second], vptr, summary, ref);
        calcInDelCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//Trans/////
//  1 2'  //
// -> ->  //
//   x    //
// -> ->  //
//  2 1'  //
////////////
void varsearch::makeTrsList(int i, int j) {
    auto &low = summary->variants[5*i+3], &high = summary->variants[5*j+3];
    auto &list = sr_trs[i][3*j];
    auto &lidx = summary->trsidx[i][j], &hidx = summary->trsidx[j][i];
    if (!vpar->srv_par.detect_comp_var[2] || low.empty() || high.empty()) return;
    auto lbeg = low.begin()+lidx[0];
    auto hbeg = high.begin()+hidx[0];
    if (lbeg->pos2.idx != j || hbeg->pos2.idx != i) return;
    for (auto lit = lbeg; lit < low.end(); ++lit) {
        if (j < lit->pos2.idx) break;
        auto limit = lit->pos2.pos+vpar->max_dist;
        for (auto hit = hbeg; hit < high.end(); ++hit) {
            if (i < hit->pos2.idx || limit < hit->pos1.pos) break;
            if (srVarCheck1(3, *lit, vpar) &&
                srVarCheck1(3, *hit, vpar) &&
                srVarCheck2(2, *lit, *hit, vpar) &&
                lit->pos1.pos-vpar->max_dist <= hit->pos2.pos)
                list.add(std::make_pair(lit-low.begin(), hit-high.begin()));
        }
    }
}
inline void makeTrsInfo(variant_t *var1, variant_t *var2, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var1;
    vi->variant.pos1.len = var2->pos2.pos-var1->pos1.pos-1;
    if (0 < vi->variant.pos1.len) vi->variant.pos1.pos++;
    else {
        vi->variant.pos1.pos = var2->pos2.pos;
        vi->variant.pos1.len = 0;
    }
    vi->variant.pos2 = var2->pos1;
    vi->variant.pos2.len = var1->pos2.pos-var2->pos1.pos-1;
    if (0 < vi->variant.pos2.len) vi->variant.pos2.pos++;
    else {
        vi->variant.pos2.pos = var1->pos2.pos;
        vi->variant.pos2.len = 0;
    }
    if (vi->variant.pos1.len || vi->variant.pos2.len)
        vi->variant.type |= SB_DELETION;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, var2);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
}
inline void calcTrsCp(varinfo_t *vi, sbamsummary *sum, sbamsummary *bg, cnanalysis *cnan, variant_param_t *par) {
    if(vi->variant.type&SB_DELETION) {
        if (vi->variant.pos1.len) {
            calcCP(&vi->variant.pos1, vi->cp.copy, vi->cp.bgcopy, sum, bg, cnan);
            if (vi->cp.copy[2] <= par->cnv_par.homo_del_cp) vi->cp.homo = true;
            if (par->cnv_par.del_cp < vi->cp.copy[2] ||
                (bg && vi->cp.bgcopy[0] < par->cnv_par.bg_depth)) vi->variant.type = 0;
        }
        if (vi->variant.type && vi->variant.pos2.len) {
            calcCP(&vi->variant.pos2, &vi->cp.copy[3], &vi->cp.bgcopy[2], sum, bg, cnan);
            if (vi->cp.copy[5] <= par->cnv_par.homo_del_cp) vi->cp.homo = true;
            if (par->cnv_par.del_cp < vi->cp.copy[5] ||
                (bg && vi->cp.bgcopy[2] < par->cnv_par.bg_depth)) vi->variant.type = 0;
        }
    }
    else if(par->cnv_par.homo_freq <= vi->cp.frequency) vi->cp.homo = true;
    if (vi->variant.type) vi->var_method |= COPY_NUM_VARIANT;
}
void varsearch::makeTrsInfoList(int i, int j, int v) {
    auto &low = summary->variants[5*i+3], &high = summary->variants[5*j+3];
    auto &list = sr_trs[i][3*j];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInfo(&low[selement.first], &high[selement.second], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeTrsCPInfoList(int i, int j, int v) {
    auto &low = summary->variants[5*i+3], &high = summary->variants[5*j+3];
    auto &list = sr_trs[i][3*j];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInfo(&low[selement.first], &high[selement.second], vptr, summary, ref);
        calcTrsCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//TrsIns///////////
//  1 2     2'1' //
// -> ->   -> -> //
//  2 1     1'2' //
// -> ->   -> -> //
///////////////////
void varsearch::makeTrsInsList1(int i, int j) {
    auto &low = summary->variants[5*i+3], &high = summary->variants[5*j+3];
    auto &list = sr_trs[i][3*j+1];
    auto &lidx = summary->trsidx[i][j], &hidx = summary->trsidx[j][i];
    if (!vpar->srv_par.detect_comp_var[0] || low.empty() || high.empty()) return;
    auto lbeg = low.begin()+lidx[0];
    auto hbeg = high.begin()+hidx[0];
    if (lbeg->pos2.idx != j || hbeg->pos2.idx != i) return;
    for (auto lit = lbeg; lit < low.end(); ++lit) {
        if (j < lit->pos2.idx) break;
        hbeg = high.begin()+hidx[lit->pos2.pos>>14];
        for (auto hit = hbeg; hit < high.end(); ++hit) {
            if (i < hit->pos2.idx) break;
            if (srVarCheck1(3, *lit, vpar) &&
                srVarCheck1(3, *hit, vpar) &&
                srVarCheck2(0, *lit, *hit, vpar) &&
                lit->pos2.pos < hit->pos1.pos &&
                lit->pos1.pos-vpar->max_dist <= hit->pos2.pos &&
                (vpar->min_length[0] <= hit->pos2.pos-lit->pos1.pos ||
                 vpar->min_length[1] <= hit->pos1.pos-lit->pos2.pos))
                list.add(std::make_pair(lit-low.begin(), hit-high.begin()));
        }
    }
}
void varsearch::makeTrsInsList2(int i, int j) {
    auto &low = summary->variants[5*j+3], &high = summary->variants[5*i+3];
    auto &list = sr_trs[i][3*j+1];
    auto &lidx = summary->trsidx[j][i], &hidx = summary->trsidx[i][j];
    if (!vpar->srv_par.detect_comp_var[0] || low.empty() || high.empty()) return;
    auto lbeg = low.begin()+lidx[0];
    auto hbeg = high.begin()+hidx[0];
    if (lbeg->pos2.idx != i || hbeg->pos2.idx != j) return;
    for (auto hit = hbeg; hit < high.end(); ++hit) {
        if (j < hit->pos2.idx) break;
        lbeg = low.begin()+lidx[hit->pos2.pos>>14];
        for (auto lit = lbeg; lit < low.end(); ++lit) {
            if (j < lit->pos2.idx) break;
            if (srVarCheck1(3, *lit, vpar) &&
                srVarCheck1(3, *hit, vpar) &&
                srVarCheck2(0, *lit, *hit, vpar) &&
                hit->pos2.pos < lit->pos1.pos &&
                hit->pos1.pos-vpar->max_dist <= lit->pos2.pos &&
                (vpar->min_length[0] <= lit->pos2.pos-hit->pos1.pos ||
                 vpar->min_length[1] <= lit->pos1.pos-hit->pos2.pos))
                list.add(std::make_pair(hit-high.begin(), lit-low.begin()));
        }
    }
}
inline void makeTrsInsInfo(variant_t *var1, variant_t *var2, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var1;
    vi->variant.type |= SB_INSERTION;
    vi->variant.pos1.len = var2->pos2.pos-var1->pos1.pos-1;
    if (0 < vi->variant.pos1.len) {
        vi->variant.pos1.pos++;
        vi->variant.type |= SB_DELETION;
    }
    else {
        vi->variant.pos1.pos = var2->pos2.pos;
        vi->variant.pos1.len = 1;
    }
    vi->variant.pos2.len = var2->pos1.pos-var1->pos2.pos+1;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, var2);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][vi->variant.pos2.pos/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len-1)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
    oneStart(vi);
}
void varsearch::makeTrsInsInfoList(int i, int j, int v, int k) {
    auto &trs1 = summary->variants[5*i+3], &trs2 = summary->variants[5*j+3];
    auto &list = sr_trs[i][3*j+k];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInfo((k==1?&trs1[selement.first]:&trs2[selement.first]),
                    (k==1?&trs2[selement.second]:&trs1[selement.second]), vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeTrsInsCPInfoList(int i, int j, int v, int k) {
    auto &trs1 = summary->variants[5*i+3], &trs2 = summary->variants[5*j+3];
    auto &list = sr_trs[i][3*j+k];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInfo((k==1?&trs1[selement.first]:&trs2[selement.first]),
                    (k==1?&trs2[selement.second]:&trs1[selement.second]), vptr, summary, ref);
        calcInDelCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//TrsInv///////////////////////////
//          1 2     2 1     1 2  //
//         -> <-   -> <- = -> <- //
//           x       x           //
// <- -> = <- ->   <- ->         //
//  1'2'    2'1'    1'2'         //
///////////////////////////////////
void varsearch::makeTrsInvList(int i, int j) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j];
    auto &fidx = summary->trinvidx[i][2*j], &ridx = summary->trinvidx[i][2*j+1];
    if (!vpar->srv_par.detect_comp_var[3] || tinv.empty()) return;
    auto fbeg = tinv.begin()+fidx[0];
    auto rbeg = tinv.begin()+ridx[0];
    if (fbeg->pos2.idx != j || rbeg->pos2.idx != j || fbeg->pos1.dir || rbeg->pos2.dir) return;
    for (auto tit = fbeg; tit < tinv.end(); ++tit) {
        if (j < tit->pos2.idx || tit->pos1.dir) break;
        rbeg = tinv.begin()+ridx[(tit->pos1.pos-vpar->max_dist)>>14];
        for (auto tit_ = rbeg; tit_ < tinv.end(); ++tit_) {
            if (j < tit_->pos2.idx || tit_->pos2.dir) break;
            if (srVarCheck1(4, *tit, vpar) &&
                srVarCheck1(4, *tit_, vpar) &&
                srVarCheck2(3, *tit, *tit_, vpar) &&
                tit->pos1.pos-vpar->max_dist <= tit_->pos1.pos &&
                tit->pos2.pos-vpar->max_dist <= tit_->pos2.pos)
                list.add(std::make_pair(tit-tinv.begin(), tit_-tinv.begin()));
        }
    }
}
inline void makeTrsInvInfo(variant_t *var1, variant_t *var2, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    vi->variant = *var1;
    variant_t var_ = *var2;
    var_.comp();
    vi->variant.pos1.len = var_.pos2.pos-var1->pos1.pos-1;
    if (0 < vi->variant.pos1.len) vi->variant.pos1.pos++;
    else {
        vi->variant.pos1.pos = var_.pos2.pos;
        vi->variant.pos1.len = 0;
    }
    vi->variant.pos2.len = var_.pos1.pos-var1->pos2.pos-1;
    if (0 < vi->variant.pos2.len) vi->variant.pos2.pos++;
    else {
        vi->variant.pos2.pos = var_.pos1.pos;
        vi->variant.pos2.len = 0;
    }
    if (vi->variant.pos1.len || vi->variant.pos2.len)
        vi->variant.type |= SB_DELETION;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, &var_);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
}
void varsearch::makeTrsInvInfoList(int i, int j, int v) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInvInfo(&tinv[selement.first], &tinv[selement.second], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeTrsInvCPInfoList(int i, int j, int v) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInvInfo(&tinv[selement.first], &tinv[selement.second], vptr, summary, ref);
        calcTrsCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}

//TrsInvIns//////////////////
//  1 2'    2 1'      1'2  //
// -> <-   <- ->  =  <- -> //
//  1'2     2 1'      1 2' //
// -> <- = -> <-     <- -> //
/////////////////////////////
void varsearch::makeTrsInvInsList1(int i, int j) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j+1];
    auto &fidx = summary->trinvidx[i][2*j], &ridx = summary->trinvidx[i][2*j+1];
    if (!vpar->srv_par.detect_comp_var[0] || tinv.empty()) return;
    auto fbeg = tinv.begin()+fidx[0];
    auto rbeg = tinv.begin()+ridx[0];
    if (fbeg->pos2.idx != j || rbeg->pos2.idx != j || fbeg->pos1.dir || rbeg->pos2.dir) return;
    for (auto tit = fbeg; tit < tinv.end(); ++tit) {
        if (j < tit->pos2.idx || tit->pos1.dir) break;
        rbeg = tinv.begin()+ridx[(tit->pos1.pos-vpar->max_dist)>>14];
        for (auto tit_ = rbeg; tit_ < tinv.end(); ++tit_) {
            if (j < tit_->pos2.idx || tit_->pos2.dir) break;
            if (srVarCheck1(4, *tit, vpar) &&
                srVarCheck1(4, *tit_, vpar) &&
                srVarCheck2(0, *tit, *tit_, vpar) &&
                tit_->pos2.pos < tit->pos2.pos &&
                tit->pos1.pos-vpar->max_dist <= tit_->pos1.pos &&
                (vpar->min_length[0] <= tit_->pos1.pos-tit->pos1.pos ||
                 vpar->min_length[1] <= tit->pos2.pos-tit_->pos2.pos+1))
                list.add(std::make_pair(tit-tinv.begin(), tit_-tinv.begin()));
        }
    }
}
void varsearch::makeTrsInvInsList2(int i, int j) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j+2];
    auto &fidx = summary->trinvidx[i][2*j], &ridx = summary->trinvidx[i][2*j+1];
    if (!vpar->srv_par.detect_comp_var[0] || tinv.empty()) return;
    auto fbeg = tinv.begin()+fidx[0];
    auto rbeg = tinv.begin()+ridx[0];
    if (fbeg->pos2.idx != j || rbeg->pos2.idx != j || fbeg->pos1.dir || rbeg->pos2.dir) return;
    for (auto tit = fbeg; tit < tinv.end(); ++tit) {
        if (j < tit->pos2.idx || tit->pos1.dir) break;
        for (auto tit_ = rbeg; tit_ < tinv.end(); ++tit_) {
            if (j < tit_->pos2.idx || tit_->pos2.dir || tit->pos1.pos < tit_->pos1.pos) break;
            if (srVarCheck1(4, *tit, vpar) &&
                srVarCheck1(4, *tit_, vpar) &&
                srVarCheck2(0, *tit, *tit_, vpar) &&
                tit->pos2.pos-vpar->max_dist <= tit_->pos2.pos &&
                (vpar->min_length[0] <= tit_->pos2.pos-tit->pos2.pos ||
                 vpar->min_length[1] <= tit->pos1.pos-tit_->pos1.pos+1))
                list.add(std::make_pair(tit_-tinv.begin(), tit-tinv.begin()));
        }
    }
}
void makeTrsInvInsInfo(variant_t *var1, variant_t *var2, varinfo_t *vi, sbamsummary *sum, sdnafile *ref) {
    vi->var_method = SPLIT_READ_VARIANT;
    variant_t var_;
    if (var1->pos1.dir) {
        var_ = *var1;
        vi->variant = *var2;
        vi->variant.comp();
    }
    else {
        vi->variant = *var1;
        var_ = *var2;
        var_.comp();
    }
    vi->variant.type |= SB_INSERTION;
    vi->variant.pos1.len = var_.pos2.pos-vi->variant.pos1.pos-1;
    if (0 < vi->variant.pos1.len) {
        vi->variant.pos1.pos++;
        vi->variant.type |= SB_DELETION;
    }
    else {
        vi->variant.pos1.pos = var_.pos2.pos;
        vi->variant.pos1.len = 1;
    }
    vi->variant.pos2.len = vi->variant.pos2.pos-var_.pos1.pos+1;
    vi->variant.pos2.pos = var_.pos1.pos;
    vi->variant.pos1.setName(ref->name(vi->variant.pos1.idx));
    vi->variant.pos2.setName(ref->name(vi->variant.pos2.idx));
    combinate(vi, var2);
    vi->cp.frequency = 2.0*(vi->variant.nread+vi->variant.pread)/
    (sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos-1)/sum->bin]+
     sum->depth[vi->variant.pos1.idx][(vi->variant.pos1.pos+vi->variant.pos1.len)/sum->bin]+
     sum->depth[vi->variant.pos2.idx][vi->variant.pos2.pos/sum->bin]+
     sum->depth[vi->variant.pos2.idx][(vi->variant.pos2.pos+vi->variant.pos2.len-1)/sum->bin]);
    if (1.0f < vi->cp.frequency) vi->cp.frequency = 1.0;
    oneStart(vi);
    oneStart(vi);
}
void varsearch::makeTrsInvInsInfoList(int i, int j, int v, int k) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j+k];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInvInsInfo(&tinv[selement.first], &tinv[selement.second], vptr, summary, ref);
        ++vptr;
    }
}
void varsearch::makeTrsInvInsCPInfoList(int i, int j, int v, int k) {
    auto &tinv = summary->variants[5*i+4];
    auto &list = sr_trinvs[i][3*j+k];
    if (list.empty()) return;
    auto vptr = &vlist->at(v);
    sforeach(list) {
        makeTrsInvInsInfo(&tinv[selement.first], &tinv[selement.second], vptr, summary, ref);
        calcInDelCp(vptr, summary, background, cnan, vpar);
        ++vptr;
    }
}


void varsearch::searchSRVar() {
    sforeachi(i, *ref) {
        threads->addTask(&varsearch::makeDelList, this, i);
        threads->addTask(&varsearch::makeDupList, this, i);
        threads->addTask(&varsearch::makeInDelList1, this, i);
        threads->addTask(&varsearch::makeInDelList2, this, i);
        threads->addTask(&varsearch::makeInvList, this, i);
        threads->addTask(&varsearch::makeInvInsList1, this, i);
        threads->addTask(&varsearch::makeInvInsList2, this, i);
        sforin(j, 0, i) threads->addTask(&varsearch::makeTrsInsList2, this, i, j);
        sforin(j, i+1, ref->size()) {
            threads->addTask(&varsearch::makeTrsList, this, i, j);
            threads->addTask(&varsearch::makeTrsInsList1, this, i, j);
            threads->addTask(&varsearch::makeTrsInvList, this, i, j);
            threads->addTask(&varsearch::makeTrsInvInsList1, this, i, j);
            threads->addTask(&varsearch::makeTrsInvInsList2, this, i, j);
        }
    }
    threads->complete();
}
void varsearch::makeSRVInfo() {
    int voffset = 0;
    sforeachi(i, *ref) {
        threads->addTask(&varsearch::makeDelInfoList, this, i, voffset); voffset += sr_dels[i].size();
        threads->addTask(&varsearch::makeDupInfoList, this, i, voffset); voffset += sr_dups[i].size();
        threads->addTask(&varsearch::makeInDelInfoList, this, i, voffset, 0); voffset += sr_indels[2*i].size();
        threads->addTask(&varsearch::makeInDelInfoList, this, i, voffset, 1); voffset += sr_indels[2*i+1].size();
        threads->addTask(&varsearch::makeInvInfoList, this, i, voffset); voffset += sr_invs[3*i].size();
        threads->addTask(&varsearch::makeInvInsInfoList, this, i, voffset, 1); voffset += sr_invs[3*i+1].size();
        threads->addTask(&varsearch::makeInvInsInfoList, this, i, voffset, 2); voffset += sr_invs[3*i+2].size();
        sforeachi(j, *ref) {
            threads->addTask(&varsearch::makeTrsInfoList, this, i, j, voffset); voffset += sr_trs[i][3*j].size();
            threads->addTask(&varsearch::makeTrsInsInfoList, this, i, j, voffset, 1); voffset += sr_trs[i][3*j+1].size();
            threads->addTask(&varsearch::makeTrsInsInfoList, this, i, j, voffset, 2); voffset += sr_trs[i][3*j+2].size();
            threads->addTask(&varsearch::makeTrsInvInfoList, this, i, j, voffset); voffset += sr_trinvs[i][3*j].size();
            threads->addTask(&varsearch::makeTrsInvInsInfoList, this, i, j, voffset, 1); voffset += sr_trinvs[i][3*j+1].size();
            threads->addTask(&varsearch::makeTrsInvInsInfoList, this, i, j, voffset, 2); voffset += sr_trinvs[i][3*j+2].size();
        }
    }
    threads->complete();
}
void varsearch::makeCNSRVInfo() {
    int voffset = 0;
    sforeachi(i, *ref) {
        threads->addTask(&varsearch::makeDelCPInfoList, this, i, voffset); voffset += sr_dels[i].size();
        threads->addTask(&varsearch::makeDupCPInfoList, this, i, voffset); voffset += sr_dups[i].size();
        threads->addTask(&varsearch::makeInDelCPInfoList, this, i, voffset, 0); voffset += sr_indels[2*i].size();
        threads->addTask(&varsearch::makeInDelCPInfoList, this, i, voffset, 1); voffset += sr_indels[2*i+1].size();
        threads->addTask(&varsearch::makeInvCPInfoList, this, i, voffset); voffset += sr_invs[3*i].size();
        threads->addTask(&varsearch::makeInvInsCPInfoList, this, i, voffset, 1); voffset += sr_invs[3*i+1].size();
        threads->addTask(&varsearch::makeInvInsCPInfoList, this, i, voffset, 2); voffset += sr_invs[3*i+2].size();
        sforeachi(j, *ref) {
            threads->addTask(&varsearch::makeTrsCPInfoList, this, i, j, voffset); voffset += sr_trs[i][3*j].size();
            threads->addTask(&varsearch::makeTrsInsCPInfoList, this, i, j, voffset, 1); voffset += sr_trs[i][3*j+1].size();
            threads->addTask(&varsearch::makeTrsInsCPInfoList, this, i, j, voffset, 2); voffset += sr_trs[i][3*j+2].size();
            threads->addTask(&varsearch::makeTrsInvCPInfoList, this, i, j, voffset); voffset += sr_trinvs[i][3*j].size();
            threads->addTask(&varsearch::makeTrsInvInsCPInfoList, this, i, j, voffset, 1); voffset += sr_trinvs[i][3*j+1].size();
            threads->addTask(&varsearch::makeTrsInvInsCPInfoList, this, i, j, voffset, 2); voffset += sr_trinvs[i][3*j+2].size();
        }
    }
    threads->complete();
}
                             
void varsearch::search(int mode) {
    //if(mode&CNV_DETECTION || mode&CNSRV_DETECTION) cnan->search();
    if(mode&SRV_DETECTION || mode&CNSRV_DETECTION) {
        //if(background) summary->subtract(background, vpar);
        //else summary->varindex(&par->var_par);
        searchSRVar();
        int total = 0;
        sforeachi(i, *ref) {
            total += sr_dels[i].size();
            total += sr_dups[i].size();
            total += sr_indels[2*i].size();
            total += sr_indels[2*i+1].size();
            total += sr_invs[3*i].size();
            total += sr_invs[3*i+1].size();
            total += sr_invs[3*i+2].size();
            sforeachi(j, *ref) {
                total += sr_trs[i][3*j].size();
                total += sr_trs[i][3*j+1].size();
                total += sr_trs[i][3*j+2].size();
                total += sr_trinvs[i][3*j].size();
                total+= sr_trinvs[i][3*j+1].size();
                total += sr_trinvs[i][3*j+2].size();
            }
        }
        vlist->resize(total);
        if (mode&CNSRV_DETECTION) makeCNSRVInfo();
        else makeSRVInfo();
    }
    if(mode&CNV_DETECTION) cnan->makelist(vlist);
    vlist->reorder();
}

void varsearch::init() {
    //if(cnan) cnan->init();
    sforeachi(i, *ref) {
        sr_dels[i].clear();
        sr_dups[i].clear();
        sr_indels[2*i].clear();
        sr_indels[2*i+1].clear();
        sr_invs[3*i].clear();
        sr_invs[3*i+1].clear();
        sr_invs[3*i+2].clear();
        sforeachi(j, *ref) {
            sr_trs[i][3*j].clear();
            sr_trs[i][3*j+1].clear();
            sr_trs[i][3*j+2].clear();
            sr_trinvs[i][3*j].clear();
            sr_trinvs[i][3*j+1].clear();
            sr_trinvs[i][3*j+2].clear();
        }
    }
}