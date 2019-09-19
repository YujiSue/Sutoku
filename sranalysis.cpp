#include "analyzer.hpp"

using namespace stk;

split_read_t::split_read_t() { init(); }
split_read_t::split_read_t(const split_read_t &sr) {
    align = sr.align;
    head_clip = sr.head_clip;
    read_dir = sr.read_dir;
    seq = sr.seq;
}
split_read_t::~split_read_t() {}

split_read_t &split_read_t::operator=(const split_read_t &sr) {
    align = sr.align;
    head_clip = sr.head_clip;
    read_dir = sr.read_dir;
    seq = sr.seq;
    return *this;
}
void split_read_t::setHead(salign_t *a) {
    align = a;
    head_clip = false;
    size_t qend = a->pos+a->len;
    memmove(seq.ptr(), seq.ptr(qend), seq.length()-qend);
    seq.setSize(seq.length()-qend);
}
void split_read_t::setTail(salign_t *a) {
    align = a;
    head_clip = true;
    seq.setSize(a->pos);
}
inline void split_read_t::init() {
    align = nullptr;
    head_clip = false;
    read_dir = false;
    seq = sbioseq(SB_DNA_SEQ, "");
}

inline void toVar(variant_t *var, split_read_t *read, salign_t *align, int &&score, int &total) {
    if (read->read_dir) var->nread = 1;
    else var->pread = 1;
    if (score < 0) {
        std::cout<<"score error"<<std::endl;
    }
    var->qual = 1.0-(double)score/total;
    int beg, end;
    if (read->head_clip) {
        var->pos1 = align->ref;
        var->pos2 = read->align->ref;
        if (var->pos1.dir) beg = read->seq.length()-align->pos;
        else beg = align->pos+align->len;
        end = read->seq.length();
    }
    else {
        var->pos2 = align->ref;
        var->pos1 = read->align->ref;
        beg = 0;
        if (var->pos2.dir)
            end = read->seq.length()-(align->pos+align->len);
        else end = align->pos;
    }
    if (beg < end)
        var->alt = read->seq.raw().substring(beg, end-beg);
    if (var->pos1.idx == var->pos2.idx) {
        if(!align->ref.dir) {
            var->pos1.pos += var->pos1.len-1;
            if(var->pos1.pos < var->pos2.pos) {
                var->type = SB_DELETION;
                var->pos1.len = var->pos2.pos-var->pos1.pos-1;
            }
            else {
                var->type = SB_DUPLICATION;
                var->pos1.len = var->pos1.pos-var->pos2.pos+1;
            }
        }
        else {
            var->type = SB_INVERSION;
            if (!(read->head_clip)) {
                var->pos1.pos += var->pos1.len-1;
                var->pos2.pos += var->pos2.len-1;
            }
            var->pos1.len = 0;
            var->pos2.len = 0;
        }
    }
    else {
        if(!align->ref.dir) {
            var->type = SB_TRANSLOCATION;
            var->pos1.pos += var->pos1.len-1;
            var->pos1.len = 0;
            var->pos2.len = 0;
        }
        else {
            var->type = SB_INVERSION|SB_TRANSLOCATION;
            if (!(read->head_clip)) {
                var->pos1.pos += var->pos1.len-1;
                var->pos2.pos += var->pos2.len-1;
            }
            var->pos1.len = 0;
            var->pos2.len = 0;
        }
    }
    if(var->type&SB_INVERSION &&
       (var->pos2.idx < var->pos1.idx ||
        (var->pos1.idx == var->pos2.idx &&
         var->pos2.pos < var->pos1.pos))) var->comp();
}
inline salign_t *bestAlign(split_read_t *sr, slist<salign_t *> *aligns) {
    salign_t *best = aligns->at(0);
    if (1 < aligns->size()) {
        auto ait = aligns->begin()+1;
        while (ait < aligns->end()) {
            if ((*ait)->ref.idx == sr->align->ref.idx) {
                if (best->ref.idx == sr->align->ref.idx) {
                    if ((*ait)->ref.dir == best->ref.dir) {
                        int dist1, dist2;
                        if ((*ait)->ref.dir) {
                            if (sr->head_clip) {
                                dist1 = smath::abs(sr->align->ref.pos - (*ait)->ref.pos);
                                dist2 = smath::abs(sr->align->ref.pos - best->ref.pos);
                            }
                            else {
                                dist1 = smath::abs((*ait)->ref.pos+(*ait)->ref.len-
                                                   (sr->align->ref.pos+sr->align->ref.len));
                                dist2 = smath::abs(best->ref.pos+best->ref.len-
                                                   (sr->align->ref.pos+sr->align->ref.len));
                            }
                        }
                        else {
                            if (sr->head_clip) {
                                dist1 = smath::abs(sr->align->ref.pos - ((*ait)->ref.pos+(*ait)->ref.len));
                                dist2 = smath::abs(sr->align->ref.pos - (best->ref.pos+best->ref.len));
                            }
                            else {
                                dist1 = smath::abs((*ait)->ref.pos - (sr->align->ref.pos+sr->align->ref.len));
                                dist2 = smath::abs(best->ref.pos - (sr->align->ref.pos+sr->align->ref.len));
                            }
                        }
                        if (dist1 < dist2) best = *ait;
                    }
                    else if (!(*ait)->ref.dir) best = *ait;
                }
                else best = *ait;
            }
            ++ait;
        }
    }
    return best;
}

sranalysis::sranalysis(analyzer *a) {
    par = a->par;
    ref = par->ref;
    bam = &a->bam;
    summary = a->summary;
    status = a->status;
    threads = &a->threads;
    rivs.resize(2*ref->size());
    srvs.resize(ref->size());
    ques.resize(ref->size());
    searches.resize(ref->size());
    extends.resize(ref->size());
    sforeachi(i, *ref) {
        rivs[2*i].resize(par->an_par.max_clip_count);
        rivs[2*i+1].resize(par->an_par.max_clip_count);
        srvs[i].resize(par->an_par.max_clip_count);
        ques[i] = new sbquery(&par->ls_par);
        searches[i] = new slocalsearch(&par->ls_par);
        extends[i] = new slocalextend(&par->ls_par);
    }
    size_t size;
    if (par->bam_par.parallel) size = ref->size();
    else size = 1;
    mtx = std::vector<std::mutex>(ref->size());
}
sranalysis::~sranalysis() {
    sforeachi(i, *ref) {
        if(ques[i]) delete ques[i];
        if(searches[i]) delete searches[i];
        if(extends[i]) delete extends[i];
    }
}

int stk::srvtype(variant_t *var) {
    switch (var->type) {
        case SB_DUPLICATION:
            return 1;
            break;
        case SB_INVERSION:
            return 2;
            break;
        case SB_TRANSLOCATION:
            return 3;
            break;
        case SB_TRANSLOCATION|SB_INVERSION:
            return 4;
            break;
        default:
            return 0;
            break;
    }
}
inline int alignScore(salign_t *al, int32_t &len, salignment_param_t *par) {
    auto score = al->score;
    if (al->pos) score += par->gap_score+(2<al->pos?(al->pos-1)*par->gap2_score:0);
    if (al->pos+al->len < len) {
        auto l = len - (al->pos+al->len);
        score += par->gap_score+(2<l?(l-1)*par->gap2_score:0);
    }
    return score<=0?1:score;
}
void makesrv(int r, sr_vec *sv, sbquery *sbq, slocalsearch *sls, slocalextend *sle,
                    sbamsummary *summary, stk_param_t *par) {
    auto count = sbq->count();
    sbq->makeTrie();
    sls->search(par->ref, sbq, false);
    slist<salign_t *> aligns;
    int score = 0;
    sforin(q, 0, count) {
        aligns.clear();
        auto &searched = sls->aligns[q];
        sforeachi(i, *(par->ref)) {
            auto &result = searched[i];
            sforeach(result) {
                score += alignScore(&selement, sbq->seq(q, false).length(), &par->ls_par.aln_par);
                if (aligns.empty() || aligns[0]->score == selement.score) aligns.add(&selement);
                else {
                    aligns.clear();
                    aligns.add(&selement);
                }
            }
        }
        /*
        for (int j = 0; j < par->ref->size(); ++j) {
            auto &res = sls->result()[q][j];
            for (auto ait = res.begin(); ait < res.end(); ++ait) {
                sle->extend(par->ref->at(j), &sbq->seq(q, ait->ref.dir), &(*ait));
                score += ait->score;
                if (aligns.empty() || aligns[0]->score == ait->score) aligns.add(&(*ait));
                else {
                    aligns.clear();
                    aligns.add(&(*ait));
                }
            }
        }
         */
        if (!aligns.empty()) {
            salign_t *best = bestAlign(&sv->at(q), &aligns);
            variant_t var;
            toVar(&var, &sv->at(q), best,
                  alignScore(best, sbq->seq(q, false).length(), &par->ls_par.aln_par), score);
            mapping(summary, best, true);
            addSR(summary, &var, true);
        }
    }
    sls->reset();
    sbq->reset();
}
void clipextend(sbreadinfo *ri, split_read_t *sr, sbquery *sbq, slocalextend *sle,
                sbamsummary *summary, sdnafile *ref, status_t *status, int length) {
    try {
        
        if (ri->align.ref.idx == 0 && srange(183364, 208465).include(ri->align.ref.pos)) {
            std::cout<<ri->align.cigars.toString()<<":"<<ri->sequence.raw()<<std::endl;
        }
        
        if (ri->headclip(length)) {
            ri->align.pos += ri->align.cigars.first().length;
            ri->align.len -= ri->align.cigars.first().length;
            ri->align.cigars.trimhead();
        }
        if (ri->tailclip(length)) {
            ri->align.len -= ri->align.cigars.last().length;
            ri->align.cigars.trimtail();
        }
        sr->read_dir = ri->flag&0x10;
        sr->seq.setSize(ri->sequence.length());
        DNA_CONVERTER[2][1](ri->sequence.ptr(), 0, ri->sequence.length(), sr->seq.ptr());
        sle->extend(ref->at(ri->align.ref.idx), &sr->seq, &ri->align);
        mapping(summary, &ri->align, true);
        if (length <= ri->align.pos) {
            sr->setTail(&ri->align);
            sbq->addQuery(sr->seq);
        }
        else if(ri->align.pos+ri->align.len+length <= ri->sequence.length()) {
            sr->setHead(&ri->align);
            sbq->addQuery(sr->seq);
        }
        else sbq->addQuery("");
    } catch (sbioexception be) {
        status->abort(&be);
    }
}
void clipreads(int r, sbamfile *bam, ri_vec *rv, clip_vec *cv, status_t *status, int32_t max) {
    sint &current = status->current_task[r];
    int count = 0;
    while (status->analysis && current < cv->size() && count < max) {
        try {
            bam->setVOffset(cv->at(current));
            bam->next(&rv->at(count));
            ++current;
            ++count;
        } catch (sbioexception be) {
            status->abort(&be);
            break;
        }
    }
    if (count < max) rv->resize(count);
}

void sranalysis::sranalyze(int r, clip_vec *cv) {
    auto sbq = ques[r];
    auto sls = searches[r];
    auto sle = extends[r];
    stpool bufthread(1);
    auto file = bam->at(r);
    int max = par->an_par.max_clip_count, len = par->an_par.min_clip_length;
    ri_vec *reads = &rivs[2*r], *buf_reads = &rivs[2*r+1];
    sr_vec *sreads = &srvs[r];
    clipreads(r, file, reads, cv, status, max);
    bufthread.addTask(clipreads, r, file, buf_reads, cv, status, max);
    while (status->analysis && status->current_task[r] < cv->size()) {
        sreads->resize(reads->size());
        auto rit = reads->begin();
        auto sit = sreads->begin();
        while (status->analysis && rit < reads->end()) {
            clipextend(&(*rit), &(*sit), sbq, sle, summary, par->ref, status, len);
            ++rit; ++sit;
        }
        makesrv(r, sreads, sbq, sls, sle, summary, par);
        bufthread.complete();
        auto tmp = reads; reads  = buf_reads; buf_reads = tmp;
        bufthread.addTask(clipreads, r, file, buf_reads, cv, status, max);
    }
}

void sranalysis::analyze(slist<clip_vec> &cv) {
    sforeachi(i, cv) { sranalyze(i, &cv[i]); }
    /*
    status->current_task.resize(ref->size(), 0);
    if (par->bam_par.parallel) {
        sforeachi(i, cv) threads->addTask(&sranalysis::sranalyze, this, i, &cv[i]);
        threads->complete();
    }
    else sforeachi(i, cv) { sranalyze(i, &cv[i]); }
     */
}

void pair2var(int i, int j, sbamsummary *summary) {
    auto &vec = summary->variants[5*i+j];
    sforeach(vec) {
        if (selement.type&SB_INVERSION) {
            if (selement.pos1.dir) {
                selement.nread = 1;
                selement.pos2.pos += selement.pos2.len-1;
            }
            else {
                selement.pread = 1;
                selement.pos1.pos += selement.pos2.len-1;
            }
            selement.pos1.len = 0;
            selement.pos2.len = 0;
            if ((selement.type&SB_TRANSLOCATION && selement.pos2.idx < selement.pos1.idx) ||
                (selement.type==SB_INVERSION && selement.pos2.pos < selement.pos1.pos)) {
                auto p = selement.pos1;
                selement.pos1 = selement.pos2;
                selement.pos2 = p;
                selement.pos1.dir = !selement.pos1.dir;
                selement.pos2.dir = !selement.pos2.dir;
                if (selement.nread) { selement.nread = 0; selement.pread = 1; }
                else { selement.pread = 0; selement.nread = 1; }
            }
        }
        else {
            if (selement.pos1.dir) {
                selement.nread = 1;
                selement.pos1.dir = false;
                selement.pos2.dir = false;
            }
            else selement.pread = 1;
            selement.pos1.pos += selement.pos1.len;
            if (selement.type&SB_TRANSLOCATION) {
                selement.pos1.len = 0;
                selement.pos2.len = 0;
            }
            else {
                if (selement.type&SB_DELETION)
                    selement.pos1.len = selement.pos2.pos-selement.pos1.pos;
                else {
                    selement.pos1.len = selement.pos1.pos-selement.pos2.pos;
                    --selement.pos1.pos;
                }
                selement.pos2.len = 0;
            }
        }
    }
    auto size = vec.size();
    std::sort(vec.begin(), vec.end(), [](const variant_t &v1, const variant_t &v2) { return v1.pos1 < v2.pos1; });
    sforeach(vec) {
        if (!(selement.type)) continue;
        auto it_ = it+1;
        while (it_ < vec.end()) {
            if (it_->type != 0) {
                if (it->pos1.idx != it_->pos1.idx || it->pos1.dir != it_->pos1.dir ||
                    it->pos2.idx != it_->pos2.idx || it->pos2.dir != it_->pos2.dir ||
                    it->pos1.pos+200 < it_->pos1.pos) break;
                if (abs(it->pos1.pos-it_->pos1.pos) < 200 && abs(it->pos2.pos-it_->pos2.pos) < 200) {
                    it->pread+=it_->pread; it->nread += it_->nread; --size; it_->type = 0;
                }
            }
            ++it_;
        }
    }
    std::sort(vec.begin(), vec.end(), [](const variant_t &v1, const variant_t &v2) {
        if (!(v1.type)) return false;
        if (!(v2.type)) return true;
        return v1 < v2;
    });
    vec.resize(size);
}

void sranalysis::analyze() {
    for (int i = 0; i < ref->size(); ++i) {
        for (int j = 0; j < 5; ++j) pair2var(i, j, summary);
    }
    //threads->complete();
}
void sranalysis::init() {
    sforeachi(i, *ref)  {
        rivs[i].resize(par->an_par.max_clip_count);
        srvs[i].resize(par->an_par.max_clip_count);
        ques[i]->reset();
        searches[i]->reset();
        extends[i]->reset();
    }
}