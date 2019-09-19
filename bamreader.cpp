#include "analyzer.hpp"

using namespace stk;

bamreader::bamreader(analyzer *a) {
    par = a->par;
    status = a->status;
    ref = par->ref;
    clips.resize(ref->size());
    bam = &a->bam;
    summary = a->summary;
    threads = &a->threads;
}

bamreader::~bamreader() {}

void bamreader::read(slist<sbreadinfo> *list) {
    auto b = bam->first();
    sforeachi(i, *ref) {
        if (par->target.empty())
            b->reads(list, i, slib::srange(0, ref->length(i)));
        else b->reads(list, i, par->target[i]);
    }
}

void searchSRFromPaired(sbamfile *bam, sbamsummary *summary, slist<clip_vec> *clips,
                status_t *status, analyze_param_t *par) {
    try {
        sbreadinfo ri;
        auto init = bam->voffset().file_offset;
        while (status->analysis && bam->next(&ri)) {
            status->current_task[0] = bam->voffset().file_offset-init;
            if(ri.flag==4) {
                ++summary->total_reads;
                continue;
            }
            int32_t &r = ri.align.ref.idx;
            ++summary->read_count[r];
            summary->read_length[r]+=ri.sequence.length();
            if (!(ri.flag&0x02) && (ri.flag&0x01) && (ri.flag&0x40)) {
                auto inv_flag = (ri.flag&0x10)|(ri.flag&0x20);
                if (ri.next_refid != r) {
                    if (inv_flag == 0 || inv_flag == 0x30) {
                        variant_t v;
                        v.type = SB_TRANSLOCATION|SB_INVERSION;
                        v.pos1 = ri.align.ref;
                        v.pos1.len = ri.sequence.length();
                        v.pos2.idx = ri.next_refid;
                        v.pos2.pos = ri.next_pos;
                        v.pos2.dir = !(ri.flag&0x20);
                        summary->variants[5*r+4].add(v);
                    }
                    else {
                        variant_t v;
                        v.type = SB_TRANSLOCATION;
                        v.pos1 = ri.align.ref;
                        v.pos1.len = ri.sequence.length();
                        v.pos2.idx = ri.next_refid;
                        v.pos2.pos = ri.next_pos;
                        summary->variants[5*r+3].add(v);
                    }
                }
                else if (inv_flag == 0 || inv_flag == 0x30) {
                    variant_t v;
                    v.type = SB_INVERSION;
                    v.pos1 = ri.align.ref;
                    v.pos1.len = ri.sequence.length();
                    v.pos2.idx = ri.next_refid;
                    v.pos2.pos = ri.next_pos;
                    v.pos2.dir = !(ri.flag&0x20);
                    summary->variants[5*r+2].add(v);
                }
                else {
                    variant_t v;
                    v.pos1 = ri.align.ref;
                    v.pos1.len = ri.sequence.length();
                    v.pos2.idx = ri.next_refid;
                    v.pos2.pos = ri.next_pos;
                    
                    if (v.pos1.pos < v.pos2.pos) {
                        v.type = SB_DELETION;
                        if (v.pos2.pos-v.pos1.pos-ri.sequence.length() > 200) {
                            summary->variants[5*r].add(v);
                        }
                    }
                    else {
                        v.type = SB_DUPLICATION;
                        summary->variants[5*r+1].add(v);
                    }
                    
                }
            }
            mapping(summary, &ri.align);
        }
    }
    catch(sbioexception be) {
        be.print();
        status->abort(&be);
    }
}

void searchClip(sbamfile *bam, sbamsummary *summary, slist<clip_vec> *clips,
                status_t *status, analyze_param_t *par) {
    try {
        sbreadinfo ri;
        auto init = bam->voffset().file_offset;
        while (status->analysis && bam->next(&ri)) {
            status->current_task[0] = bam->voffset().file_offset-init;
            if(ri.flag==4) {
                ++summary->total_reads;
                continue;
            }
            int32_t &r = ri.align.ref.idx;
            ++summary->read_count[r];
            summary->read_length[r]+=ri.sequence.length();
            if(ri.headclip(par->min_clip_length)||
               ri.tailclip(par->min_clip_length)) {
                clips->at(r).add(ri.offset.intOffset());
                continue;
            }
            mapping(summary, &ri.align);
        }
    }
    catch(sbioexception be) {
        be.print();
        status->abort(&be);
    }
}

void multiSearchClip(int r, sbvoffset *end, sbamfile *bam, sbamsummary *summary,
                     clip_vec *clips, status_t *status, analyze_param_t *par) {
    try {
        sbreadinfo ri;
        auto init = bam->voffset().file_offset;
        while (status->analysis && bam->voffset() < *end && bam->next(&ri)) {
            status->current_task[r] = bam->voffset().file_offset-init;
            if(ri.flag==4) {
                ++summary->total_reads;
                continue;
            }
            ++summary->read_count[r];
            summary->read_length[r]+=ri.sequence.length();
            if(ri.headclip(par->min_clip_length)||
               ri.tailclip(par->min_clip_length)) {
                clips->add(ri.offset.intOffset());
                continue;
            }
            mapping(summary, &ri.align);
        }
    }
    catch(sbioexception be) {
        be.print();
        status->abort(&be);
    }
}

void bamreader::summarize(bool paired) {
    if (!paired) {
        if (par->bam_par.parallel && bam->at(0)->hasIndex()) {
            status->current_task.resize(ref->size());
            status->current_task.reset(0);
            auto index = bam->first()->index;
            slist<sbvoffset> end(ref->size());
            sforeachi(i, *ref)  {
                if (0 < i) {
                    bam->at(i)->open(bam->first()->path());
                    bam->at(i)->setVOffset(index->lin_offset[i][0]);
                }
                if (i < ref->size()-1) end[i] = index->lin_offset[i+1][0];
            }
            end.last() = sbvoffset(bam->at(0)->size(), 0);
            sforeachi(j, *ref)
            threads->addTask(multiSearchClip, j, &end[j], bam->at(j), summary,
                             &clips[j], status, &par->an_par);
            threads->complete();
        }
        else {
            status->current_task.resize(1);
            status->current_task[0] = 0;
            searchClip(bam->first(), summary, &clips, status, &par->an_par);
        }
    }
    else {
        searchSRFromPaired(bam->first(), summary, &clips, status, &par->an_par);
        
        
    }
}

void bamreader::init() {
    sforeachi(i, clips) clips[i].clear();
}