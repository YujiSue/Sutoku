#include "analyzer.hpp"

using namespace slib;
using namespace slib::sio;
using namespace slib::sbio;
using namespace stk;

inline void loadRegion(const char *path, slist<sregion> *region, sbioseqlist *ref) {
    sfile file(path, SIO_READ);
    if (!file.size()) return;
    sstring row;
    region->resize(ref->size());
    while (!file.eof()) {
        file.readLine(row);
        auto vals = row.split("\t");
        if(vals.size() != 3) continue;
        srange range(vals[1].intValue(), vals[2].intValue());
        if(-1 < ref->indexOf(vals[0]))
            region->at(ref->indexOf(vals[0])).append(range);
    }
}

analyze_param_t::analyze_param_t() {
    thread_count = DEFAULT_THREAD_COUNT;
    cp_bin = DEFAULT_CP_BIN;
    max_clip_count = DEFAULT_CLIP_COUNT;
    min_clip_length = DEFAULT_CLIP_LENGTH;
}
analyze_param_t::~analyze_param_t() {}

void analyze_param_t::set(sobj_ptr &obj) {
    if (obj["thread"]) thread_count = obj["thread"];
    if (obj["bin"]) cp_bin = obj["bin"];
    if (obj["clip"]) max_clip_count = obj["clip"];
    if (obj["clen"]) min_clip_length = obj["clen"];
}

sobj_ptr analyze_param_t::toObj() { return sdict(
    {
        kv("thread", thread_count),
        kv("bin", cp_bin),
        kv("clip", max_clip_count),
        kv("clen", min_clip_length)
    });
}

annot_param_t::annot_param_t() {
    memset(var_filter, 0, 6*sizeof(bool));
}
annot_param_t::~annot_param_t() {}

void annot_param_t::set(sobj_ptr &obj) {
    if(!(obj["cds_var"].isNull())) var_filter[0] = obj["cds_var"];
    if(!(obj["exon_var"].isNull())) var_filter[1] = obj["exon_var"];
    if(!(obj["intr_var"].isNull())) var_filter[2] = obj["intr_var"];
    if(!(obj["utr_var"].isNull())) var_filter[3] = obj["utr_var"];
    if(!(obj["intg_var"].isNull())) var_filter[4] = obj["intg_var"];
}
sobj_ptr annot_param_t::toObj() { return sdict(
    {
        kv("cds_var", var_filter[0]),
        kv("exon_var", var_filter[1]),
        kv("intr_var", var_filter[2]),
        kv("utr_var", var_filter[3]),
        kv("intg_var", var_filter[4])
    });
}

stk_param_t::stk_param_t() {
    ref = nullptr;
    db = nullptr;
    bg = nullptr;
    outformat = "plist";
    ls_par.ref_type = COMPRESS4|SB_DNA_SEQ;
    ls_par.aln_par.seq_type = SB_DNA_SEQ;
    ls_par.aln_par.makeScore();
    
    bam_par.parallel = true;
    an_par.max_clip_count = 5000;
}
stk_param_t::~stk_param_t() {
    if(ref) delete ref;
    if(db) delete db;
    if(bg) delete bg;
}

void stk_param_t::set(sobj_ptr &obj) {
    if(obj["analyze"]) an_par.set(obj["analyze"]);
    if(obj["bam"]) bam_par.set(obj["bam"]);
    if(obj["search"]) ls_par.set(obj["search"]);
    if(obj["variant"]) var_par.set(obj["variant"]);
    if(obj["annot"]) ant_par.set(obj["annot"]);
}

sobj_ptr stk_param_t::toObj() { return sdict(
    {
        kv("analyze", an_par.toObj()),
        kv("bam", bam_par.toObj()),
        kv("search", ls_par.toObj()),
        kv("variant", var_par.toObj()),
        kv("annot", ant_par.toObj())
    });
}
void stk_param_t::load(const char *path) {
    sdict dic;
    dic->load(path);
    set(dic);
}
void stk_param_t::saveA(const char *path) {
    sdictionary dic = {
        kv("analyze", an_par.toObj()),
        kv("bam", bam_par.toObj()),
        kv("search", ls_par.toObj())
    };
    dic.save(path);
}
void stk_param_t::saveV(const char *path) {
    sdictionary dic = {
        kv("variant", var_par.toObj()),
        kv("annot", ant_par.toObj())
    };
    dic.save(path);
}
void stk_param_t::setRefPath(const char *path) {
    if (!ref) ref = new sdnafile();
    ref->load(path);
}
void stk_param_t::setDBPath(const char *path) {
    if (!db) db = new sbannotdb(TRANSCRIPT_INFO|MUT_INFO);
    db->load(path);
}
void stk_param_t::setBGPath(const char *path) {
    if (!bg) bg = new sbamsummary();
    bg->load(path);
}
void stk_param_t::setCtrlPath(const char *path) {
    if (!ctrl) ctrl = new varlist();
    ctrl->load(path);
}
void stk_param_t::setTgtPath(const char *path) {
    if (!ref) {
        std::cout<<"Reference file is required for setting target region."<<std::endl;
        return;
    }
    loadRegion(path, &target, ref);
}
void stk_param_t::setNtgtPath(const char *path) {
    if (!ref) {
        std::cout<<"Reference file is required for setting ntarget region."<<std::endl;
        return;
    }
    loadRegion(path, &ntarget, ref);
}
void stk::mapping(sbamsummary *sum, salign_t *al, bool sync) {
    int32_t refid = al->ref.idx, rpos = al->ref.pos, &bin = sum->bin;
    if (sync) sum->lock(refid);
    float *dp = &sum->depth[refid][0];
    for (auto it = al->cigars.begin(); it < al->cigars.end(); ++it) {
        int beg = rpos/bin, end = (rpos+it->length-1)/bin;
        if(it->option == SB_CIGAR_MATCH || it->option == SB_CIGAR_PMATCH || it->option == SB_CIGAR_MMATCH) {
            for (int i = beg; i < end; ++i) dp[i]+=1.0;
            dp[beg]-=(float)(rpos%bin)/bin;
            dp[end]+=(float)(rpos+it->length-end*bin)/bin;
            rpos+=it->length;
        }
        else if(it->option == SB_CIGAR_DELETION || it->option == SB_CIGAR_SKIP) {
            dp[beg]+=(float)(rpos%bin)/bin;
            dp[end]+=1.0-(float)(rpos+it->length-end*bin)/bin;
            rpos += it->length;
        }
    }
    if (sync) sum->unlock(refid);
}
void stk::addSR(sbamsummary *sum, variant_t *var, bool sync) {
    int refid = var->pos1.idx, vtype = srvtype(var);
    if (sync) sum->lock(var->pos1.idx, vtype);
    sum->variants[5*refid+vtype].add(*var);
    if (sync) sum->unlock(var->pos1.idx, vtype);
}

status_t::status_t() { init(); }
status_t::~status_t() {}

inline void status_t::abort(sexception *ex) {
    std::unique_lock<std::mutex> (mtx);
    err = ex->err;
    src = ex->src;
    msg = ex->msg;
    analysis = false;
}
inline double status_t::progress() {
    return (double)smath::sum(current_task.begin(), current_task.end())/total_task;
}
inline void status_t::init() {
    analysis = true;
    err = 0;
    msg = "";
}

analyzer::analyzer() {
    status = new status_t();
}
analyzer::analyzer(stk_param_t *p) : analyzer() {
    setParam(p);
    reader = new bamreader(this);
    sran = new sranalysis(this);
    cnan = new cnanalysis(this);
    search = new varsearch(this);
}

analyzer::~analyzer() {
    for (auto it = bam.begin(); it < bam.end(); ++it) {
        if(*it) delete *it;
    }
    if(status) delete status;
    if(summary) delete summary;
    if(variants) delete variants;
    if(reader) delete reader;
    if(sran) delete sran;
    if(cnan) delete cnan;
    if(search) delete search;
}

void analyzer::readview(const sfile &infile, sfile &outfile) {
    try {
        slist<sbreadinfo> list;
        std::cout<<"Loading bam : '"<<infile<<"'"<<std::endl;
        bam[0]->load(infile);
        reader->read(&list);
        if (status->err) return;
        std::cout<<"Exporting result : '"<<outfile<<"'"<<std::endl;
        for (auto it = list.begin(); it < list.end(); ++it) {
            outfile<<it->toString()<<NEW_LINE;
            outfile.flush();
        }
    }
    catch (sioexception ie) { status->abort(&ie); }
    catch (sbioexception be) { status->abort(&be); }
}
inline void summarizeThread(bamreader *reader, bool paired) {
    reader->summarize(paired);
}
inline void searchSRThread(sranalysis *sran, slist<clip_vec> *cv) {
    //sran->analyze(*cv);
    sran->analyze();
}
void analyzer::summarize(const sfile &infile, sfile &outfile) {
    try {
        std::cout<<"Loading bam : '"<<infile<<"'"<<std::endl;
        bam[0]->load(infile);
        status->total_task = bam[0]->size();
        sfile bai = infile+".bai";
        if (bai.exist()) {
            std::cout<<"Loading bam index : '"<<bai<<"'"<<std::endl;
            bam[0]->loadIndex(bai);
            if (!bam[0]->hasIndex())
                std::cout<<" Failed."<<std::endl;
        }
        std::cout<<"Reading...    0%"<<std::flush;
        stpool thread(1);
        thread.addTask(summarizeThread, reader, true);
        while (thread.isWorking()) {
            std::cout<<"\b\b\b\b"<<snumber((int)(status->progress()*100.0)).filled(3, ' ')<<"%"<<std::flush;
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
        thread.complete();
        if (status->err) {
            std::cout<<std::endl;
            return;
        }
        else std::cout<<"\b\b\b\b100%"<<std::endl;
        sint cl_num = 0;
        for (auto cit = reader->clips.begin(); cit < reader->clips.end(); ++cit) cl_num += cit->size();
        std::cout<<"Clip read count : "<<cl_num<<std::endl;
        status->total_task = cl_num;
        status->current_task.reset(0);
        std::cout<<"Clip read analyzing...    0%"<<std::flush;
        thread.addTask(searchSRThread, sran, &reader->clips);
        while (thread.isWorking()) {
            std::cout<<"\b\b\b\b"<<snumber((int)(status->progress()*100.0)).filled(3, ' ')<<"%"<<std::flush;
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
        thread.complete();
        if (status->err) {
            std::cout<<std::endl;
            return;
        }
        else std::cout<<"\b\b\b\b100%"<<std::endl;
        if (status->err) return;
        summary->integrate(&par->var_par);
        std::cout<<"Exporting summary : '"<<outfile<<"'"<<std::endl;
        summary->save(outfile);
    }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::sumintegrate(const slist<sfile> &infiles, sfile &outfile) {
    try {
        for (auto it = infiles.begin(); it < infiles.end(); ++it) {
            if (it->extension() != "bsm") continue;
            std::cout<<"Loading summary : '"<<*it<<"'"<<std::endl;
            summary->load(*it); ++it;
            int count = 0;
            while (it < infiles.end()) {
                if (lower((*it).extension()) == "bsm") {
                    sbamsummary sum;
                    std::cout<<"Loading summary : '"<<*it<<"'"<<std::endl;
                    sum.load((*it));
                    std::cout<<"Merging..."<<std::endl;
                    summary->merge(&sum, &par->var_par);
                    ++count;
                }
                ++it;
            }
            std::cout<<"Merged "<<count<<" files."<<std::endl;
        }
        std::cout<<"Exporting : '"<<outfile<<"'"<<std::endl;
        summary->save(outfile);
    }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::copynum(const sfile &infile, sfile &outfile) {
    try {
        std::cout<<"Loading summary : '"<<infile<<"'"<<std::endl;
        summary->load(infile);
        outfile.make();
        auto &bin = par->an_par.cp_bin;
        auto max = 0;
        slist<vecNd<float>> cpvalues;
        slist<sbpos_t> locations;
        if (par->target.size()) {
            for (int i = 0; i <  par->ref->size(); ++i) {
                if(par->target[i].empty()) continue;
                for (auto it = par->target[i].begin(); it < par->target[i].end(); ++it) {
                    sbpos_t pos(i, it->begin, it->length(), false);
                    pos.setName(par->ref->name(i));
                    locations.add(pos);
                    size_t len = (pos.len-1)/bin+1;
                    cpvalues.add(vecNd<float>(len, 0.0f));
                    if (max < len) max = len;
                    auto offset = pos.pos-1;
                    for (int l = 0; l < len; ++l) {
                        if (!summary->average_depth) cpvalues[-1][l] = -1;
                        else {
                            float bcp = 1.0;
                            if (par->bg && par->bg->average_depth)
                                bcp = par->bg->depthIn(i, offset, bin)/par->bg->average_depth;
                            cpvalues[-1][l] = summary->depthIn(i, offset, bin)/summary->average_depth/bcp;
                        }
                        offset += bin;
                    }
                }
            }
        }
        else {
            cpvalues.resize(par->ref->size());
            locations.resize(par->ref->size());
            for (int i = 0; i < par->ref->size(); ++i) {
                size_t len = par->ref->length(i);
                locations[i] = sbpos_t(i, 1, len, false);
                locations[i].setName(par->ref->name(i));
                len = (len-1)/bin+1;
                cpvalues[i].resize(len);
                if (max < len) max = len;
                if (!summary->average_depth) for (int l = 0; l < len; ++l) { cpvalues[i][l] = -1; }
                else {
                    
                    
                }
                auto offset = 0;
                for (int l = 0; l < len-1; ++l) {
                    if (!summary->average_depth) cpvalues[i][l] = -1;
                    else {
                        float bcp = 1.0;
                        if (par->bg && par->bg->average_depth)
                            bcp = par->bg->depthIn(i, offset, bin)/par->bg->average_depth;
                        cpvalues[i][l] = summary->depthIn(i, offset, bin)/summary->average_depth/bcp;
                    }
                    offset += bin;
                }
                
            }
        }
        std::cout<<"Exporting : '"<<outfile<<"'"<<std::endl;
        for (auto it = locations.begin(); it < locations.end(); ++it)
            outfile<<it->pos_id<<":"<<it->pos<<".."<<it->pos+it->len-1<<"\t\t";
        outfile<<NEW_LINE;
        outfile.flush();
        for (int r = 0; r < max; ++r) {
            for (int c = 0; c < locations.size(); ++c) {
                if (locations[c].len <= 0) outfile<<"\t\t";
                else {
                    outfile<<locations[c].pos<<".."<<locations[c].pos+bin<<"\t"<<
                    snumber(cpvalues[c][r]).precised(2)<<"\t";
                    locations[c].pos+=bin;
                    locations[c].len-=bin;
                }
            }
            outfile<<NEW_LINE;
        }
    }
    catch (sioexception ie) { status->abort(&ie); }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::splitread(const sfile &infile, sfile &outfile) {
    try {
        
        par->var_par.srv_par.min_sr[0] = 3;
        par->var_par.srv_par.min_sr[1] = 3;
        par->var_par.srv_par.min_sr[2] = 3;
        par->var_par.srv_par.min_sr[3] = 3;
        par->var_par.srv_par.min_sr[4] = 3;
        par->var_par.min_qual = 5.0;
        par->var_par.srv_par.max_fr_bias = 1.0;
        
        std::cout<<"Loading summary : '"<<infile<<"'"<<std::endl;
        summary->load(infile);
        outfile.make();
        std::cout<<"Exporting : '"<<outfile<<"'"<<std::endl;
        outfile<<"Type\tChr1\tPos1\tDir1\tChr2\tPos2\tDir2\tAlt\tRead(+)\tRead(-)\tQual"<<NEW_LINE;
        for (int i = 0; i < summary->ref_num; ++i) {
            for (int j = 0; j < 5; ++j) {
                auto &variants = summary->variants[5*i+j];
                if (variants.empty()) continue;
                for (auto it = variants.begin(); it < variants.end(); ++it) {
                    auto available = false;
                    if (par->target.size()) {
                        if ((it->type==SB_DELETION &&
                            par->target[it->pos1.idx].overlap(srange(it->pos1.pos+1, it->pos2.pos+1))) ||
                            (it->type==SB_DUPLICATION &&
                             par->target[it->pos1.idx].overlap(srange(it->pos2.pos+1, it->pos1.pos+1))) ||
                            (par->target[it->pos1.idx].include(it->pos1.pos+1) ||
                             par->target[it->pos2.idx].include(it->pos2.pos+1))
                            ) available = true;
                    }
                    else available = true;
                    if (available && par->ntarget.size()) {
                        if ((it->type==SB_DELETION &&
                            par->ntarget[it->pos1.idx].overlap(srange(it->pos1.pos+1, it->pos2.pos+1))) ||
                            (it->type==SB_DUPLICATION &&
                             par->ntarget[it->pos1.idx].overlap(srange(it->pos2.pos+1, it->pos1.pos+1))) ||
                            (par->ntarget[it->pos1.idx].include(it->pos1.pos+1) ||
                             par->ntarget[it->pos2.idx].include(it->pos2.pos+1))
                            ) available = false;
                    }
                    if (available &&
                        par->var_par.srv_par.min_sr[j] <= (it->pread+it->nread) &&
                        readBias(it->pread, it->nread) <= par->var_par.srv_par.max_fr_bias &&
                        par->var_par.min_qual <= phredVal(it->qual))
                        outfile<<varTypeStr(it->type)<<"\t"<<par->ref->name(it->pos1.idx)<<"\t"<<
                        it->pos1.pos+1<<"\t"<<(it->pos1.dir?"-":"+")<<"\t"<<
                        par->ref->name(it->pos2.idx)<<"\t"<<
                        it->pos2.pos+1<<"\t"<<(it->pos2.dir?"-":"+")<<"\t"<<
                        it->alt<<"\t"<<it->pread<<"\t"<<it->nread<<"\t"<<
                        snumber(phredVal(it->qual)).precised(2)<<NEW_LINE;
                    outfile.flush();
                }
            }
        }
    }
    catch (sioexception ie) { status->abort(&ie); }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::vsearch(int method, const sfile &infile, sfile &outfile) {
    try {
        std::cout<<"Loading bsm : '"<<infile<<"'"<<std::endl;
        summary->load(infile);
        search->search(method);
        if(par->db) {
            std::cout<<"Gene annotation."<<std::endl;
            variants->annotate(par->db, nullptr);
        }
        std::cout<<"Exporting variant list : '"<<outfile<<"'"<<std::endl;
        variants->name() = outfile.filename(false);
        if (par->outformat == "xls")
            variants->exportTable(outfile);
        else variants->save(outfile);
    }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::analyze(int method, const sfile &infile, sfile &outfile) {
    try {
        sfile bsm = infile+".bsm";
        summarize(infile, bsm);
        search->search(method);
        std::cout<<"Gene annotating."<<std::endl;
        for (auto it = variants->begin(); it < variants->end(); ++it) {
            if(-1 < it->variant.pos1.idx) it->variant.pos1.setName(par->ref->name(it->variant.pos1.idx));
            if(-1 < it->variant.pos2.idx) it->variant.pos2.setName(par->ref->name(it->variant.pos2.idx));
        }
        if(par->db) variants->annotate(par->db, &threads);
        sfile vout = outfile.parent();
        vout.downstair("analyzed");
        if (!vout.exist()) vout.make();
        vout.downstair(infile.filename()+"."+par->outformat);
        if (outfile.isDir()) outfile.downstair(infile.filename(false)+"."+par->outformat);
        std::cout<<"Exporting variant list : "<<outfile<<std::endl;
        variants->name() = infile.filename(false);
        if (par->outformat == "xls")
            variants->exportTable(outfile);
        else variants->save(outfile);
    }
    catch (sioexception ie) { status->abort(&ie); }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::vlmerge(const slist<sfile> &infiles, sfile &outfile) {
    try {
        for (auto it = infiles.begin(); it < infiles.begin(); ++it) {
            
            
            
            
        }
        if (outfile.isDir()) outfile.downstair(variants->name()+"_merged."+par->outformat);
        variants->save(outfile);
    }
    catch (sbioexception be) { status->abort(&be); }
}
void analyzer::subtract(const sfile &infile, sfile &outfile) {
    try {
        std::cout<<"Loading variant list : "<<infile<<std::endl;
        variants->load(infile);
        if (par->ctrl) {
            variants->subtract(*par->ctrl);
            if (outfile.isDir()) outfile.downstair(infile.filename(false)+"_subtracted."+par->outformat);
            std::cout<<"Exporting list : "<<infile<<std::endl;
            variants->save(outfile);
        }
        else std::cout<<"No control list."<<std::endl;
    } catch (sbioexception be) {
        status->abort(&be);
    }
}
void analyzer::common(const slist<sfile> &infiles, sfile &outfile) {
    try {
        for (auto it = infiles.begin(); it < infiles.begin(); ++it) {
            
            
            
            
        }
        if (outfile.isDir()) outfile.downstair(variants->name()+"_common."+par->outformat);
        variants->save(outfile);
    }
    catch (sbioexception be) { status->abort(&be); }
}

void analyzer::setParam(stk_param_t *p) {
    par = p;
    threads.setSize(par->an_par.thread_count);
    summary = new sbamsummary(&par->bam_par, p->ref);
    if (par->bam_par.parallel) {
        bam.resize(par->ref->size());
        for (int r = 0; r < par->ref->size(); ++r) {
            bam[r] = new sbamfile();
        }
    }
    else {
        bam.resize(1);
        bam[0] = new sbamfile();
    }
    variants = new varlist();
    variants->setParam(&par->var_par);
    variants->setReference(par->ref);
}

void analyzer::init() {
    if (status) status->init();
    sforeach(bam) { if (selement) selement->init(); }
    if (summary) summary->reset();
    if (variants) variants->init();
    if (reader) reader->init();
    if (cnan) cnan->init();
    if (search) search->init();
}
