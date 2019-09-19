#include "analyzer.hpp"

using namespace slib;
using namespace slib::sio;
using namespace slib::sdb;
using namespace slib::sbio;
using namespace slib::sngs;

using namespace stk;

extern inline sstring usage() {
    sstring str;
    str<<"usage:"<<NEW_LINE<<"  Sutoku <mode> \\\\"<<NEW_LINE<<
    "    --BAM files(.bam)/directory\\\\"<<NEW_LINE<<
    "    --BSM files(.bsm)/directory\\\\"<<NEW_LINE<<
    "    --reference file\\\\"<<NEW_LINE<<
    "    --annotation file\\\\"<<NEW_LINE<<
    "    --control file(.bsm)\\\\"<<NEW_LINE<<
    "    --target file(.txt)\\\\"<<NEW_LINE<<
    "    --method {CN or SR or CNSR}(+..)(+..)\\\\"<<NEW_LINE<<
    "    --param file(.json/plist)\\\\"<<NEW_LINE<<
    "    --output directory\\\\"<<NEW_LINE<<
    "    --otype {txt or json or xml(plist) or xls(TSV)}\\\\"<<NEW_LINE<<
    "    --bin integer\\\\"<<NEW_LINE<<
    "    --thread integer"<<NEW_LINE<<NEW_LINE<<
    "mode:"<<NEW_LINE<<
    "    analyze"<<NEW_LINE<<
    "        Perform \"summary\" and \"vsearch\" sequentially."<<NEW_LINE<<
    "    copynum"<<NEW_LINE<<
    "        Make copy number array from the input and control BSM files."<<NEW_LINE<<
    "        *Require --BSM option."<<NEW_LINE<<
    "    readinfo"<<NEW_LINE<<
    "        Show inforamtion of the reads in the target regions."<<NEW_LINE<<
    "        *Require --BAM option."<<NEW_LINE<<
    "    summary"<<NEW_LINE<<
    "        Make BSM files from input BAM data."<<NEW_LINE<<
    "        *Require --BAM option."<<NEW_LINE<<
    "    srlist"<<NEW_LINE<<
    "        Make split read list from input BSM data."<<NEW_LINE<<
    "        *Require --BSM option."<<NEW_LINE<<
    "    vsearch"<<NEW_LINE<<
    "        Search variants based on the input BSM files."<<NEW_LINE<<
    "        *Require --BSM option.";
    
    
    
    
    
    return str;
}

//readview, copyinfo, splitinfo, summary, vsearch, analyze, integrate, merge, subract, common

int main(int argc, const char * argv[]) {
    short mode = 0;
    stk_param_t par;
    sdate date;
    slist<sfile> inputs;
    sfile output;
    analyzer *an;
    
    if(argc < 2) {
        std::cout<<usage()<<std::endl;
        return 1;
    }
    sstring option(argv[1]);
    sdictionary params;
    sforin(i, 2, argc) {
        if (argv[i][0] == '-') {
            if (i == argc-1 || argv[i+1][0] == '-')
                params[&argv[i][1]] = true;
            else params[&argv[i][1]] = argv[i+1];
            ++i;
        }
        else params["_args"].add(argv[i]);
    }
    if (option == "readview") mode = READ_INFO_VIEW;
    else if (option == "summary") mode = MAKE_SUMMARY;
    else if (option == "analyze") mode = MAKE_SUMMARY|VAR_SEARCH;
    else if (option == "integrate") mode = SUMMARY_INTEGRATE;
    else if (option == "cpinfo") mode = COPY_NUM_VIEW;
    else if (option == "srinfo") mode = SPLIT_READ_VIEW;
    else if (option == "vsearch") mode = VAR_SEARCH;
    else if (option == "merge") mode = MERGE_VLIST;
    else if (option == "subtract") mode = SUBTRACT_VLIST;
    else if (option == "common") mode = COMMON_VLIST;
    else if (option == "template") mode = TEMPLATE_EXPORT;
    else {
        std::cout<<"'"<<option<<"' command is not defined."<<std::endl;
        return 1;
    }
    sfile refpath, dbpath, parpath, bgpath, ctrlpath, trgtpath, ntrgtpath;
    try {
        std::cout<<params<<std::endl;
        
        if (params["param"]) par.load(params["param"]);
        if (params["reference"]) par.setRefPath(params["reference"]);
        if (params["annotation"]) par.setDBPath(params["annotation"]);
        if (params["control"]) par.setBGPath(params["control"]);
        if (params["vcontrol"]) par.setCtrlPath(params["vcontrol"]);
        if (params["target"]) par.setTgtPath(params["target"]);
        if (params["ntarget"]) par.setNtgtPath(params["ntarget"]);
        if (params["bam"]) {
            auto list = params["bam"].split(",");
            sforeach(list) {
                sfile file(selement);
                if(file.isDir()) inputs.append(file.filelist({ "bam" }));
                else if(file.extension() == "bam") inputs.add(file);
            }
        }
        if (params["bsm"]) {
            auto list = params["bsm"].split(",");
            sforeach(list) {
                sfile file(selement);
                if(file.isDir()) inputs.append(file.filelist({ "bsm" }));
                else if(file.extension() == "bsm") inputs.add(file);
            }
        }
        if (params["vlist"]) {
            auto list = params["vlist"].split(",");
            sforeach(list) {
                sfile file(selement);
                if(file.isDir()) inputs.append(file.filelist({ "plist", "json", "vcf", "soml" }));
                else inputs.add(file);
            }
        }
        if (params["method"]) {
            sstringarray arr = params["method"].split("+");
            sforeach(arr) {
                if (selement == "CN") mode |= CNV_DETECTION;
                else if(selement == "SR") mode |= SRV_DETECTION;
                else if(selement == "CNSR") mode |= CNSRV_DETECTION;
            }
        }
        if (params["outdir"]) output = params["outdir"];
        if (params["otype"]) par.outformat = params["otype"];
        if (params["bin"]) par.an_par.cp_bin = params["bin"];
        if (inputs.empty()) {
            std::cout<<"No input."<<std::endl;
            return 1;
        }
        std::cout<<date.formed("YYYY/MM/DD HH:mm:ss")+" started."<<NEW_LINE;
        an = new analyzer(&par);
        if (!an) {
            std::cout<<"Analyzer initialization failed."<<std::endl;
            return 1;
        }
        
        sforeach(inputs) {
            sfile out;
            if (output.empty()) out = it->parent();
            else out = output;
            if (!out.isDir()) out.upstair();
            if (!(it->exist())) {
                std::cout<<"File was not found. : "<<*it<<std::endl;
                continue;
            }
            if (mode == READ_INFO_VIEW) {
                if (lower((*it).extension()) != "bam") continue;
                out = out.child("read_info.txt");
                out.make();
                an->readview(*it, out);
            }
            else if (mode&MAKE_SUMMARY) {
                if (lower((*it).extension()) != "bam") continue;
                out = out.child(it->filename()+".bsm");
                if (mode&VAR_SEARCH)
                    an->analyze(mode&0xF000, *it, out);
                else an->summarize(*it, out);
            }
            else if (mode == SUMMARY_INTEGRATE) {
                out = out.child("integrated.bsm");
                an->sumintegrate(inputs, out);
                break;
            }
            else if (mode == COPY_NUM_VIEW) {
                if (lower((*it).extension()) != "bsm") continue;
                out = out.child(it->filename(false)+"_cp.txt");
                an->copynum(*it, out);
            }
            else if (mode == SPLIT_READ_VIEW) {
                if (lower((*it).extension()) != "bsm") continue;
                out = out.child(it->filename(false)+"_sr.txt");
                an->splitread(*it, out);
            }
            else if ((mode&0xFF) == VAR_SEARCH) {
                if (lower((*it).extension()) != "bsm") continue;
                out = out.child(it->filename(false)+"_variant."+par.outformat);
                an->vsearch(mode&0xF000, *it, out);
            }
            else if (mode == MERGE_VLIST) {
                an->vlmerge(inputs, out);
                break;
            }
            else if (mode == SUBTRACT_VLIST) an->subtract(*it, out);
            else if (mode == COMMON_VLIST) {
                an->common(inputs, out);
                break;
            }
            else if (mode == TEMPLATE_EXPORT) {
                
                
            }
            if(an->status->err)
                std::cout<<sdate().formed("YYYY/MM/DD HH:mm:ss  ")<<
                "[Error]: "<<an->status->err<<" @"<<an->status->src<<NEW_LINE<<an->status->msg<<std::endl;
            else std::cout<<"Completed."<<std::endl;
            an->init();
        }
        
    }
    catch (sioexception ie) { ie.print(); }
    catch (sbioexception be) { be.print(); }
    if (an) delete an;
    return 0;
}
