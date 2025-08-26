#include "../src/analyzer.h"
///////////////////////////////////////////// slib namespace ////////////////////////////////////////////
using namespace slib;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sutil;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" {
    enum MUTAGEN {
        ENU = 1,
        EMS = 2,
    };

    splugin vfilter(VarList& vl, const stk::Param& par, const SDictionary& pref) {
        try {
            /////////////////////////////////////////////////

            int bin = pref.hasKey("bin") ? pref["bin"] : (par.reference.total() / (560 * snum::PI));
            int b = 0;
            while (0 < bin) {
                bin >>= 1; ++b;
            }
            bin = (1 << b);

            int bin2 = bin / 2;
            int mincount = pref.hasKey("mcount") ? pref["mcount"] : 5;

            int mg = 0;
            if (pref.hasKey("mutagen")) {
                if (pref["mutagen"] == "ENU") mg = (int)MUTAGEN::ENU;
                else if (pref["mutagen"] == "EMS") mg = (int)MUTAGEN::EMS;
            }

            String opath;
            if (pref.hasKey("output")) {}
            else {
                auto &vlpath = vl.path();
                opath = sfs::joinPath(sfs::splitPath(vlpath).first, sfs::fileName(vlpath, false) + ".mutagen_homo.filtered.html");
            }
            //
            slib::SFile f(opath, slib::sio::MAKE);

            //
            f << "<!DOCTYPE html><html><head>\n<meta charset=\"utf-8\">\n";
            f << "<script src=\"https://d3js.org/d3.v4.min.js\"></script>\n<script src=\"https://unpkg.com/circos@2.0.2/dist/circos.min.js\"></script>\n";
            f << "</head>\n";
            f << "<body><div id=\"chart\"></div>\n";
            f << "<script>\n";
            f << "const circos = new Circos({container:'#chart', width:640, height:640 });\n";
            f << "circos.layout([\n";
            sfor(par.reference) {
                f << "{id:'" << $_.name << "', len:" << $_.length() << ", label:'" << $_.name << "'},\n";
            }
            f << "], {innerRadius: 275, outerRadius: 285, labels: {display: true, position : 'center', size : 12, color : '#000', radialOffset : 20 }});\n";
            //

            intarray2d counts[2];
            srange hrange;
            bool match;
            //
            if (par.vcontrol.size()) {
                counts[0].resize(par.reference.size());
                sfor2(counts[0], par.reference) {
                    $_1.resize((($_2.length() - 1) / bin2) + 1, 0);
                }
                sforeach(var, par.vcontrol) {
                    if (!(var->genotype == sbio::HOMO_VAR)) continue;
                    match = false;
                    if (var->type == sbio::SNV) {
                        if (mg == (int)MUTAGEN::EMS && (
                            (var->attribute["_ref_"] == "G" && var->alt == "A") ||
                            (var->attribute["_ref_"] == "C" && var->alt == "T"))) match = true;
                        else if (mg == (int)MUTAGEN::ENU && (
                            (var->attribute["_ref_"] == "A" && (var->alt == "T" || var->alt == "G")) ||
                            (var->attribute["_ref_"] == "T" && (var->alt == "A" || var->alt == "C")) ||
                            (var->attribute["_ref_"] == "G" && var->alt == "A") ||
                            (var->attribute["_ref_"] == "C" && var->alt == "T"))) match = true;

                        if (match) {
                            hrange = srange((var->pos[0].begin - 1) / bin2, (var->pos[0].end - 1) / bin2 + 1);
                            if (0 < hrange.begin) --hrange.begin;
                            sforin(h, hrange.begin, hrange.end) ++counts[0][var->pos[0].idx][h];
                        }
                    }
                }
            }

            //
            f << "circos.stack('variants', [\n";
            counts[1].resize(par.reference.size());
            sfor2(counts[1], par.reference) {
                $_1.resize(($_2.length() - 1) / bin2, 0);
            }
            sforeach(var, vl) {
                //
                f << "{block_id: '" << par.reference[var->pos[0].idx].name << "', start: " << (var->pos[0].begin - 1) << ", end: " << min((var->pos[0].end - 1), par.reference[var->pos[0].idx].length()) << "},\n";
                //
                if (var->genotype != sbio::HOMO_VAR) {
                    var->attribute["filter"] = "Not Homozygous";
                    var->flag = UNAVAILABLE_FLAG;
                    continue;
                }
                match = false;
                if (var->type == sbio::SNV) {
                    if (mg == (int)MUTAGEN::EMS && (
                        (var->attribute["_ref_"] == "G" && var->alt == "A") ||
                        (var->attribute["_ref_"] == "C" && var->alt == "T"))) match = true;
                    else if (mg == (int)MUTAGEN::ENU && (
                        (var->attribute["_ref_"] == "A" && (var->alt == "T" || var->alt == "G")) ||
                        (var->attribute["_ref_"] == "T" && (var->alt == "A" || var->alt == "C")) ||
                        (var->attribute["_ref_"] == "G" && var->alt == "A") ||
                        (var->attribute["_ref_"] == "C" && var->alt == "T"))) match = true;
                    if (match) {
                        var->attribute["filter"] += S(",") + (mg == (int)MUTAGEN::EMS ? "EMS" : (mg == (int)MUTAGEN::ENU ? "ENU" : "")) + "-type";
                        hrange = srange((var->pos[0].begin - 1) / bin2, (var->pos[0].end - 1) / bin2 + 1);
                        if (0 < hrange.begin) --hrange.begin;
                        sforin(h, hrange.begin, hrange.end) ++counts[1][var->pos[0].idx][h];
                    }
                }
            }
            f << "], {innerRadius: 150, outerRadius : 180, color: 'red', tooltipContent: function(d){return `${d.block_id}:${d.start} `}});\n";
            
            //
            f << "circos.histogram('hist', [\n";
            sforin(r, 0, counts[1].size()) {
                sfori(counts[1][r]) {
                    f << "{block_id: '" << par.reference[r].name << "', start: " << (i * bin2) << ",end : " << min((i * bin2 + bin - 1), par.reference[r].length()) << ", value: " << counts[1][r][i] << "},\n";
                }
            }
            f << "], {innerRadius: 180, outerRadius : 250, color : 'steelblue',tooltipContent: function(d){return `${d.block_id}:${d.start+1} - ${d.end+1}`}});\n";

            //
            if (par.vcontrol.size()) {
                f << "circos.histogram('control',[\n";
                sforin(r, 0, counts[0].size()) {
                    sfori(counts[0][r]) {
                        f << "{block_id: '" << par.reference[r].name << "', start: " << (i * bin2) << ",end : " << min((i * bin2 + bin - 1), par.reference[r].length()) << ", value: " << counts[0][r][i] << "},\n";
                    }
                }
                f << "], {innerRadius: 180, outerRadius : 250, color : 'tomato', opacity : 0.6,tooltipContent: function(d){return `${d.block_id}:${d.start+1} - ${d.end+1}`}});\n";
            }

            //
            f << "circos.render();\n";
            f << "</script></body></html>";

            /////////////////////////////////////////////////
            return 0;
        }
        catch (Exception ex) {
            ex.print();
            return ex.code;
        }
    }
}