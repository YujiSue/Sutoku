#include "../src/analyzer.h"
///////////////////////////////////////////// slib namespace ////////////////////////////////////////////
using namespace slib;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sutil;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" {
    splugin vfilter(VarList& vl, const stk::Param& par, const SDictionary& pref) {
        try {
            sforeach(var, vl) {
                if (var->flag == sbio::SMALL_VARIANT) {
                    if (var->type == sbio::SNV || var->type == sbio::MNV) {
                        var->attribute["filter"] = "PASS";
                        if(var->attribute["info"]["QD"].doubleValue() < 2.0) {
                            var->flag = UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Low Qual/Depth";
                        }
                        if(var->attribute["info"]["FS"].doubleValue() > 60.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Biased strand (Fisher)";
                        }
                        if(var->attribute["info"]["SOR"].doubleValue() > 3.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Biased strand (Odds)";
                        }
                        if(var->attribute["info"]["MQ"].doubleValue() < 40.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Low Mapping Qual";
                        }
                        if(var->attribute["info"]["MQRankSum"].doubleValue() < -12.5) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Alt Low MQ ";
                        }
                        if(var->attribute["info"]["ReadPosRankSum"].doubleValue() < -8.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Alt Biased Pos";
                        }
                    }
                    else {
                        var->attribute["filter"] = "PASS";
                        if(var->attribute["info"]["QD"].doubleValue() < 2.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Low Qual/Depth";
                        }
                        if(var->attribute["info"]["FS"].doubleValue() > 200.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Biased strand (Fisher)";
                        }
                        if(var->attribute["info"]["ReadPosRankSum"].doubleValue() < -20.0) {
                            var->flag = sbio::UNAVAILABLE_FLAG;
                            var->attribute["filter"] = "Alt Biased Pos";
                        }
                    }
                }
            }
            //vl.tidyUp();
            /////////////////////////////////////////////////
            return 0;
        }
        catch (Exception ex) {
            ex.print();
            return ex.code;
        }
    }
}