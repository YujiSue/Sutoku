#include "param.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sbio;
using namespace stk;

inline void loadTarget(const char* path, Array<sregion> &regions, SeqList &ref) {
	// Open
	slib::sbio::BEDFile target(path);
	// Ref. name => index
	target.setRef(ref);
	target >> regions;
	// Set zero-based position
	sfor(target) $_.shift(-1);
}
stk::Param::Param() {
	async_load = false;
	detect_sv = false;
	ignore_dp = false;
	clipped = false;
	depth_bin = 1;
	min_clip = 20;
	clip_buffer = (1 << 10);
	min_err_prob = 1E-6;

	sforin(i, 0, 10) detect[i] = true;
	
	mlcheck = false;

	annotation = false;
	var_site = 0xFFFF;

	max_thread = 4;
	
	seqtype = sngs::SEQ_TYPE::SINGLE;
	seqp.setRefType(DNA_SEQ4);
	seqp.setSeed((((int)min_clip / 4) - 1) * 4);
	seqp.min_match = (int)min_clip * 3 / 4;
	seqp.max_gap = 2;
	seqp.max_miss = 2;
	seqp.ext_threshold = 0.75;
}
stk::Param::Param(slib::SDictionary& pref) : stk::Param() { set(pref); }
stk::Param::~Param() {}
void stk::Param::set(SDictionary &pref) {
	command = pref["_cmd_"];
	try {
		if (pref["param"]) load(pref["param"]);
		if (pref["bam"]) {
			if (command == "readinfo" || command == "summary" || command == "analyze") {
				if (sfs::isDir(pref["bam"])) 
					inputs.append(SDirectory(pref["bam"]).fileList({ "bam" }));
				else {
					auto bams = pref["bam"].split(",");
					sforeach(in, bams) { if (sfs::exist(in)) inputs.add(in); }
				}
			}
		}
		if (pref["bsm"]) {
			if (command == "depth" || command == "copynum" || command == "splitread" || 
				command == "vsearch" || command == "integrate" || command == "subtract") {
				if (sfs::isDir(pref["bsm"])) 
					inputs.append(SDirectory(pref["bsm"]).fileList({ "bsm" }));
				else {
					auto bsms = pref["bsm"].split(",");
					sforeach(in, bsms) { if (sfs::exist(in)) inputs.add(in); }
				}
			}
		}
		if (pref["vlist"]) {
			if (command == "merge" || command == "unique" || command == "common") {
				if (sfs::isDir(pref["vlist"])) 
					inputs.append(SDirectory(pref["vlist"]).fileList({ "txt", "tsv", "json", "vcf" }));
				else {
					auto vls = pref["vlist"].split(",");
					sforeach(in, vls) { if (sfs::exist(in)) inputs.add(in); }
				}
			}
		}
		if (pref["reference"]) {
			reference.load(pref["reference"]);
			logger.log(S("Load reference data.") << LF <<
				" > Sequence count: " << reference.size() << LF <<
				" > Total size: " << reference.total() << LF <<
				" > Species: " << (reference.attribute.hasKey("species") ? reference.attribute["species"] : "") << LF <<
				" > Version: " << (reference.attribute.hasKey("version") ? reference.attribute["version"] : ""));
			mlock.resize(reference.size());
		}
		if (pref["annotdb"]) {
			annotation = true;
			annotdb.open(pref["annotdb"]);
			annotdb.loadGenes();
			logger.log(S("Load annotation dataset.") << LF <<
				" > Version: " << (reference.attribute.hasKey("version") ? reference.attribute["version"] : ""));
		}
		if (pref["control"]) {
			logger.log("Load control data.");
			control.load(pref["control"]);
		}
		if (pref["vcontrol"]) {
			logger.log("Load control variants.");
			vcontrol.load(pref["vcontrol"], &reference);
		}
		if (pref["target"]) {
			logger.log("Load target regions.");
			loadTarget(pref["target"], target, reference);
		}
		if (pref["outdir"]) {
			if (!sfs::exist(pref["outdir"])) sfs::makeDir(pref["outdir"]);
			outdir = pref["outdir"];
			if (!sfs::exist(outdir)) sfs::makeDir(outdir);
		}
		if (pref["oformat"]) oformat = pref["oformat"];
		if (pref["paired"]) seqtype = sngs::SEQ_TYPE::PAIRED;
		if (pref["async-load"]) async_load = true;
		if (pref["detect-sv"]) detect_sv = true;
		if (pref["ignore-pcrdp"]) ignore_dp = true;
		if (pref["clip-only"]) clipped = true;
		if (pref["depth-bin"]) depth_bin = pref["depth-bin"];
		if (pref["cliplen"]) min_clip = pref["cliplen"];
		if (pref["realign-param"]) seqp.set(pref["realign-param"]);
		if (pref["variant-param"]) varp.set(pref["variant-param"]);
		//
		if (pref["thread"]) max_thread = pref["thread"];
		threads.setSize(max_thread);
		//
		if (pref["log-dir"]) _logpath = sfs::joinPath(pref["log-dir"], "sutoku.log");
		else if (sfs::exist(outdir)) _logpath = sfs::joinPath(outdir, "sutoku.log");
		else if (inputs.size()) _logpath = sfs::joinPath(sfs::splitPath(inputs[0]).first, "sutoku.log");
		else _logpath = sfs::joinPath(ssys::home(), "sutoku.log");
		logger.log(S("The log file path is set to '") << _logpath << "'");
		logger.open(_logpath);
	}
	catch (Exception ex) {
		logger.log(ex);
	}
}
void stk::Param::load(const char* path) { 
	auto pars = sjson::load(path);
	seqtype = pars["paired"] ? sngs::SEQ_TYPE::PAIRED : sngs::SEQ_TYPE::SINGLE;
	async_load = pars["async-load"] ? true : false;
	detect_sv = pars["detect-sv"] ? true : false;
	ignore_dp = pars["ignore-pcrdp"] ? true : false;
	annotation = pars["annotation"] ? true : false;
	if (pars["depth-bin"]) depth_bin = pars["depth-bin"];
	if (pars["cliplen"]) min_clip = pars["cliplen"];
	if (pars["clipnum"]) clip_buffer = pars["clipnum"];
	if (pars["minerr"]) min_err_prob = pars["minerr"];
	if (pars["vdetect"]) {
		sforin(i, 0, 10) detect[i] = pars["vdetect"][i];
	}
	if (pars["vsite"]) var_site = pars["vsite"];
	if (pars["thread"]) max_thread = pars["thread"];
	seqp.set(pars["realign-param"]);
	varp.set(pars["variant-param"]);
}
void stk::Param::save(const char* path) { sjson::save(toObj(), path); }
sobj stk::Param::toObj() {
	SArray detect_filters;
	sforin(i, 0, 10) detect_filters.add(detect[i]);
	return {
		D_("reference", sobj({
			D_("path", reference.path()),
			D_("version", reference.attribute["version"]),
			D_("species", reference.attribute["species"])
			})),
		D_("annotdb", sobj({
			D_("path", annotdb.path)
			})),
		D_("control", control.path()),
		D_("vcontrol", vcontrol.path()),

		D_("paired", seqtype == sngs::SEQ_TYPE::PAIRED),
		D_("async-load", async_load),
		D_("detect-sv", detect_sv),
		D_("ignore-pcrdp", ignore_dp),
		D_("clip-only", clipped),

		D_("depth-bin", depth_bin),
		D_("cliplen", min_clip),
		D_("clipnum", clip_buffer),
		D_("minerr", min_err_prob),

		D_("vdetect", detect_filters),

		D_("annotation", annotation),
		D_("vsite", var_site),

		D_("thread", max_thread),

		D_("realign-param", seqp.toObj()),
		D_("variant-param", varp.toObj())
	};
}