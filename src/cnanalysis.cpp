#include "analyzer.h"

inline size_t maxLength(SeqList& ref) {
	size_t sz = 0;
	sfor(ref) {
		if (sz < $_.length()) sz = $_.length();
	}
	return sz;
}
/**
* CNAnalysis implementation
*/
stk::CNAnalysis::CNAnalysis() {
	par = nullptr;
	status = nullptr;
	logger = nullptr;
	smpl = nullptr;
	ctrl = nullptr;
};
stk::CNAnalysis::CNAnalysis(Analyzer* an) : CNAnalysis() { setParam(&an->par); }
stk::CNAnalysis::~CNAnalysis() {}
void stk::CNAnalysis::setParam(stk::Param* p) {
	par = p;
	status = &par->status;
	logger = &par->logger;
}
void stk::CNAnalysis::setData(NGSData* data, NGSData* control) {
	smpl = data;
	ctrl = control;
}
/**
* Get depth info
*/
void stk::CNAnalysis::depth(smath::Vector<svecf> &values) {
	auto it = values.begin();
	srange range;
	//
	if (par->target.empty()) {
		sfori(par->reference) {
			$_.resize((par->reference[i].length() - 1) / par->depth_bin + 1);
			range = srange(0, par->depth_bin - 1);
			sforeach(val, $_) {
				if (par->reference[i].length() <= range.end) range.end = par->reference[i].length() - 1;
				val = count(i, range);
				range.shift(par->depth_bin);
			}
			$NEXT;
		}
	}
	else {
		sfori(par->target) {
			sforeach(target, par->target[i]) {
				$_.resize(target.length() / par->depth_bin + 1);
				range = srange(target.begin, target.begin + par->depth_bin - 1);
				sforeach(val, $_) {
					if (target.end < range.end) range.end = target.end;
					val = count(i, range);
					range.shift(par->depth_bin);
				}
				$NEXT;
			}
		}
	}
}

/**
* Get copy number
*/
void stk::CNAnalysis::copynum(smath::Vector<svecf>& values) {
	auto it = values.begin();
	srange range;
	//
	if (par->target.empty()) {
		sfori(par->reference) {
			$_.resize((par->reference[i].length() - 1) / par->depth_bin + 1);
			range = srange(0, par->depth_bin - 1);
			sforeach(val, $_) {
				if (par->reference[i].length() <= range.end) range.end = par->reference[i].length() - 1;
				val = count(i, range) / (float)smpl->summary.avedp /
					(ctrl ? (count(i, range, true) / (float)ctrl->summary.avedp) : 1.f);
				range.shift(par->depth_bin);
			}
			$NEXT;
		}
	}
	else {
		sfori(par->target) {
			sforeach(target, par->target[i]) {
				$_.resize(target.length() / par->depth_bin + 1);
				range = srange(target.begin, target.begin + par->depth_bin - 1);
				sforeach(val, $_) {
					if (target.end < range.end) range.end = target.end;
					val = count(i, range) / smpl->summary.avedp /
						(ctrl ? (count(i, range, true) / (float)ctrl->summary.avedp) : 1.f);
					range.shift(par->depth_bin);
				}
				$NEXT;
			}
		}
	}
}

float _count(int bin, svecf* depth, const srange& range) {
	auto rng = srange(range.begin / bin, range.end / bin);
	if (rng.begin == rng.end) return depth->at(rng.begin);
	else {
		double value = 0;
		auto dp = depth->data() + rng.begin;
		value += (double)(*dp) * ((rng.begin + 1) * bin - range.begin);
		++rng.begin; ++dp;
		while (rng.begin < rng.end) {
			value += (double)(*dp) * bin; ++rng.begin; ++dp;
		}
		value += (double)(*dp) * (range.end - rng.begin * bin + 1);
		return (float)(value / range.length(true));
	}
}
float stk::CNAnalysis::count(int idx, const srange &range, bool control) {
	if (!smpl) throw slib::Exception(slib::NULL_ERROR, "Null value error.", slib::String("@") << __func__ << " l. " << 119 << " in '" << "H:\\マイドライブ\\dev\\github\\Sutoku\\src\\cnanalysis.cpp" << "'." << slib::NL << nullErrorText("NGSData object"));
	if (control && !ctrl) throw NullException(nullErrorText("Control data object"));
	return _count(smpl->summary.bin, &(control ? ctrl : smpl)->depth[idx], range);
}
float stk::CNAnalysis::ncount(int idx, const srange& range, bool control) {
	if (!smpl) throw NullException(nullErrorText("NGSData object"));
	if (control && !ctrl) throw NullException(nullErrorText("Control data object"));
	return _count(smpl->summary.bin, &normalized[control?1:0][idx], range);
}
float stk::CNAnalysis::copy(int idx, const srange& range) {
	return _count(smpl->summary.bin, &copies[idx], range);
}
void _normalize(svecf* dp, svecf *ndp, sngs::Summary *sum, stk::Param *par) {
	if (dp) { 
		if (sum->avedp < snum::D_EPS) throw DivZeroException("Average depth was zero or too small.");
		sfor2(*dp, *ndp) $_2 = $_1 / sum->avedp; 
	}
	else sfor(*ndp) $_ = 1.f;
}
void _copy(svecf *ndp1, svecf *ndp2, svecf *cp, CN_METHOD method) {
	auto p1 = ndp1->data(),
		p2 = ndp2->data(),
		v = cp->data();
	auto n = ndp1->size();
	sforin(i, 0, n) { 
		if (*p2 < snum::F_EPS) *v = (*p1);
		else *v = (*p1) / (*p2); 
		++v; ++p1; ++p2; 
	}
	switch (method) {
	case CN_METHOD::RAW:
		return;
	default:
		break;
	}
}
void stk::CNAnalysis::analyze() {
	auto& num = smpl->summary.refnum;
	normalized[0].resize(num);
	normalized[1].resize(num);
	copies.resize(num);
	// Noramalize
	sforin(i, 0, num) {
		normalized[0][i].resize(smpl->depth[i].size());
		normalized[1][i].resize(smpl->depth[i].size());
		par->threads.addTask(_normalize, &smpl->depth[i], &normalized[0][i], &smpl->summary, par);
		//_normalize(&smpl->depth[i], &normalized[0][i], &smpl->summary, par);		
		if (ctrl) {
			par->threads.addTask(_normalize, &ctrl->depth[i], &normalized[1][i], &ctrl->summary, par);
			//_normalize(&ctrl->depth[i], &normalized[1][i], &ctrl->summary, par);
		}
		else _normalize(nullptr, &normalized[1][i], nullptr, nullptr);
	}
	par->threads.complete();
	// Copynumber calc.
	sforin(i, 0, num) {
		copies[i].resize(smpl->depth[i].size());
		par->threads.addTask(_copy, &normalized[0][i], &normalized[1][i], &copies[i], par->varp.cnvp.method);
		//_copy(&normalized[0][i], &normalized[1][i], &copies[i], par->varp.cnvp.method);
	}
	par->threads.complete();
}
void stk::CNAnalysis::analyze(CNData& cnd, int r, srange range) {
	// Correction for AT/GC ratio etc...
	float correct = 1.f;

	cnd.depth[0] = count(r, range, smpl);
	cnd.ndepth[0] = cnd.depth[0] / (float)smpl->summary.avedp * correct;
	if (ctrl) {
		cnd.depth[1] = count(r, range, ctrl);
		cnd.ndepth[1] = cnd.depth[1] / (float)ctrl->summary.avedp * correct;
	}
	else {
		cnd.depth[1] = 1.f;
		cnd.ndepth[1] = 1.f;
	}
	// Copy number value
	switch (par->varp.cnvp.method) {
	case CN_METHOD::RAW:
		cnd.copy = cnd.ndepth[0] / cnd.ndepth[1];
		break;
	default:
		break;
	}
}
void stk::CNAnalysis::reset() {
	sfor(normalized[0]) $_.reset(0.f);
	sfor(normalized[1]) $_.reset(0.f);
	sfor(copies) $_.reset(0.f);
}
