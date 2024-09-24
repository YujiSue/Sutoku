#include "analyzer.h"

inline bool check1(SVar& var, VarParam* vpar) {
	if (var.total() < vpar->svp.min_read || // read count
		vpar->svp.read_bias < var.bias() || // fwd./rev. bias
		sbio::sutil::phredVal(var.qual) < vpar->svp.min_qual) { // qual.
		return false;
	}
	return true;
}
inline bool check2(SVar& var1, SVar& var2, VarParam* vpar) {
	if (var1.total() < vpar->svp.min_read ||
		var2.total() < vpar->svp.min_read ||
		var1.total() + var2.total() < vpar->svp.comp_min_read || // read count
		vpar->svp.read_bias < var1.bias() || // fwd./rev. bias
		vpar->svp.read_bias < var2.bias() || // fwd./rev. bias
		vpar->svp.comb_bias < sbio::sutil::combBias(var1.read, var2.read) || // left/right bias
		sbio::sutil::phredVal(1.0 - ((1.0 - var1.qual) * (1.0 - var2.qual))) < vpar->svp.min_qual) { // qual.
		return false;
	}
	return true;
}
//
void stk::selectDelCandidate(int t, Array<SVar*>* candidate, Array<SVar>* variants, stk::Param* par) {
	if (variants->empty()) return;
	if (par->detect[0]) {
		sfor(*variants) {
			if (check1($_, &par->varp) && 
				par->varp.svp.size_range[0].include($_.pos[1].begin - $_.pos[0].end - 1)) {
				candidate->add($.ptr());
			}
			++par->status.current_task[t];
		}
	}
	else par->status.current_task[t] += variants->size();
}
void stk::getDelCpy(Array<SPointer<sbio::Variant>> *primary, Array<SVar*>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp, bg, freq;
	sushort cnv;
	srange delrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		delrange.begin = $_->pos[0].end + 1;
		delrange.end = $_->pos[1].begin - 1;
		
		// Freq. check
		freq = (2.f * (float)$_->total()) /
			(cna->count($_->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		
		// Copy value check
		cp = cna->copy($_->pos[0].idx, delrange);
		cnv = par->varp.cnvp.evaluate(cp);
		if ((cnv & 0xFF) != DELETION) continue;
		
		// Add to primary variant list
		Variant* var = new Variant(*$_);
		var->flag = SR_VARIANT;
		var->pos[0] = sbpos($_->pos[0].idx, delrange.begin + 1, delrange.end + 1, false);
		var->genotype = (cnv >> 8) & 0xFF;
		var->frequency = freq;
		var->copy[0] = cp;
		var->depth[0][0] = cna->count($_->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		//
		primary->add(var);
	}
}
//
void stk::selectDupCandidate(int t, Array<SVar*>* candidate, Array<SVar>* variants, stk::Param* par) {
	if (variants->empty()) return; 
	if (par->detect[1]) {
		sfor(*variants) {
			if (check1($_, &par->varp) && 
				par->varp.svp.size_range[1].include($_.pos[0].end - $_.pos[1].begin + 1)) {
				candidate->add($.ptr());
			}
			++par->status.current_task[t];
		}
	}
	else par->status.current_task[t] += variants->size();
}
void stk::getDupCpy(Array<SPointer<sbio::Variant>>* primary, Array<SVar*>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp, bg, freq;
	sushort cnv;
	srange duprange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		duprange.begin = $_->pos[1].begin;
		duprange.end = $_->pos[0].end;
		// Freq. check
		if (duprange.length(true) < par->varp.svp.break_site_len) {
			freq = (float)$_->total() /
				cna->count($_->pos[0].idx, srange(duprange.begin, duprange.end));
		}
		else {
			freq = (2.f * (float)$_->total()) /
				(cna->count($_->pos[0].idx, srange(duprange.begin, duprange.begin  + par->varp.svp.break_site_len)) +
					cna->count($_->pos[0].idx, srange(duprange.end - par->varp.svp.break_site_len + 1, duprange.end)));
		}
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_->pos[0].idx, duprange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp = cna->copy($_->pos[0].idx, duprange);
		cnv = par->varp.cnvp.evaluate(cp);
		if ((cnv & 0xFF) != DUPLICATION && (cnv & 0xFF) != MULTIPLICATION) continue;
		// Add to primary variant list
		Variant* var = new Variant(*$_);
		var->flag = SR_VARIANT;
		var->pos[0] = sbpos($_->pos[0].idx, duprange.begin + 1, duprange.end + 1, false);
		var->genotype = (cnv >> 8) & 0xFF;
		var->frequency = freq;
		var->copy[0] = cp;
		var->depth[0][0] = cna->count($_->pos[0].idx, duprange, false);
		var->depth[0][1] = bg;
		//
		primary->add(var);
	}
}
//
void stk::selectInsCandidate(int t, Array<SVar*>* candidate, Array<SVar>* variants, stk::Param* par) {
	if (variants->empty()) return; 
	if (par->detect[1]) {
		sfor(*variants) {
			if (check1($_, &par->varp) && 
				par->varp.svp.size_range[2].include((int)$_.alt.size())) {
				if (par->target.size()) {
					if (par->target[$_.pos[0].idx].overlap(srange($_.pos[0].end, $_.pos[1].begin))) candidate->add($.ptr());
				}
				else candidate->add($.ptr());
			}
			++par->status.current_task[t];
		}
	}
	else par->status.current_task[t] += variants->size();
}
void stk::getInsCpy(Array<SPointer<sbio::Variant>>* primary, Array<SVar*>* candidates, CNAnalysis* cna, stk::Param* par) {}

// Insertion from downstream : DEL(left) + DUP(right)
//    lp0 --->|*  rp1 *:--->
//                  lp1 |+++> ...  rp0 +++>:   
//
// Insertion from upstream : DUP(left) + DEL(right)
//                    lp0 --->|    rp1 :--->
//  lp1 *|+++> ...  rp0 +++>:*   
void stk::selectComplexCandidate(int t, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants1, Array<SVar>* variants2, stk::Param* par) {
	if (variants1->empty()) return;
	if (par->detect[3]) {
		auto left = variants1->begin(), right = variants2->begin(), lend = variants1->end(), rend = variants2->end();
		while (left < lend) {
			right = variants2->begin();
			while (right < rend) {
				if (right->pos[1].begin < left->pos[0].end - par->varp.svp.break_site_len ||
					right->pos[0].end < left->pos[0].begin) {
					++right; continue;
				}
				if (check2(*left, *right, &par->varp) &&
					left->pos[1].begin <= right->pos[0].begin &&
					left->pos[1].end <= right->pos[0].end &&
					left->pos[0].end < right->pos[1].begin) {
					candidate->add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
				}
				++right;
			}
			++par->status.current_task[t];
			++left;
		}
	}
	else par->status.current_task[t] += variants1->size();
}
void stk::getComplexCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.first->pos[0].end + 1 < $_.second->pos[1].begin) {
			delrange.begin = $_.first->pos[0].end + 1;
			delrange.end = $_.second->pos[1].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.first->pos[0].end;
			delrange.end = $_.second->pos[1].begin;
		}
		//
		insrange.begin = $_.first->pos[1].begin;
		insrange.end = $_.second->pos[0].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.first->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[0].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[0].idx, insrange.begin + 1, insrange.end + 1);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[0].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[0].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}

// Simple inversion
//   lp0 --->|<--- lp1
//                rp0 <---|---> rp1
// 
// Inverted insertion from upstream
// rp0 <---|---> rp1 <complement> rp1' <---|---> rp0' 
//   lp0 --->| :---> rp0'
//                   lp1 |<+++    rp1' <+++:
//
// Inverted insertion from downstream
// lp0 --->|<--- lp1 <complement> lp1' --->|<--- lp0' 
//                   lp1' --->| :---> rp1
//  lp0' |<+++    rp0 <+++:
void stk::selectInvCandidate(int t, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants, stk::Param* par) {
	if (variants->empty()) return;
	auto left = variants->begin(), rbeg = variants->begin(), right = variants->begin(), lend = variants->end(), rend = variants->end();
	while (rbeg < rend && !rbeg->pos[0].dir) { ++rbeg; }
	lend = rbeg;
	if (lend == rend || lend == left) {
		par->status.current_task[t] = variants->size();
		return;
	}
	else par->status.current_task[t] = rend - rbeg;
	while (left < lend) {
		right = rbeg;
		while (right < rend) {
			if (check2(*left, *right, &par->varp)) {
				// Simple inversion
				if (par->detect[4] &&
					right->pos[0].begin <= left->pos[1].begin &&
					right->pos[0].end <= left->pos[1].end &&
					left->pos[0].end < right->pos[1].begin &&
					left->pos[0].end <= right->pos[0].begin &&
					left->pos[1].end <= right->pos[1].begin) {
					candidate[0].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
				}
				else if (par->detect[5]) {
					// Inverted insertion from upstream
					if (right->pos[1].begin <= left->pos[1].end &&
						left->pos[0].end < right->pos[0].begin &&
						right->pos[0].begin < left->pos[1].end) {
						candidate[1].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
					}
					// Inverted insertion from downstream
					else if (right->pos[0].begin <= left->pos[0].end &&
						left->pos[1].end < right->pos[1].begin &&
						right->pos[0].begin < left->pos[1].end) {
						candidate[2].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
					}
				}
			}
			++right;
		}
		++par->status.current_task[t];
		++left;
	}
}
void stk::getInvCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv;
	srange outrange, inrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		outrange.begin = $_.first->pos[0].end + 1;
		outrange.end = $_.second->pos[1].begin - 1;
		inrange.begin = $_.second->pos[0].begin;
		inrange.end = $_.first->pos[1].end;
		if (inrange.length(true) + 1 < outrange.length(true)) del = true;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(outrange.begin - par->varp.svp.break_site_len, outrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(outrange.end + 1, outrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth		
		if (hasCtrl) {
			bg = cna->count($_.first->pos[0].idx, outrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;

		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, outrange);
		cp2 = cna->copy($_.first->pos[0].idx, inrange);
		if (del) {
			cp1 = (cp1 * outrange.length(true) - cp2 * inrange.length(true)) / (outrange.length(true) - inrange.length(true));
			cnv = par->varp.cnvp.evaluate(cp1);
			if ((cnv & 0xFF) != DELETION) continue;
		}
		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = INVERSION | (del ? DELETION : 0);
		var->genotype = (del ? (cnv >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, outrange.begin + 1, outrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[0].idx, inrange.begin + 1, inrange.end + 1, true);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, outrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[0].idx, inrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[0].idx, inrange, true) : 1.0f;
		//
		primary->add(var);
	}
}
void stk::getInvInsCpy1(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.first->pos[0].end + 1 < $_.second->pos[0].begin) {
			delrange.begin = $_.first->pos[0].end + 1;
			delrange.end = $_.second->pos[0].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.first->pos[0].end;
			delrange.end = $_.second->pos[0].begin;
		}
		//
		insrange.begin = $_.second->pos[1].begin;
		insrange.end = $_.first->pos[1].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.first->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[0].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;
		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[0].idx, insrange.begin + 1, insrange.end + 1, true);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[0].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[0].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}
void stk::getInvInsCpy2(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.first->pos[1].end + 1 < $_.second->pos[1].begin) {
			delrange.begin = $_.first->pos[1].end + 1;
			delrange.end = $_.second->pos[1].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.first->pos[1].end;
			delrange.end = $_.second->pos[0].begin;
		}
		//
		insrange.begin = $_.second->pos[0].begin;
		insrange.end = $_.first->pos[0].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.first->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[0].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;
		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[0].idx, insrange.begin + 1, insrange.end + 1, true);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[0].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[0].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}

// Simple translocation
//   lp0 --->|===> lp1
//   rp0 ===>:---> rp1
//
// Insertion from different chrom. with higher index
//    lp0 --->| :---> rp1
//       lp1 |===> ... ===>: rp0
//
// Insertion from different chrom. with lower index
//    rp1 :---> ... --->| lp0
//        rp0 ===>: |===> lp1
void stk::selectTrsCandidate(int t, int i1, int i2, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants1, Array<SVar>* variants2, stk::Param* par) {
	if (variants1->empty()) return;
	//
	auto left = variants1->begin(), rbeg = variants2->begin(), right = variants2->begin(), lend = variants1->end(), rend = variants2->end();
	//
	while (left < lend && left->pos[1].idx != i2) { ++left; }
	if (left == lend) return;
	--lend;
	while (left < lend && lend->pos[1].idx != i2) { --lend; }
	++lend;
	//
	while (rbeg < rend && rbeg->pos[1].idx != i1) { ++rbeg; }
	--rend;
	while (rbeg < rend && rend->pos[1].idx != i1) { --rend; }
	++rend;
	//
	if (rbeg == rend) {
		par->status.current_task[t] = lend - left;
		return;
	}
	//
	while (left < lend) {
		right = rbeg;
		while (right < rend) {
			if (left->pos[0].idx != right->pos[1].idx || left->pos[1].idx != right->pos[0].idx) break;
			if (check2(*left, *right, &par->varp)) {
				// Simple translocation
				if (par->detect[6] &&
					left->pos[0].end < right->pos[1].begin && right->pos[0].begin < left->pos[1].begin) {
					candidate[0].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
				}
				else if (par->detect[7]) {
					//  Insertion from different chrom. with higher index
					if (left->pos[0].end < right->pos[1].begin &&
						left->pos[1].begin <= right->pos[0].begin && left->pos[1].end <= right->pos[0].end) {
						candidate[1].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
					}
					// Insertion from different chrom.  with lower index
					else if (right->pos[0].end < left->pos[1].begin &&
						right->pos[1].begin <= left->pos[0].begin && right->pos[1].end <= left->pos[0].end) {
						candidate[2].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
					}
				}
			}
			++right;
		}
		++par->status.current_task[t];
		++left;
	}
}
void stk::getTrsCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg1, bg2, freq;
	sushort cnv1, cnv2;
	srange lowrange, highrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del[2] = { false, false };
		if ($_.first->pos[0].end + 1 < $_.second->pos[1].begin) {
			lowrange.begin = $_.first->pos[0].end + 1;
			lowrange.end = $_.second->pos[1].begin - 1;
			del[0] = true;
		}
		else {
			lowrange.begin = $_.first->pos[0].end;
			lowrange.end = $_.second->pos[1].begin;
		}
		if ($_.second->pos[0].end + 1 < $_.first->pos[1].begin) {
			highrange.begin = $_.second->pos[0].end + 1;
			highrange.end = $_.first->pos[1].begin - 1;
			del[1] = true;
		}
		else {
			highrange.begin = $_.second->pos[0].end;
			highrange.end = $_.first->pos[1].begin;
		}
		// Freq. check
		freq = (2.f * (float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(lowrange.begin - par->varp.svp.break_site_len, lowrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(lowrange.end + 1, lowrange.end + par->varp.svp.break_site_len)) +
			cna->count($_.second->pos[0].idx, srange(highrange.begin - par->varp.svp.break_site_len, highrange.begin - 1)) +
				cna->count($_.second->pos[0].idx, srange(highrange.end + 1, highrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth		
		if (hasCtrl) {
			bg1 = cna->count($_.first->pos[0].idx, lowrange, true);
			if (bg1 < par->varp.cnvp.min_bg) continue;
			bg2 = cna->count($_.second->pos[0].idx, highrange, true);
			if (bg2 < par->varp.cnvp.min_bg) continue;
		}
		else {
			bg1 = 1.f; bg2 = 1.f;
		}
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, lowrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.second->pos[0].idx, highrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del[0] && (cnv1 & 0xFF) != DELETION) continue;
		if (del[1] && (cnv2 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = TRANSLOCATION | ((del[0] || del[1]) ? DELETION : 0);
		var->genotype = ((del[0] || del[1]) ? 
			((del[0] ? (cnv1 >> 8) : 0) | (del[1] ? (cnv1 >> 8) : 0)) : 
			(par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, lowrange.begin + 1, lowrange.end + 1);
		var->pos[1] = sbpos($_.second->pos[0].idx, highrange.begin + 1, highrange.end + 1);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, lowrange, false);
		var->depth[0][1] = bg1;
		var->depth[1][0] = cna->count($_.second->pos[0].idx, highrange, false);
		var->depth[1][1] = bg2;
		//
		primary->add(var);
	}
}
void stk::getTrsInsCpy1(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.first->pos[0].end + 1 < $_.second->pos[1].begin) {
			delrange.begin = $_.first->pos[0].end + 1;
			delrange.end = $_.second->pos[1].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.first->pos[0].end;
			delrange.end = $_.second->pos[1].begin;
		}
		//
		insrange.begin = $_.first->pos[1].begin;
		insrange.end = $_.second->pos[0].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.first->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.second->pos[0].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = TRANSLOCATION | INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.second->pos[0].idx, insrange.begin + 1, insrange.end + 1);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.second->pos[0].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.second->pos[0].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}
void stk::getTrsInsCpy2(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.second->pos[0].end + 1 < $_.first->pos[1].begin) {
			delrange.begin = $_.second->pos[0].end + 1;
			delrange.end = $_.first->pos[1].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.second->pos[0].end;
			delrange.end = $_.first->pos[1].begin;
		}
		//
		insrange.begin = $_.second->pos[1].begin;
		insrange.end = $_.first->pos[0].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.second->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.second->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.second->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.second->pos[0].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[0].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second);
		var->flag = SR_VARIANT;
		var->type = TRANSLOCATION | INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.second->pos[0].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[0].idx, insrange.begin + 1, insrange.end + 1);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.second->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[0].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[0].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}

// Simple inverted translocation
//  rp0 <---|===> rp1 <complement> rp1' <===|---> rp0' 
//  lp0  --->|<=== lp1
//  rp1' <===:---> rp0'
//  
// Inverted insertion from different chrom. with higher index
//  rp0 <---|===> rp1 <complement> rp1' <===|---> rp0'   
//    lp0 --->| :---> rp0'
//          lp1 |<=== ... <===: rp1'
//
// Inverted insertion from different chrom.  with lower index
//  lp0  --->|<=== lp1 <complement> lp1' ===>|<--- lp0'  
//    lp0' |<--- ... <---: rp0
//        lp1' ===>| :===> rp1
void stk::selectTrInvCandidate(int t, int i1, int i2, Array<Pair<SVar*, SVar*>>* candidate, Array<SVar>* variants, stk::Param* par) {
	if (variants->empty()) return;
	variants->sort(sbio::sutil::svsorter(variants->at(0).type & 0xFF));
	auto left = variants->begin(), rbeg = variants->begin(), right = variants->begin(), lend = variants->end(), rend = variants->end();	
	//
	while (left < lend && left->pos[1].idx != i2) { ++left; }
	if (left == lend) return;
	--rend;
	while (left < rend && rend->pos[1].idx != i2) { --rend; }
	++rend;
	rbeg = left;
	while (rbeg < rend && !rbeg->pos[0].dir) { ++rbeg; }
	lend = rbeg;

	if (left->pos[0].dir || lend == rend) {
		par->status.current_task[t] = rend - left;
		return;
	}
	else par->status.current_task[t] = rend - right;

	while (left < lend) {
		right = rbeg;
		while (right < rend) {
			if (left->pos[0].idx != right->pos[0].idx || left->pos[1].idx != right->pos[1].idx) break;
			if (check2(*left, *right, &par->varp)) {
				// Simple inverted translocation
				if (par->detect[8] &&
					left->pos[0].end < right->pos[0].begin && left->pos[1].end < right->pos[1].begin) {
					candidate[0].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
				}
				//  Inverted insertion from different chrom. with higher index
				else if (par->detect[9]) {
					if (left->pos[0].end < right->pos[0].begin &&
						right->pos[1].begin <= left->pos[1].begin && right->pos[1].end <= left->pos[1].end) {
						candidate[1].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
					}
					// Inverted insertion from different chrom.  with lower index
					else if (left->pos[1].end < right->pos[1].begin &&
						right->pos[0].begin <= left->pos[0].begin && right->pos[0].end <= left->pos[0].end) {
						candidate[2].add(Pair<SVar*, SVar*>(left.ptr(), right.ptr()));
					}
				}
			}
			++right;
		}
		++par->status.current_task[t];
		++left;
	}
}
void stk::getTrsInvCpy(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg1, bg2, freq;
	sushort cnv1, cnv2;
	srange lowrange, highrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del[2] = { false, false };
		if ($_.first->pos[0].end + 1 < $_.second->pos[0].begin) {
			lowrange.begin = $_.first->pos[0].end + 1;
			lowrange.end = $_.second->pos[0].begin - 1;
			del[0] = true;
		}
		else {
			lowrange.begin = $_.first->pos[0].end;
			lowrange.end = $_.second->pos[0].begin;
		}
		if ($_.first->pos[1].end + 1 < $_.second->pos[1].begin) {
			highrange.begin = $_.first->pos[1].end + 1;
			highrange.end = $_.second->pos[1].begin - 1;
			del[1] = true;
		}
		else {
			highrange.begin = $_.first->pos[1].end;
			highrange.end = $_.second->pos[1].begin;
		}
		// Freq. check
		freq = (2.f * (float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(lowrange.begin - par->varp.svp.break_site_len, lowrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(lowrange.end + 1, lowrange.end + par->varp.svp.break_site_len)) +
				cna->count($_.first->pos[1].idx, srange(highrange.begin - par->varp.svp.break_site_len, highrange.begin - 1)) +
				cna->count($_.first->pos[1].idx, srange(highrange.end + 1, highrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg1 = cna->count($_.first->pos[0].idx, lowrange, true);
			if (bg1 < par->varp.cnvp.min_bg) continue;
			bg2 = cna->count($_.first->pos[1].idx, highrange, true);
			if (bg2 < par->varp.cnvp.min_bg) continue;
		}
		else { bg1 = 1.f; bg2 = 1.f; }
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, lowrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[1].idx, highrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del[0] && (cnv1 & 0xFF) != DELETION) continue;
		if (del[1] && (cnv2 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second, 1);
		var->flag = SR_VARIANT;
		var->type = INVERSION | TRANSLOCATION | ((del[0] || del[1]) ? DELETION : 0);
		var->genotype = ((del[0] || del[1]) ?
			((del[0] ? (cnv1 >> 8) : 0) | (del[1] ? (cnv1 >> 8) : 0)) :
			(par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, lowrange.begin + 1, lowrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[1].idx, highrange.begin + 1, highrange.end + 1, true);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, lowrange, false);
		var->depth[0][1] = bg1;
		var->depth[1][0] = cna->count($_.first->pos[1].idx, highrange, false);
		var->depth[1][1] = bg2;
		//
		primary->add(var);
	}
}
void stk::getTrsInvInsCpy1(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.first->pos[0].end + 1 < $_.second->pos[0].begin) {
			delrange.begin = $_.first->pos[0].end + 1;
			delrange.end = $_.second->pos[0].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.first->pos[0].end;
			delrange.end = $_.second->pos[0].begin;
		}
		//
		insrange.begin = $_.second->pos[1].begin;
		insrange.end = $_.first->pos[1].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[0].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.first->pos[0].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.first->pos[0].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.first->pos[0].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[1].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second, 1);
		var->flag = SR_VARIANT;
		var->type = INVERSION | TRANSLOCATION | INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[0].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[1].idx, insrange.begin + 1, insrange.end + 1, true);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[0].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[1].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[1].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}
void stk::getTrsInvInsCpy2(Array<SPointer<sbio::Variant>>* primary, Array<Pair<SVar*, SVar*>>* candidates, CNAnalysis* cna, stk::Param* par) {
	if (candidates->empty()) return;
	float cp1, cp2, bg, freq;
	sushort cnv1, cnv2;
	srange delrange, insrange;
	bool hasCtrl = par->control.isLoaded();
	sfor(*candidates) {
		//
		bool del = false;
		if ($_.first->pos[1].end + 1 < $_.second->pos[1].begin) {
			delrange.begin = $_.first->pos[1].end + 1;
			delrange.end = $_.second->pos[1].begin - 1;
			del = true;
		}
		else {
			delrange.begin = $_.first->pos[1].end;
			delrange.end = $_.second->pos[1].begin;
		}
		//
		insrange.begin = $_.second->pos[0].begin;
		insrange.end = $_.first->pos[0].end;
		// Freq. check
		freq = ((float)($_.first->total() + $_.second->total())) /
			(cna->count($_.first->pos[1].idx, srange(delrange.begin - par->varp.svp.break_site_len, delrange.begin - 1)) +
				cna->count($_.first->pos[1].idx, srange(delrange.end + 1, delrange.end + par->varp.svp.break_site_len)));
		if (freq < par->varp.svp.min_freq) continue;
		// Background depth
		if (hasCtrl) {
			bg = cna->count($_.first->pos[1].idx, delrange, true);
			if (bg < par->varp.cnvp.min_bg) continue;
		}
		else bg = 1.f;
		// Copy value check
		cp1 = cna->copy($_.first->pos[1].idx, delrange);
		cnv1 = par->varp.cnvp.evaluate(cp1);
		cp2 = cna->copy($_.first->pos[0].idx, insrange);
		cnv2 = par->varp.cnvp.evaluate(cp2);
		if (del && (cnv1 & 0xFF) != DELETION) continue;

		// Add to primary variant list
		Variant* var = new Variant(*$_.first, *$_.second, 2);
		var->flag = SR_VARIANT;
		var->type = INVERSION | TRANSLOCATION | INSERTION | (del ? DELETION : 0) |
			(cnv2 & DUPLICATION ? DUPLICATION : 0) | (cnv2 & MULTIPLICATION ? MULTIPLICATION : 0);
		var->genotype = (del ? (cnv1 >> 8) : (par->varp.svp.homo_freq <= freq ? HOMO_VAR : HETERO_VAR));
		var->pos[0] = sbpos($_.first->pos[1].idx, delrange.begin + 1, delrange.end + 1);
		var->pos[1] = sbpos($_.first->pos[0].idx, insrange.begin + 1, insrange.end + 1, true);
		var->frequency = freq;
		var->copy[0] = cp1;
		var->copy[1] = cp2;
		var->depth[0][0] = cna->count($_.first->pos[1].idx, delrange, false);
		var->depth[0][1] = bg;
		var->depth[1][0] = cna->count($_.first->pos[0].idx, insrange, false);
		var->depth[1][1] = hasCtrl ? cna->count($_.first->pos[0].idx, insrange, true) : 1.f;
		//
		primary->add(var);
	}
}