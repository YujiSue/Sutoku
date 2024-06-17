#include "analyzer.h"
using namespace stk;

stk::VarSearch::VarSearch(Analyzer* an) {
	par = &an->par;
	status = &par->status;
	logger = &par->logger;
	variants.setReference(&an->par.reference);
	threads = &par->threads;
	filter = VarFilter(&par->reference, &par->annotdb, &par->varp, &par->target);
}
stk::VarSearch::~VarSearch() {}
/**
* Get primary SV candidates
*/
void stk::VarSearch::selectCandidate(NGSData* data) {
	auto num = data->summary.refnum;
	auto cand1 = candidates1.begin();
	auto cand2 = candidates2.begin();
	auto task = 0;
	sforin(r, 0, num) {
		threads->addTask(selectDelCandidate, task, cand1.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DEL], par);
		++cand1; ++task;
		threads->addTask(selectDupCandidate, task, cand1.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DUP], par);
		++cand1; ++task;
		threads->addTask(selectInsCandidate, task, cand1.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::INS], par);
		++cand1; ++task;
		threads->addTask(selectComplexCandidate, task, cand2.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DEL], &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DUP], par);
		++cand2; ++task;
		threads->addTask(selectComplexCandidate, task, cand2.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DUP], &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DEL], par);
		++cand2; ++task;
		threads->addTask(selectInvCandidate, task, cand2.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::INV], par);
		cand2 += 3; ++task;
		//
		sforin(s, r + 1, num) {
			threads->addTask(selectTrsCandidate, task, r, s, cand2.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::TRS], &data->variants[s * SV_TYPE_COUNT + (int)SVREAD_TYPE::TRS], par);
			cand2 += 3; ++task;
			//
			threads->addTask(selectTrInvCandidate, task, r, s, cand2.ptr(), &data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::TRS], par);
			cand2 += 3;
		}
	}
	threads->complete();
}
/**
* Filter by copy number
*/
void stk::VarSearch::copyCheck(NGSData* data) {
	primary.resize(candidates1.size() + candidates2.size());
	auto prim = primary.begin();
	auto cand1 = candidates1.begin();
	auto cand2 = candidates2.begin();
	sforin(r, 0, data->summary.refnum) {
		//threads->addTask(getDelCpy, prim.ptr(), cand1.ptr(), &cna, par);
		getDelCpy(prim.ptr(), cand1.ptr(), &cna, par);
		++cand1; ++prim;
		//threads->addTask(getDupCpy, prim.ptr(), cand1.ptr(), &cna, par);
		getDupCpy(prim.ptr(), cand1.ptr(), &cna, par);
		++cand1; ++prim;
		//threads->addTask(getInsCpy, prim.ptr(), cand1.ptr(), &cna, par);
		getInsCpy(prim.ptr(), cand1.ptr(), &cna, par);
		++cand1; ++prim;
		//threads->addTask(getComplexCpy, prim.ptr(), cand2.ptr(), &cna, par);
		getComplexCpy(prim.ptr(), cand2.ptr(), &cna, par);
		++cand2; ++prim;
		//threads->addTask(getComplexCpy, prim.ptr(), cand2.ptr(), &cna, par);
		getComplexCpy(prim.ptr(), cand2.ptr(), &cna, par);
		++cand2; ++prim;
		//threads->addTask(getInvCpy, prim.ptr(), cand2.ptr(), &cna, par);
		getInvCpy(prim.ptr(), cand2.ptr(), &cna, par);
		++cand2; ++prim;
		//threads->addTask(getInvInsCpy1, prim.ptr(), cand2.ptr(), &cna, par);
		getInvInsCpy1(prim.ptr(), cand2.ptr(), &cna, par);
		++cand2; ++prim;
		//threads->addTask(getInvInsCpy2, prim.ptr(), cand2.ptr(), &cna, par);
		getInvInsCpy2(prim.ptr(), cand2.ptr(), &cna, par);
		++cand2; ++prim;
		sforin(s, r + 1, data->summary.refnum) {
			//threads->addTask(getTrsCpy, prim.ptr(), cand2.ptr(), &cna, par);
			getTrsCpy(prim.ptr(), cand2.ptr(), &cna, par);
			++cand2; ++prim;
			//threads->addTask(getTrsInsCpy1, prim.ptr(), cand2.ptr(), &cna, par);
			getTrsInsCpy1(prim.ptr(), cand2.ptr(), &cna, par);
			++cand2; ++prim;
			//threads->addTask(getTrsInsCpy2, prim.ptr(), cand2.ptr(), &cna, par);
			getTrsInsCpy2(prim.ptr(), cand2.ptr(), &cna, par);
			++cand2; ++prim;
			//threads->addTask(getTrsInvCpy, prim.ptr(), cand2.ptr(), &cna, par);
			getTrsInvCpy(prim.ptr(), cand2.ptr(), &cna, par);
			++cand2; ++prim;
			//threads->addTask(getTrsInvInsCpy1, prim.ptr(), cand2.ptr(), &cna, par);
			getTrsInvInsCpy1(prim.ptr(), cand2.ptr(), &cna, par);
			++cand2; ++prim;
			//threads->addTask(getTrsInvInsCpy2, prim.ptr(), cand2.ptr(), &cna, par);
			getTrsInvInsCpy2(prim.ptr(), cand2.ptr(), &cna, par);
			++cand2; ++prim;
		}
	}
	//threads->complete();
}
/*
* Filter conflict variant
*/
inline bool synVar(Variant *v1, Variant* v2, slib::sbio::SeqList *ref, stk::Param *par) {
	String s1, s2;
	if (v1->type == DELETION) {
		s1 = ref->at(v1->pos[0].idx).raw(v1->pos[0].begin - par->varp.svp.break_site_len, par->varp.svp.break_site_len) + v1->alt;
		s2 = ref->at(v2->pos[0].idx).raw(v2->pos[0].begin - par->varp.svp.break_site_len, par->varp.svp.break_site_len) + v2->alt;
		if (s1.size() != par->varp.svp.break_site_len) s1.clip(s1.size() - par->varp.svp.break_site_len);
		if (s2.size() != par->varp.svp.break_site_len) s2.clip(s2.size() - par->varp.svp.break_site_len);
		if (par->varp.svp.max_dist < slib::smath::levenshtein(&s1[0], s1.size(), &s2[0], s2.size())) return false;
		s1 = v1->alt + ref->at(v1->pos[0].idx).raw(v1->pos[0].end - 1, par->varp.svp.break_site_len);
		s2 = v2->alt + ref->at(v2->pos[0].idx).raw(v2->pos[0].end - 1, par->varp.svp.break_site_len);
		if (s1.size() != par->varp.svp.break_site_len) s1.resize(par->varp.svp.break_site_len);
		if (s2.size() != par->varp.svp.break_site_len) s2.resize(par->varp.svp.break_site_len);
		if (par->varp.svp.max_dist < slib::smath::levenshtein(s1.cstr(), s1.size(), s2.cstr(), s2.size())) return false;
		return true;
	}
	else if (v1->type == DUPLICATION) {
		s1 = ref->at(v1->pos[0].idx).raw(v1->pos[0].end - par->varp.svp.break_site_len, par->varp.svp.break_site_len) + v1->alt;
		s2 = ref->at(v2->pos[0].idx).raw(v2->pos[0].end - par->varp.svp.break_site_len, par->varp.svp.break_site_len) + v2->alt;
		if (s1.size() != par->varp.svp.break_site_len) s1.clip(s1.size() - par->varp.svp.break_site_len);
		if (s2.size() != par->varp.svp.break_site_len) s2.clip(s2.size() - par->varp.svp.break_site_len);
		if (par->varp.svp.max_dist < slib::smath::levenshtein(s1.cstr(), s1.size(), s2.cstr(), s2.size())) return false;
		s1 = v1->alt + ref->at(v1->pos[0].idx).raw(v1->pos[0].begin - 1, par->varp.svp.break_site_len);
		s2 = v2->alt + ref->at(v2->pos[0].idx).raw(v2->pos[0].begin - 1, par->varp.svp.break_site_len);
		if (s1.size() != par->varp.svp.break_site_len) s1.resize(par->varp.svp.break_site_len);
		if (s2.size() != par->varp.svp.break_site_len) s2.resize(par->varp.svp.break_site_len);
		if (par->varp.svp.max_dist < slib::smath::levenshtein(s1.cstr(), s1.size(), s2.cstr(), s2.size())) return false;
		return true;
	}
	return false;
}
inline bool transHetero(Variant* v1, Variant* v2, stk::Param* par) {
	if (v1->type == DELETION || v1->type == DUPLICATION) {
		if (v1->pos[0].include(v2->pos[0]) &&
			(v1->genotype & HETERO_VAR)&&
			(v2->genotype & HOMO_VAR)) return true;
		else if (v2->pos[0].include(v1->pos[0]) &&
			(v2->genotype & HETERO_VAR) &&
			(v1->genotype & HOMO_VAR)) return true;
		else if ((v1->genotype & HETERO_VAR) &&
			(v2->genotype & HETERO_VAR)) return true;
	}
	return false;
}
void stk::VarSearch::conflictCheck() {
	// For germline variant
	Array<FreePointer<Variant>> tmp;
	sfor(variants) {
		if (($_->flag & NOT_USE_FLAG) || ($_->flag & UNAVAILABLE_FLAG)) continue;
		Variant* var = $_;
		tmp.clear();
		sforin(nxt, $ + 1, variants.end()) {
			if (((*nxt)->flag & NOT_USE_FLAG) || ((*nxt)->flag & UNAVAILABLE_FLAG)) continue;
			if (var->pos[0].idx != (*nxt)->pos[0].idx || var->pos[0].end < (*nxt)->pos[0].begin) break;
			if (var->type == (*nxt)->type && var->pos[0].overlap((*nxt)->pos[0])) {
				if (synVar(var, *nxt, &par->reference, par)) {
					if (var->qual < (*nxt)->qual) {
						var->flag |= NOT_USE_FLAG; break;
					}
					else (*nxt)->flag |= NOT_USE_FLAG;
				}
				else tmp.add(*nxt);
			}
		}
		// 
		if (1 < tmp.size()) {
			auto maxv = var;
			sforeach(candidate, tmp) {
				if (maxv->qual < candidate->qual) {
					maxv->flag |= NOT_USE_FLAG; maxv = candidate;
				}
			}
		}
		else if (tmp.size() == 1) {
			if (transHetero(var, tmp[0], par)) {
				var->genotype = TRANS_HETERO_VAR;
				tmp[0]->genotype = TRANS_HETERO_VAR;
			}
			else {
				if (var->qual < tmp[0]->qual) var->flag |= NOT_USE_FLAG;
				else tmp[0]->flag |= NOT_USE_FLAG;
			}
		}
	}
}
/*
* Gene and reported variant annotation
*/
void stk::VarSearch::annotate() {
	sforeach(var, variants) {
		if ((var->flag & NOT_USE_FLAG) || (var->flag & UNAVAILABLE_FLAG)) continue;
		par->annotdb.annotate(*var, par->reference, par->varp);
	}
}

inline void resizeCandidates(NGSData* data, Array<Array<SVar*>>* candidates1, Array<Array<Pair<SVar*, SVar*>>>* candidates2, stk::Param *par) {
	int n = (int)par->reference.size(),
		n2 = smath::combination(n, 2);
 	candidates1->resize(n * 3);
	candidates2->resize(n * 5 + n2 * 6);
	size_t total = 0;
	sforin(r, 0, n) {
		total += data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DEL].size() * 2;
		total += data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::DUP].size() * 2;
		total += data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::INS].size();
		total += data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::INV].size();
		total += data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::TRS].size();
		total += data->variants[r * SV_TYPE_COUNT + (int)SVREAD_TYPE::TRINV].size();
	}
	par->status.setTask(total, n * 8 + n2 * 6);
}
inline void setVarID(slib::sbio::VarList& vlist) {
	auto idx = 1;
	sfor(vlist) {
		if ($_->varid.empty()) $_->varid = "variant_" + S(idx);
		++idx;
	}
}
void stk::VarSearch::detect(NGSData* data) {
	try {
		// Status => Running
		status->setState(stk::RUNNING);
		//
		variants.attribute["detect-type"] = "SV";
		variants.attribute["_prog_"] = { "slib", "Sutoku" };
		// Log
		logger->log("Started to detect structure variants.");
		//
		//VarFilter filter(&par->reference, &par->annotdb, &par->varp, &par->target);
		bool hasCtrl = par->control.isLoaded();
		if (hasCtrl) {
			// Log
			logger->log("Subtraction of background variants.");
			// Subtraction
			sfor2(data->variants, par->control.variants) threads->addTask(stk::subtract, &$_1, &$_2, par);
			threads->complete();
			logger->log("Completed.");
		}
		// Log
		logger->log("Primary detection.");
		//
		cna.setParam(par);
		cna.setData(data, (hasCtrl ? &par->control : nullptr));
		cna.analyze();
		//
		resizeCandidates(data, &candidates1, &candidates2, par);
		// Display progress
		SWrite(SP * 4, "> Progress: ", SP * 4);
		std::thread thread1(stk::showProgress, status);
		//
		selectCandidate(data);
		//
		stk::closeThread(&thread1, status);
		// Log
		logger->log("Copy number check.");
		// Display progress
		SWrite(SP * 4, "> Progress: ", SP * 4);
		std::thread thread2(stk::showProgress, status);
		//
		copyCheck(data);
		//
		stk::closeThread(&thread2, status);
		// Log
		logger->log("Integration.");
		// 
		size_t count = 0;
		sfor(primary) count += $_.size();
		variants.reserve(count + 1);
		sfor(primary) variants.append($_);
		//
		conflictCheck();
		//
		if (par->annotation) {
			// Log
			logger->log("Annotation.");
			annotate();
			logger->log("Completed.");
		}
		//
		filter.filter(variants);
		//if (par->filters)
		// 
		variants.tidyUp();
		setVarID(variants);
	}
	catch (Exception ex) {
		logger->log(ex);
		status->setState(stk::ERRORED);
	}
}

template<class Content>
inline void clearCandidates(Array<Content>* array) { array->clear(); }
void stk::VarSearch::reset() {
	if (threads->isWorking()) threads->complete();
	cna.reset();
	sfor(candidates1) threads->addTask(clearCandidates<SVar*>, $.ptr());
	sfor(candidates2) threads->addTask(clearCandidates<Pair<SVar*, SVar*>>, $.ptr());
	sfor(primary) threads->addTask(clearCandidates<SPointer<Variant>>, $.ptr());
	threads->complete();
	variants.clearAll();
}
