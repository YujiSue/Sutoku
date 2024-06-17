#include "smath/stat.h"
#include "param.h"
stk::Status::Status() : state(0), total_task(0) {
/*
#ifdef SUB_PROCESS
	smem = SharedMemory(1024, "stkm.dat");
#endif
*/
}
stk::Status::~Status() {}
void stk::Status::reset() {
	message.clear();
	total_task = 0;
	current_task.clear();
}
void stk::Status::setState(const int st) {
	SAutoLock al(_lock);
	state = st;
}
void stk::Status::setTask(const size_t total, const size_t count) {
	total_task = total;
	current_task.resize(count);
	current_task.reset(0);
}
void stk::Status::setMessage(const char* s) {
	SAutoLock al(_lock);
	message = s;
}
double stk::Status::progress() {
	if (total_task) return (double)sstat::sum(current_task) / total_task;
	else return 0.0;
}
/*
void Status::share() {
#ifdef SUB_PROCESS
	size_t off = 0;
	smem.set(&state, sizeof(sint), off); off += sizeof(sint);
	double p = progress();
	smem.set(&p, sizeof(double), off); off += sizeof(double);
	smem.set(msg.cstr(), msg.length(), off);
#endif
}
*/

void stk::showProgress(stk::Status* status) {
	int prog;
	while (status->state <= stk::RUNNING) {
		prog = (int)(status->progress() * 100.0);
		SWrite(DEL * 4, sstr::lfill(String(prog), ' ', 3), "%");
		slib::sutil::sleep(1000);
	}
}
void stk::closeThread(std::thread* thread, stk::Status* status) {
	if (status->state == stk::RUNNING) status->setState(stk::FINISHED);
	thread->join();
	if (status->state != stk::ERRORED) {
		SPrint(DEL * 4, "100% ... Completed.");
	}
}
void stk::interruptProc(stk::Status* status, SLogger* logger, Exception* ex) {
	if (status->state == stk::RUNNING) {
		status->setState(stk::ERRORED);
		SPrint(" ... Interrupted.");
	}
	logger->log(*ex);
}
bool stk::headClip(CigarArray* cigars, size_t sz) {
	if (cigars->empty()) return false;
	auto& cig = cigars->at(0);
	return ((cig.option == scigar::HCLIP || cig.option == scigar::SCLIP) && sz <= cig.length);
}
bool stk::tailClip(CigarArray* cigars, size_t sz) {
	if (cigars->empty()) return false;
	auto& cig = cigars->at(-1);
	return ((cig.option == scigar::HCLIP || cig.option == scigar::SCLIP) && sz <= cig.length);
}
double stk::pairQual(sbam::ReadInfo* read, sbam::ReadInfo* pair, stk::Param* par) {
	auto q1 = sbio::sutil::rawVal(read->mapq, par->min_err_prob), q2 = sbio::sutil::rawVal(pair->mapq, par->min_err_prob);
	auto pq = 1.0 - ((1.0 - q1) * (1.0 - q2));
	return pq < par->min_err_prob ? par->min_err_prob : pq;
}
// 
void stk::countLowCover(int r, NGSData* data, suinteger* count) {
	auto dp = data->depth[r].data();
	auto size = data->depth[r].size();
	auto bin = data->summary.bin;
	sforin(i, 0, size - 1) {
		if ((*dp) < 1.f) (*count) += (int)((1.f - (*dp)) * bin);
		++dp;
	}
	auto rest = data->summary.reflen[r] - (int)((size - 1) * bin);
	if ((*dp) < 1.f) (*count) += (int)((1.f - (*dp)) * rest);
}
// 
void stk::checkVariants(Array<SVar>* variants, stk::Param* par) {
	if (variants->empty()) return;
	auto sz = variants->size();
	auto sorter = slib::sbio::sutil::svsorter(variants->at(0).type & 0xFF);
	variants->sort(sorter);
	sfor(*variants) {
		if ($_.type == 0) continue;
		auto current = $, compare = $ + 1;
		while (compare < variants->end()) {
			if (!current->comparable(*compare) || current->lt(*compare, &par->varp.svp)) break;
			if (current->equal(*compare, &par->reference, &par->varp.svp)) {
				current->merge(*compare, &par->varp.svp); --sz;
				if (current->type == 0) current = compare;
			}
			++compare;
		}
	}
	variants->sort(sorter);
	variants->resize(sz);
}
inline void addDepth(int r, slib::sbio::NGSData* data1, slib::sbio::NGSData* data2) {
	auto sz = data1->depth[r].size();
	auto dp1 = data1->depth[r].data(), dp2 = data2->depth[r].data();
	sforin(i, 0, sz) { (*dp1) += (*dp2); ++dp1; ++dp2; }
}
void stk::integrate(NGSData* data1, NGSData* data2, stk::Param* par) {
	data1->summary.avelen *= sstat::sum(data1->summary.count);
	data1->summary.avelen += data2->summary.avelen * sstat::sum(data2->summary.count);
	data1->summary.total += data2->summary.total;
	auto refnum = data1->summary.refnum;
	svecu uncov(refnum, 0);
	sforin(r, 0, refnum) {
		data1->summary.bases[r] *= data1->summary.count[r];
		data1->summary.bases[r] += data2->summary.bases[r] * data2->summary.count[r];
		data1->summary.count[r] += data2->summary.count[r];
		data1->summary.bases[r] /= data1->summary.count[r];
		//
		addDepth(r, data1, data2);
		// 
		stk::countLowCover(r, data1, &uncov[r]);
		//
		sforin(v, 0, SV_TYPE_COUNT) {
			data1->variants[SV_TYPE_COUNT * r + v].append(data2->variants[SV_TYPE_COUNT * r + v]);
			par->threads.addTask(stk::checkVariants, &data1->variants[SV_TYPE_COUNT * r + v], par);
		}
	}
	size_t tlen = par->reference.total(), len = 0;
	if (par->target.size()) {
		sforeach(target, par->target) len += target.length(true);
	}
	else len = tlen;
	data1->summary.avedp = data1->summary.avelen / tlen;
	data1->summary.avelen /= sstat::sum(data1->summary.count);
	par->threads.complete();
	data1->summary.cover = 1.0 - ((double)sstat::sum(uncov) / tlen);
}

void stk::subtract(Array<SVar>* variants1, Array<SVar>* variants2, stk::Param* par) {
	if (variants1->empty() || variants2->empty()) return;
	size_t sz = variants1->size();
	auto sorter = slib::sbio::sutil::svsorter(variants1->at(0).type & 0xFF);
	variants1->sort(sorter);
	variants2->sort(sorter);
	//
	sfor(*variants1) {
		auto compare = variants2->begin();
		if ($_.pos[1].idx != compare->pos[1].idx ||
			$_.pos[0].dir != compare->pos[0].dir) {
			while (compare < variants2->end() && !$_.comparable(*compare)) { ++compare; }
		}
		while (compare < variants2->end()) {
			if ($_.lt(*compare, &par->varp.svp)) break;
			else if($_.comparable(*compare) && 
				$_.equal(*compare, &par->reference, &par->varp.svp)) { 
				$_.type = 0; --sz; break;
			}
			++compare;
		}
	}
	variants1->sort(sorter);
	variants1->resize(sz);
}