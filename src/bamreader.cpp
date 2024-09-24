#include "analyzer.h"
#include "smath/stat.h"
using namespace slib::smath;
//
stk::ReadCounter::ReadCounter() : bin(1), depths(nullptr) {}
stk::ReadCounter::ReadCounter(NGSData* data, Param* par) : ReadCounter() { set(data, par); }
stk::ReadCounter::~ReadCounter() {}
void stk::ReadCounter::set(NGSData* data, Param* par) {
	bin = data->summary.bin; depths = &data->depth;
}
void stk::ReadCounter::count(sbpos* pos, CigarArray* cigars) {
	auto rpos = pos->begin;
	float* depth = depths->at(pos->idx).data(), * dp;
	int beg, end;
	if (cigars) {
		sfor(*cigars) {
			beg = rpos / bin, end = (rpos + $_.length - 1) / bin;
			dp = depth + beg;
			if ($_.option == scigar::MATCH || $_.option == scigar::PMATCH || $_.option == scigar::MMATCH) {
				if (beg == end) { *dp += (float)$_.length / bin; }
				else {
					*dp += (float)(bin - (rpos % bin)) / bin; ++beg; ++dp;
					while (beg < end) {
						*dp += 1.f; ++beg; ++dp;
					}
					*dp += (float)((rpos + $_.length) - (end * bin)) / bin;
				}
				rpos += $_.length;
			}
			else if ($_.option == scigar::DELETION || $_.option == scigar::SKIP) rpos += $_.length;
		}
	}
	else {
		beg = rpos / bin, end = pos->end / bin;
		dp = depth + beg;
		if (beg == end) { *dp += (float)pos->length(true) / bin; }
		else {
			*dp += (float)(bin - (rpos % bin)) / bin; ++beg; ++dp;
			while (beg < end) {
				*dp += 1.f; ++beg; ++dp;
			}
			*dp += (float)((rpos + pos->length(true)) - (end * bin)) / bin;
		}
	}
}
void stk::ReadCounter::count(sbam::ReadInfo* read) { count(&read->ref, &read->cigars); }

/**
* BamReader implementation
*/
stk::BamReader::BamReader(stk::Analyzer *an) {
	par = &an->par;
	status = &par->status;
	logger = &par->logger;
	threads = &par->threads;
}
stk::BamReader::~BamReader() {}
/**
* Get read info as SAM format (brief check). For details, use samtools and other softwares.
*/
void stk::BamReader::readinfo(BamFile* bam, IOStream& stream) {
	sbam::ReadInfo* read;
	// All
	if (par->target.empty()) {
		while (read = bam->next()) { 
			if (!par->clipped || (par->clipped && read->cigars.clipped(par->min_clip)))
				stream << read->toString() << NL; stream.flush();
		}
	}
	// Targeted
	else {
		sfori(par->target) {
			auto &trgt = par->target[i];
			auto bins = sbam::getBins(trgt);
			vchunks chunks;
			sforeach(bin, bins) {
				chunks.append(bam->index.chunks[i][bam->index.bin_map[i][bin]]);
			}
			chunks.unique();
			sforeach(chunk, chunks) {
				auto current = chunk.begin;
				bam->setVOffset(current);
				while (current < chunk.end && (read = bam->next())) { 
					if (par->target[i].overlap(read->ref)) {
						if (!par->clipped || (par->clipped && read->cigars.clipped(par->min_clip))) {
							stream << read->toString() << NL; stream.flush();
						}
					}
					current = bam->voffset();
				}
			}
		}
	}
	logger->log("Finished.");
}
/**
* Extend alignment before re-alignment.
*/
void _extendRead(sbam::ReadInfo* read, Sequence* ref, AlignExtend* extender) {
	if (read->cigars[0].option == scigar::SCLIP) {
		read->query.begin += read->cigars[0].length;
		read->cigars.remove(0, 1);
	}
	if (read->cigars[-1].option == scigar::SCLIP) {
		read->query.end -= read->cigars[-1].length;
		read->cigars.resize(read->cigars.size() - 1);
	}
	bool rev = read->ref.dir;
	read->ref.dir = false;
	extender->extend(ref, &read->seq, read);
	if (read->query.begin) read->cigars.add(Cigar(scigar::SCLIP, read->query.begin), DIRECTION::HEAD);
	if (read->query.end < read->seq.size() - 1) read->cigars.add(Cigar(scigar::SCLIP, (int)read->seq.size() - read->query.end - 1));
	read->ref.dir = rev;
	read->query = srange(0, (int)read->seq.size() - 1);
	extender->reset();
}
// Single read summary
void stk::BamReader::summary1(int refidx, BamFile* bam, NGSData* data, SRAnalysis* svd, vchunk range, bool async) {
	stk::ReadCounter counter(data, par);
	sbam::ReadInfo* read = nullptr;
	AlignExtend extender(&par->seqp);
	//
	auto current = range.begin;
	try {
		// Move to offset
		bam->setVOffset(current);
		// Loop
		while ((current = bam->voffset()) < range.end &&
			(read = bam->next()) &&
			status->state == stk::RUNNING) {
			// Ignore secondary alignment
			if (read->flag & sbam::SECONDARY_ALIGN_READ || read->flag & sbam::SUPPLEMENTAL) continue;
			// Increment total read count
			if (async) {
				SAutoLock al(par->slock);
				++(data->summary.total);
			}
			else ++(data->summary.total);
			// Discard unmapped read
			if (read->flag & sbam::UNMAPPED_READ) continue;
			// PCR duplication check
			if (par->ignore_dp && (read->flag & sbam::PCR_DUPLICATE)) continue;
			// Update read count and total length
			++(data->summary.count[refidx]);
			data->summary.bases[refidx] += read->seq.size();
			// Extend
			_extendRead(read, &par->reference[read->ref.idx], &extender);
			// Check read with clipped region
			if (stk::headClip(&read->cigars, par->min_clip)) {
				//extender.extendHead...
				if (par->detect_sv) svd->addClipRead(read, DIRECTION::HEAD);
			}
			else if (stk::tailClip(&read->cigars, par->min_clip)) {
				//extender.extendTail...
				if (par->detect_sv) svd->addClipRead(read, DIRECTION::TAIL);
			}
			// Depth count
			counter.count(read);
			// Progress indicate
			status->current_task[async?refidx:0] = current.file_offset - range.begin.file_offset;
		}
		// Re-align leftovers
		if (par->detect_sv && status->state == stk::RUNNING) svd->realign();
	}
	catch (Exception ex) { 
		stk::interruptProc(status, logger, &ex); 
	}
}
inline void _runSummary1(stk::BamReader* br, int i, BamFile* bam, NGSData* summary, stk::SRAnalysis* svd, vchunk range) {
	br->summary1(i, bam, summary, svd, range, true);
}
// Correct sequenced site of paired reads
void _correctOverlap(sbam::ReadInfo* read, sbam::ReadInfo* pair) {
	if (read->ref.idx != pair->ref.idx) return;
	// Re-define 1st and 2nd reads 
	auto r1 = read->ref.dir ? pair : read;
	auto r2 = read->ref.dir ? read : pair;
	// Add empty cigar 
	if (r1->ref.end + 1 == r2->ref.begin) {
		if (r1->cigars[-1].option != scigar::SCLIP) r1->cigars.add(Cigar(scigar::SCLIP, 0));
		if (r2->cigars[0].option != scigar::SCLIP) r2->cigars.add(Cigar(scigar::SCLIP, 0), DIRECTION::HEAD);
	}
	if (r1->ref.end < r2->ref.begin || r2->ref.end < r1->ref.begin) return;
	// Use only one side read
	if (r2->ref.begin <= r1->ref.begin) { r1->cigars.clear(); return; }
	else if (r2->ref.end <= r1->ref.end) { r2->cigars.clear(); return; }
	// Check whether pairs has insertion 
	bool midclip = r1->cigars[-1].option == scigar::SCLIP && r2->cigars[0].option == scigar::SCLIP;
	// Trimming overlap
	sbio::sutil::trimOver(r1->ref, r1->query, r1->cigars, r2->ref, r2->query, r2->cigars, false);
	if (midclip && r1->cigars.size() && r1->query.end < r1->seq.size() - 1) {
		r1->cigars.add(Cigar(scigar::SCLIP, (int)r1->seq.size() - r1->query.end - 1));
	}
	else {
		// Add empty cigar to mark the read is 
		if (r1->cigars[-1].option != scigar::SCLIP) r1->cigars.add(Cigar(scigar::SCLIP, 0));
		if (r2->cigars[0].option != scigar::SCLIP) r2->cigars.add(Cigar(scigar::SCLIP, 0), DIRECTION::HEAD);
	}
}
// Paired end summary 
void stk::BamReader::summary2_(int refidx, BamFile* bam, NGSData* data, Array<Map<String, slib::Pointer<sbam::ReadInfo>>>* pairs, SRAnalysis* svd, vchunk range, bool async) {
	stk::ReadCounter counter(data, par);
	sbam::ReadInfo* read = nullptr, * pair = nullptr;
	AlignExtend extender(&par->seqp);
	auto current = range.begin;
	try {
		bam->setVOffset(current);
		while ((current = bam->voffset()) < range.end &&
			(read = bam->next()) &&
			status->state == stk::RUNNING) {
			// Ignore secondary alignment
			if (read->flag & sbam::SECONDARY_ALIGN_READ || read->flag & sbam::SUPPLEMENTAL) continue;
			// Increment total read count
			if (async) {
				SAutoLock al(par->slock);
				++(data->summary.total);
			}
			else ++(data->summary.total);
			// Discard unmapped read
			if (read->flag & sbam::UNMAPPED_READ) continue;
			// PCR duplication check
			if (par->ignore_dp && (read->flag & sbam::PCR_DUPLICATE)) continue;
			// Extend
			_extendRead(read, &par->reference[read->ref.idx], &extender);
			// Check paired read is mapped
			if (read->flag & sbam::NEXT_UNMAPPED_READ) {
				// Solo read
				// Check clipped read
				if (par->detect_sv) {
					// Add to realign list
					if (stk::headClip(&read->cigars, par->min_clip))
						svd->addClipRead(read, DIRECTION::HEAD);
					if (stk::tailClip(&read->cigars, par->min_clip))
						svd->addClipRead(read, DIRECTION::TAIL);
				}
				// Update read count and total length
				++(data->summary.count[read->ref.idx]);
				// Update total aligned length
				data->summary.bases[read->ref.idx] += read->cigars.refSize();
				// Count depth
				if (read->cigars.size()) counter.count(read);
			}
			else {
				// Search pair
				if (async) {
					SAutoLock al(par->mlock[read->next.idx]);
					pair = pairs->at(read->next.idx).hasKey(read->name) ?
						(sbam::ReadInfo*)pairs->at(read->next.idx)[read->name] : nullptr;
				}
				else pair = pairs->at(read->next.idx).hasKey(read->name) ?
					(sbam::ReadInfo*)pairs->at(read->next.idx)[read->name] : nullptr;
				// Pairing
				if (pair) {
					if (read->ref.dir == pair->ref.dir) {
						// Inverted read pair
						if (par->detect_sv) svd->addClipReadPair(read, pair, true);
					}
					else {
						// Check overlap
						_correctOverlap(read, pair);
						if (par->detect_sv) {
							if (read->cigars.size() && pair->cigars.size()) svd->addClipReadPair(read, pair);
							else {
								if (read->cigars.size()) {
									if (stk::headClip(&read->cigars, par->min_clip))
										svd->addClipRead(read, DIRECTION::HEAD);
									if (stk::tailClip(&read->cigars, par->min_clip))
										svd->addClipRead(read, DIRECTION::TAIL);
								}
								if (pair->cigars.size()) {
									if (stk::headClip(&pair->cigars, par->min_clip))
										svd->addClipRead(pair, DIRECTION::HEAD);
									if (stk::tailClip(&pair->cigars, par->min_clip))
										svd->addClipRead(pair, DIRECTION::TAIL);
								}
							}
						}
					}
					// Update read count
					++(data->summary.count[read->ref.idx]);
					// Update total aligned length
					if (read->cigars.size() && read->cigars.refSize()) {
						// Update total aligned length
						data->summary.bases[read->ref.idx] += read->cigars.refSize();
						// Count depth
						counter.count(read);
					}
					if (pair->cigars.size() && pair->cigars.refSize()) {
						// Update total aligned length
						data->summary.bases[pair->ref.idx] += pair->cigars.refSize();
						// Count depth
						counter.count(pair);
					}
					// Remove from buffer
					pairs->at(read->next.idx).remove(read->name);
				}
				else {
					if (async) {
						SAutoLock al(par->mlock[read->ref.idx]);
						pairs->at(read->ref.idx).set(read->name, new sbam::ReadInfo(*read));
					}
					else pairs->at(read->ref.idx).set(read->name, new sbam::ReadInfo(*read));
				}
			}
			// Update progress
			status->current_task[async ? refidx : 0] = current.file_offset - range.begin.file_offset;
		}
		if (par->detect_sv && status->state == stk::RUNNING) svd->realign();
	}
	catch (Exception ex) { SPrint("Error read.", NL, read->toString(), NL, pair->toString()); stk::interruptProc(status, logger, &ex); }
}

void stk::BamReader::summary2(BamFile* bam, NGSData* data, Array<Map<String, SPointer<sbam::ReadInfo>>> *pairs, SRAnalysis* svd, vchunk range) {
	stk::ReadCounter counter(data, par);
	sbam::ReadInfo* read = nullptr, * pair = nullptr;
	AlignExtend extender(&par->seqp);
	auto current = range.begin;
	try {
		bam->setVOffset(current);
		while ((current = bam->voffset()) < range.end &&
			(read = bam->next()) &&
			status->state == stk::RUNNING) {
			// Ignore secondary alignment
			if (read->flag & sbam::SECONDARY_ALIGN_READ || read->flag & sbam::SUPPLEMENTAL) continue;
			// Increment total read count
			++(data->summary.total);
			// Discard unmapped read
			if (read->flag & sbam::UNMAPPED_READ) continue;
			// PCR duplication check
			if (par->ignore_dp && (read->flag & sbam::PCR_DUPLICATE)) continue;
			// Extend
			_extendRead(read, &par->reference[read->ref.idx], &extender);
			// Check paired read is mapped
			if (read->flag & sbam::NEXT_UNMAPPED_READ) {
				// Solo read
				// Check clipped read
				if (par->detect_sv) {
					// Add to realign list
					if (stk::headClip(&read->cigars, par->min_clip))
						svd->addClipRead(read, DIRECTION::HEAD);
					if (stk::tailClip(&read->cigars, par->min_clip))
						svd->addClipRead(read, DIRECTION::TAIL);
				}
				// Update read count and total length
				++(data->summary.count[read->ref.idx]);
				data->summary.bases[read->ref.idx] += read->seq.size();
				// Update total aligned length
				data->summary.bases[read->ref.idx] += read->cigars.refSize();
				// Count depth
				if (read->cigars.size()) counter.count(read);
			}
			else {
				// Search pair
				pair = pairs->at(read->next.idx).hasKey(read->name) ?
					(sbio::sbam::ReadInfo *)pairs->at(read->next.idx)[read->name] : nullptr;
				// Pairing
				if (pair) {
					if (read->ref.dir == pair->ref.dir) {
						// Inverted read pair
						if (par->detect_sv) svd->addClipReadPair(read, pair, true);
					}
					else {
						// Check overlap
						_correctOverlap(read, pair);
						if (par->detect_sv) {
							if (read->cigars.size() && pair->cigars.size()) svd->addClipReadPair(read, pair);
							else {
								if (read->cigars.size()) {
									if (stk::headClip(&read->cigars, par->min_clip))
										svd->addClipRead(read, DIRECTION::HEAD);
									if (stk::tailClip(&read->cigars, par->min_clip))
										svd->addClipRead(read, DIRECTION::TAIL);
								}
								if (pair->cigars.size()) {
									if (stk::headClip(&pair->cigars, par->min_clip))
										svd->addClipRead(pair, DIRECTION::HEAD);
									if (stk::tailClip(&pair->cigars, par->min_clip))
										svd->addClipRead(pair, DIRECTION::TAIL);
								}
							}
						}
					}
					// Update read count
					++(data->summary.count[read->ref.idx]);
					// Update total aligned length
					if (read->cigars.size() && read->cigars.refSize()) {
						// Update total aligned length
						data->summary.bases[read->ref.idx] += read->cigars.refSize();
						// Count depth
						counter.count(read);
					}
					if (pair->cigars.size() && pair->cigars.refSize()) {
						// Update total aligned length
						data->summary.bases[pair->ref.idx] += pair->cigars.refSize();
						// Count depth
						counter.count(pair);
					}
					// Remove from buffer
					pairs->at(read->next.idx).remove(read->name);
				}
				else pairs->at(read->ref.idx).set(read->name, SPointer<sbam::ReadInfo>(*read));
			}
			// Update progress
				status->current_task[0] = current.file_offset - range.begin.file_offset;
		}
		if (par->detect_sv && status->state == stk::RUNNING) svd->realign();
	} catch (Exception ex) { SPrint("Error read.", NL, read->toString(), NL, pair->toString()); stk::interruptProc(status, logger, &ex); }
}
inline void _runSummary2(stk::BamReader* br, int i, BamFile* bam, NGSData* summary, Array<Map<String, Pointer<sbam::ReadInfo>>>* pairs, stk::SRAnalysis* svd, vchunk range) {
	br->summary2_(i, bam, summary, pairs, svd, range, true);
}

void stk::BamReader::summarize(BamFile *bam, NGSData* data) {
	// Status => Running
	status->setState(stk::RUNNING);
	status->setTask(bam->filesize(), par->reference.size());
	// Log >> check all reads
	logger->log("Started to check aligned reads.");
	// Display progress
	SWrite(SP * 4, "> Progress: ", SP * 4);
	std::thread thread(stk::showProgress, status);
	// Read analysis
	// Multi-thread loading BAM
	if (par->async_load && bam->hasIndex()) {
		// Prepare objects for multi thread proc.
		Array<BamFile> bams(par->reference.size());
		sfor(bams) $_.open(bam->path());
		
		Array<Map<String, Pointer<sbam::ReadInfo>>> pairs(par->reference.size());
		
		Array<SRAnalysis> svds(par->reference.size());
		sfor(svds) { $_.setParam(par); $_.setData(data); }
		// 
		sfori(par->reference) {
			// Range
			vchunk chunk(bam->index.loffsets[i][0],
				(i == par->reference.size() - 1 ? sbam::VOffset(bam->filesize(), 0) : bam->index.loffsets[i + 1][0]));
			if (par->seqtype == sngs::SEQ_TYPE::SINGLE) 
				threads->addTask(_runSummary1, this, i, &bams[i], data, &svds[i], chunk);
			else if (par->seqtype == sngs::SEQ_TYPE::PAIRED) 
				threads->addTask(_runSummary2, this, i, &bams[i], data, &pairs, &svds[i], chunk);
		}
		threads->complete();
		if (par->seqtype == sngs::SEQ_TYPE::PAIRED) {
			// Check reamined unpaired reads
			
			sforin(i, 0, par->reference.size()) {
				SPrint("Unpaired");
				SPrint("  ", i, ":", pairs[i].size());
			}

		}
	}
	// Single-thread loading BAM
	else {
		SRAnalysis svd(data, par);
		//
		// Run
		vchunk chunk(bam->voffset(), sbam::VOffset(bam->filesize(), 0));
		if (par->seqtype == sngs::SEQ_TYPE::SINGLE)
			summary1(0, bam, data, &svd, chunk, false);
		else if (par->seqtype == sngs::SEQ_TYPE::PAIRED) {
			Array<Map<String, SPointer<sbam::ReadInfo>>> pairs(par->reference.size());
			summary2(bam, data, &pairs, &svd, chunk);
		}			
	}
	stk::closeThread(&thread, status);
	totalize(data);
}
// 
void stk::BamReader::totalize(NGSData* data) {
	size_t total = 0;
	sfor(data->variants) total += $_.size();
	status->setState(stk::RUNNING);
	status->setTask(data->variants.size(), 1);
	logger->log("Finalization.");
	try {
		svecu uncov(data->summary.refnum, 0);
		sforin(r, 0, data->summary.refnum) {
			// Correction of depth at the last bin in each LGs
			data->depth[r][-1] *= (float)data->summary.bin / (data->summary.reflen[r] - (data->depth[r].size() - 1) * data->summary.bin);
			// 
			data->summary.avelen += data->summary.bases[r];
			data->summary.bases[r] /= data->summary.count[r];
			// 
			stk::countLowCover(r, data, &uncov[r]);
			//
			sforin(v, 0, SV_TYPE_COUNT) {
				threads->addTask(stk::checkVariants, &data->variants[SV_TYPE_COUNT * r + v], par);
			}
		}
		//
		size_t tlen = par->reference.total(), len = 0;
		if (par->target.size()) {
			sforeach(target, par->target) len += target.length(true);
		}
		else len = tlen;
		// Average depth
		data->summary.avedp = data->summary.avelen / len;
		// Average length of reads
		data->summary.avelen /= sstat::sum(data->summary.count);
		threads->complete();
		data->summary.cover = 1.0 - ((double)sstat::sum(uncov) / tlen);
	}
	catch (Exception ex) { par->logger.log(ex); }
	status->setState(stk::FINISHED);
}