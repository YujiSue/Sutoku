#include "analyzer.h"

/**
* Check SV
*/ 
bool _checkSV(slib::sbio::SVar& var, stk::Param* par) {
	auto& vp = par->varp.svp;
	// Check target region
	if (par->target.size()) {
		if (!par->target[var.pos[0].idx].overlap(var.pos[0]) &&
			!par->target[var.pos[1].idx].overlap(var.pos[1])) return false;
	}
	// Verify
	return vp.min_read <= var.total() && // Reac count
		vp.min_qual <= slib::sbio::sutil::phredVal(var.qual) && // Read qual.
		var.bias() <= vp.read_bias; // Read bias
}


stk::SRAnalysis::SRAnalysis(Analyzer* an) : par(nullptr), status(nullptr), summary(nullptr), threads(nullptr), logger(nullptr) {
	if (an) {
		par = &an->par;
		status = &an->par.status;
		threads = &an->par.threads;
		logger = &an->par.logger;
	}
}
stk::SRAnalysis::SRAnalysis(NGSData* d, stk::Param* p) : SRAnalysis() { setData(d); setParam(p); }
stk::SRAnalysis::~SRAnalysis() {}
void stk::SRAnalysis::setReader(stk::BamReader* br) {
	par = br->par; 
	status = br->status;
	threads = br->threads;
	logger = br->logger;
	trie.setParam(&par->seqp); 
	search.setParam(&par->seqp);
	if (par->async_load) search.setThreads(nullptr);
	else search.setThreads(threads);
}
/**
* Export split/chimeric reads
*/
void stk::SRAnalysis::splitreads(slib::sbio::NGSData* data, slib::IOStream& stream) {
	// Writeout header
	stream.print("Read type\tAligned1\tAligned2\tUnaligned seq \tRead count (Fwd/Rev)\tQuality");

	// Check background
	if (par->control.isLoaded()) 
		sfor2(data->variants, par->control.variants) stk::subtract(&$_1, &$_2, par);

	// 
	auto it = data->variants.begin();
	sforin(i, 0, data->summary.refnum) {
		/* DEL */
		if (par->detect[0]) { sforeach(var, $_) { if (_checkSV(var, par)) { stream.print(var.toString(&par->reference)); } } }
		$NEXT;
		/* DUP */
		if (par->detect[1]) { sforeach(var, $_) { if (_checkSV(var, par)) { stream.print(var.toString(&par->reference)); } } }
		$NEXT;
		/* INS */
		if (par->detect[2]) { sforeach(var, $_) { if (_checkSV(var, par)) { stream.print(var.toString(&par->reference)); } } }
		$NEXT;
		/* INV */
		if (par->detect[4]) { sforeach(var, $_) { if (_checkSV(var, par)) { stream.print(var.toString(&par->reference)); } } }
		$NEXT;
		/* TRS */
		if (par->detect[6]) { sforeach(var, $_) { if (_checkSV(var, par)) { stream.print(var.toString(&par->reference)); } } }
		$NEXT;
		/* TRS+INV */
		if (par->detect[8]) { sforeach(var, $_) { if (_checkSV(var, par)) { stream.print(var.toString(&par->reference)); } } }
		$NEXT;
	}
}
/**
* 
*/
AlignPair* _selectPair(Array<salign*>* candidates, SVar* var) {
	//
	candidates->sort([](salign* a1, salign* a2) { return a2->score < a1->score; });
	//
	float total_score = 0.f, max_score = candidates->at(0)->score;
	salign* best = candidates->at(0);
	auto clip = var->pos[0].idx == -1 ? DIRECTION::HEAD : DIRECTION::TAIL;
	auto pos = var->pos[0].idx == -1 ? &var->pos[1] : &var->pos[0];
	
	/* DEBUG
	if (candidates->size() > 1 && candidates->at(1)->score == max_score) {
		SPrint((clip == DIRECTION::HEAD ? "*** | " : ""), pos->idx, ":", pos->begin, "-", pos->end, (clip == DIRECTION::TAIL ? " | ***" : ""));
		sforeach(al, *candidates) {
			SPrint(al->ref.idx, ":", al->ref.begin, "-", al->ref.end, ":", al->score);
		}
	}
	*/
	
	
	//
	sfor(*candidates) {
		total_score += $_->score;
		if ($_->score < max_score) continue;
		if ($_->ref.idx == best->ref.idx) {
			if (best->ref.idx == pos->idx) {
				if ($_->ref.dir == best->ref.dir) {
					if (clip == DIRECTION::HEAD) {
						if (best->ref.dir) {
							auto dist1 = std::abs(pos->begin - best->ref.begin);
							auto dist2 = std::abs(pos->begin - $_->ref.begin);
							if (dist2 < dist1) best = $_;
						}
						else {
							auto dist1 = std::abs(pos->begin - best->ref.end);
							auto dist2 = std::abs(pos->begin - $_->ref.end);
							if (dist2 < dist1) best = $_;
						}
					}
					else {
						if (best->ref.dir) {
							auto dist1 = std::abs(pos->end - best->ref.end);
							auto dist2 = std::abs(pos->end - $_->ref.end);
							if (dist2 < dist1) best = $_;
						}
						else {
							auto dist1 = std::abs(pos->end - best->ref.begin);
							auto dist2 = std::abs(pos->end - $_->ref.begin);
							if (dist2 < dist1) best = $_;
						}
					}
				}
				else if (!$_->ref.dir) best = $_;
			}
		}
		else if ($_->ref.idx == pos->idx) best = $_;
	}
	best->score = 1.f - (best->score / total_score);


	/* DEBUG
	if (candidates->size() > 1 && candidates->at(1)->score == max_score) {
		SPrint("[*] ", best->ref.idx, ":", best->ref.begin, "-", best->ref.end, ":", best->score);
	}
	*/


	return best;
}

Pair<AlignPair*, AlignPair*> _selectLinker(Array<salign*>* candidates, SVar* var, stk::Param *par) {
	auto border = var->alt.find("N");
	float total_score[2] = { 0.f, 0.f }, max_score[2] = { 0.f, 0.f };
	salign* best[2] = { nullptr, nullptr };
	//
	srange target;
	bool update = false;
	// Search best aligned for left part
	if (par->min_clip <= 
		(border == NOT_FOUND ? var->alt.size() : border)) {
		target = srange(0, (int)border - 1);
		candidates->sort([target](salign* a1, salign* a2) {
			return (a2->query & target).length() < (a1->query & target).length()
				&& a2->score < a1->score;
			});
		if ((candidates->at(0)->query & target).length()) {
			best[0] = candidates->at(0);
			max_score[0] = best[0]->score;
			sfor(*candidates) {
				if (!($_->query & target).length()) break;
				total_score[0] += $_->score;
				if ($_->score < max_score[0]) continue;
				// Check
				update = false;
				if ($_->ref.idx == var->pos[1].idx &&
					$_->ref.dir == var->pos[1].dir) {
					//
					if (best[0]->ref.idx == var->pos[1].idx &&
						best[0]->ref.dir == var->pos[1].dir) {
						if (var->pos[1].dir) {
							// update best if located downstream of right
							if (var->pos[1].begin < $_->ref.end &&
								best[0]->ref.end < var->pos[1].begin) update = true;
						}
						else {
							// update best if located upstream of right
							if ($_->ref.begin < var->pos[1].end &&
								var->pos[1].end < best[0]->ref.begin) update = true;
						}

					}
					//
					else update = true;
				}
				if (update) best[0] = $_;
			}
		}
	}
	// Search best aligned for right part
	if (par->min_clip <= 
		var->alt.size() - (border == NOT_FOUND ? var->alt.size() : (border + 1))) {
		target = srange((int)border + 1, (int)var->alt.size() - 1);
		candidates->sort([target](salign* a1, salign* a2) {
			return (a2->query & target).length() < (a1->query & target).length()
				&& a2->score < a1->score;
			});
		if ((candidates->at(0)->query & target).length()) {
			best[1] = candidates->at(0);
			max_score[1] = best[1]->score;
			sfor(*candidates) {
				if (!($_->query & target).length()) break;
				total_score[1] += $_->score;
				if ($_->score < max_score[1]) continue;
				// Check
				update = false;
				if ($_->ref.idx == var->pos[0].idx &&
					$_->ref.dir == var->pos[0].dir) {
					//
					if (best[1]->ref.idx == var->pos[0].idx &&
						best[1]->ref.dir == var->pos[0].dir) {
						if (var->pos[0].dir) {
							// update best if located upstream of left
							if ($_->ref.begin < var->pos[0].end &&
								var->pos[0].end < best[1]->ref.begin) update = true;
						}
						else {
							// update best if located downstream of left
							if (var->pos[0].begin < $_->ref.end &&
								best[1]->ref.end < var->pos[0].begin) update = true;
						}
					}
					//
					else update = true;
				}
				if (update) best[1] = $_;
			}
		}
	}
	//
	if (best[0]) {
		best[0]->score = 1.0f - (best[0]->score / total_score[0]);
		if (best[0]->score < par->min_err_prob) best[0]->score = (float)par->min_err_prob;
		if (sbio::sutil::phredVal(best[0]->score) < par->varp.svp.min_qual) best[0] = nullptr;
		else if (border != NOT_FOUND) sbio::sutil::trimAlignTail(*best[0], -1, (int)border);
	}
	//
	if (best[1]) {
		best[1]->score = 1.0f - (best[1]->score / total_score[1]);
		if (sbio::sutil::phredVal(best[1]->score) < par->varp.svp.min_qual) best[1] = nullptr;
		else if (border != NOT_FOUND) sbio::sutil::trimAlignHead(*best[1], -1, (int)border);
	}
	return Pair<AlignPair*, AlignPair*>(best[0], best[1]);
}
void _makeVariant1(AlignPair *al, SVar* var, ubytearray *que, NGSData *data, stk::Param* par) {
	// Head clip
	if (var->pos[0].idx == -1) {
		// Short insertion
		if (al->query.end + 1 < que->size()) {
			auto sz = que->size() - al->query.end - 1;
			var->alt.resize(sz);
			sdna::decode(que->data(), al->query.end + 1, sz, (subyte*)&var->alt[0]);
		}
		// set pos
		var->pos[0] = al->ref;
	}
	// Tail clip
	else if (var->pos[1].idx == -1) {
		// Short insertion
		if (0 < al->query.begin) {
			var->alt.resize(al->query.begin);
			sdna::decode(que->data(), 0, al->query.begin, (subyte*)&var->alt[0]);
		}
		// set pos
		var->pos[1] = al->ref;
	}
	// 
	
	//if (0 < var->read[1]) sna::toComplement(var->alt);
	
	// Adjust quality value
	var->qual = 1.0 - ((1.0 - al->score) * (1.0 - var->qual));
	if (var->qual > 1.0 - par->min_err_prob) var->qual = 1.0 - par->min_err_prob;
	if (var->qual < par->min_err_prob) var->qual = par->min_err_prob;
	// Decide variant type
	var->categorize();
	if ((var->type & 0xFF) != 0) data->addVariant(*var);
}

void _makeVariant2(Pair<AlignPair*, AlignPair*> &als, SVar* var, NGSData *data, stk::Param *par) {
	// Link
	bool link = false;
	if (als.first) {
		//
		if (als.first->ref.idx == var->pos[1].idx &&
			als.first->ref.dir == var->pos[1].dir) {
			//
			if (als.first->ref.dir) {
				//
				if (var->pos[1].begin < als.first->ref.end &&
					als.first->ref.begin <= var->pos[1].end + par->varp.svp.max_gap) {
					//
					if (als.second) als.second = nullptr;
					if (als.first->ref.end <= var->pos[1].end) als.first->cigars.clear();
					else if (als.first->ref.overlap(var->pos[1])) {
						sbio::sutil::trimHeadTo(als.first->ref, als.first->query, als.first->cigars, var->pos[1].end, true);
					}
					var->pos[1].end = als.first->ref.end;
					var->qual = 1.0 - ((1.0 - als.first->score) * (1.0 - var->qual));
					if (var->qual < par->min_err_prob) var->qual = par->min_err_prob;
					var->alt.resize(als.first->query.begin);
					link = true;
				}
			}
			else {
				if (als.first->ref.begin < var->pos[1].end &&
					var->pos[1].begin <= als.first->ref.end + par->varp.svp.max_gap) {
					//
					if (als.second) als.second = nullptr;
					if (var->pos[1].begin <= als.first->ref.begin) als.first->cigars.clear();
					else if (als.first->ref.overlap(var->pos[1])) {
						sbio::sutil::trimTailTo(als.first->ref, als.first->query, als.first->cigars, var->pos[1].begin, false);
					}
					var->pos[1].begin = als.first->ref.begin;
					var->qual = 1.0 - ((1.0 - als.first->score) * (1.0 - var->qual));
					if (var->qual < par->min_err_prob) var->qual = par->min_err_prob;
					var->alt.resize(als.first->query.begin);
					link = true;
				}
			}
		}
	}
	if (als.second) {
		//
		if (als.second->ref.idx == var->pos[0].idx &&
			als.second->ref.dir == var->pos[0].dir) {
			//
			if (als.second->ref.dir) {
				//
				if (als.second->ref.begin < var->pos[0].end &&
					var->pos[0].begin <= als.second->ref.end + par->varp.svp.max_gap) {
					//
					if (als.first) als.first = nullptr;
					if (var->pos[0].begin <= als.second->ref.begin) als.second->cigars.clear();
					else if (als.second->ref.overlap(var->pos[0])) {
						sbio::sutil::trimTailTo(als.second->ref, als.second->query, als.second->cigars, var->pos[0].begin, true);
					}
					var->pos[0].begin = als.second->ref.begin;
					var->qual = 1.0 - ((1.0 - als.second->score) * (1.0 - var->qual));
					if (var->qual < par->min_err_prob) var->qual = par->min_err_prob;
					if (als.second->query.end < var->alt.size() - 1) var->alt.clip(als.second->query.end + 1);
					else var->alt.clear();
					link = true;
				}
			}
			else {
				if (var->pos[0].begin < als.second->ref.end &&
					als.second->ref.begin <= var->pos[0].end + par->varp.svp.max_gap) {
					//
					if (als.first) als.first = nullptr;
					if (als.second->ref.end <= var->pos[0].end) als.second->cigars.clear();
					else if (als.second->ref.overlap(var->pos[0])) {
						sbio::sutil::trimHeadTo(als.second->ref, als.second->query, als.second->cigars, var->pos[0].end, false);
					}
					var->pos[0].end = als.second->ref.end;
					var->qual = 1.0 - ((1.0 - als.second->score) * (1.0 - var->qual));
					if (var->qual < par->min_err_prob) var->qual = par->min_err_prob;
					if (als.second->query.end < var->alt.size() - 1) var->alt.clip(als.second->query.end + 1);
					else var->alt.clear();
					link = true;
				}
			}
		}
	}
	if (link) {
		var->categorize();
		if ((var->type & 0xFF) != 0) data->addVariant(*var);
		if (als.first && als.first->cigars.empty()) als.first = nullptr;
		if (als.second && als.second->cigars.empty()) als.second = nullptr;
		return;
	}
	// Ins
	SVar var2;
	// Check overlap
	if (als.first && als.second) {
		if (als.first->ref.idx == als.second->ref.idx &&
			als.first->ref.dir == als.second->ref.dir) {
			//
			if (als.first->ref.dir) {
				if (als.second->ref.begin < als.first->ref.end &&
					als.first->ref.begin <= als.second->ref.end + par->varp.svp.max_gap) {
					//
					if (als.first->ref.begin <= als.second->ref.begin) als.second->cigars.clear();
					else if (als.first->ref.end <= als.second->ref.end) als.first->cigars.clear();
					else if (als.first->ref.overlap(als.second->ref)) {
						sbio::sutil::trimHeadTo(als.first->ref, als.first->query, als.first->cigars, als.second->ref.end, true);
					}
				}
			}
			//
			else {
				if (als.first->ref.begin < als.second->ref.end &&
					als.second->ref.begin <= als.first->ref.end + par->varp.svp.max_gap) {
					//
					if (als.second->ref.end <= als.first->ref.end) als.second->cigars.clear();
					else if (als.second->ref.begin <= als.first->ref.begin) als.first->cigars.clear();
					else if (als.first->ref.overlap(als.second->ref)) {
						sbio::sutil::trimTailTo(als.first->ref, als.first->query, als.first->cigars, als.second->ref.begin, false);
					}
				}
			}
		}
	}
	// 1st SV
	if (als.first) {
		var2 = *var;
		var2.type = 0;
		var2.pos[1] = als.first->ref;
		var2.qual = 1.0 - ((1.0 - als.first->score) * (1.0 - var->qual));
		if (var2.qual < par->min_err_prob) var2.qual = par->min_err_prob;
		var2.alt.resize(als.first->query.begin);
		var2.categorize();
		if ((var2.type & 0xFF) != 0) data->addVariant(var2);
	}
	// 2nd SV
	if (als.second) {
		var2 = *var;
		var2.pos[0] = als.second->ref;
		var2.qual = 1.0 - ((1.0 - als.second->score) * (1.0 - var->qual));
		if (var2.qual < par->min_err_prob) var2.qual = par->min_err_prob;
		if (als.second->query.end < var->alt.size() - 1) var2.alt.clip(als.second->query.end + 1);
		else var2.alt.clear();
		var2.categorize();
		if ((var2.type & 0xFF) != 0) data->addVariant(var2);
	}
	// 
	if (!als.first && !als.second) {
		if ((var->type & 0xFF) == 0) var->type |= INSERTION;
		var->alt.replace("N", "|");
		data->addVariant(*var);
	}
	if (als.first && als.first->cigars.empty()) als.first = nullptr;
	if (als.second && als.second->cigars.empty()) als.second = nullptr;
}
/**
* Single read clip
*/
void stk::SRAnalysis::addClipRead(sbam::ReadInfo* read, DIRECTION clip) {
	auto sv = &svars.next();
	if (clip == DIRECTION::HEAD) {
		sv->pos[1] = read->ref; sv->pos[1].dir = false;
		trie.addQuery(read->raw().substring(0, read->cigars[0].length));
	}	
	else {
		sv->pos[0] = read->ref; sv->pos[0].dir = false;
		trie.addQuery(read->raw().substring(read->seq.size() - read->cigars[-1].length));
	}
	sv->read[read->ref.dir ? 1 : 0] = 1;
	sv->qual = sbio::sutil::rawVal(read->mapq);
	if (sv->qual > 1.0 - par->min_err_prob) sv->qual = 1.0 - par->min_err_prob;
	if (sv->qual < par->min_err_prob) sv->qual = par->min_err_prob;
	if (par->clip_buffer <= svars.size()) realign();
}
/**
* Paired-end read clip
*/
void stk::SRAnalysis::addClipReadPair(sbam::ReadInfo* read, sbam::ReadInfo* pair, bool inv) {
	SVar* sv = nullptr;
	_buffer.clear();
	sbam::ReadInfo* r1, *r2;
	// 
	if (inv) {
		// Inverted read pair
		r1 = read->ref.idx == pair->ref.idx ?
			(read->ref.begin < pair->ref.begin ? read : pair) :
			(read->ref.idx < pair->ref.idx ? read : pair);
		r2 = read->ref.idx == pair->ref.idx ?
			(read->ref.begin < pair->ref.begin ? pair : read) :
			(read->ref.idx < pair->ref.idx ? pair : read);

		// -/- strand pair
		// Add to summary
		if (r1->ref.dir) {
			// clip
			if (stk::tailClip(&r1->cigars, par->min_clip)) addClipRead(r1, DIRECTION::TAIL);
			if (stk::tailClip(&r2->cigars, par->min_clip)) addClipRead(r2, DIRECTION::TAIL);
			// set type
			sv = &svars.next();
			sv->type = INVERSION;
			if (r1->ref.idx != r2->ref.idx) sv->type |= TRANSLOCATION;
			// set pos
			sv->pos[0] = r1->ref; sv->pos[0].dir = true;
			sv->pos[1] = r2->ref; sv->pos[1].dir = false;
			// inter clip
			if (r1->cigars[0].option == scigar::SCLIP) {
				if (0 < r1->cigars[0].length) {
					sv->alt << sna::complement(r1->raw().substring(0, r1->cigars[0].length)) << "N";
				}
			}
			else sv->type |= (UNCLEAR_TAIL << 8);
			if (r2->cigars[0].option == scigar::SCLIP) {
				if (0 < r2->cigars[0].length) {
					sv->alt << "N" << r2->raw().substring(0, r2->cigars[0].length);
				}
			}
			else sv->type |= (UNCLEAR_HEAD << 8);
			sv->alt.replace("NN", "N");
			//
			if (par->min_clip <= r1->cigars[0].length ||
				par->min_clip <= r2->cigars[0].length) _buffer = sv->alt;
		}
		// +/+ strand pair
		else {
			// clip
			if (stk::headClip(&r1->cigars, par->min_clip)) addClipRead(r1, DIRECTION::HEAD);
			if (stk::headClip(&r2->cigars, par->min_clip)) addClipRead(r2, DIRECTION::HEAD);
			// set type
			sv = &svars.next();
			sv->type = INVERSION;
			if (r1->ref.idx != r2->ref.idx) sv->type |= TRANSLOCATION;
			// set pos
			sv->pos[0] = r1->ref; sv->pos[0].dir = false;
			sv->pos[1] = r2->ref; sv->pos[1].dir = true;
			// inter clip
			if (r1->cigars[-1].option == scigar::SCLIP) {
				if (0 < r1->cigars[-1].length) {
					sv->alt << r1->raw().substring(r1->seq.size() - r1->cigars[-1].length) << "N";
				}
			}
			else sv->type |= (UNCLEAR_TAIL << 8);
			if (r2->cigars[-1].option == scigar::SCLIP) {
				if (0 < r2->cigars[-1].length) {
					sv->alt << "N" << sna::complement(r2->raw().substring(r2->seq.size() - r2->cigars[-1].length));
				}
			}
			else sv->type |= (UNCLEAR_HEAD << 8);
			sv->alt.replace("NN", "N");
			//
			if (par->min_clip <= r1->cigars[-1].length ||
				par->min_clip <= r2->cigars[-1].length) _buffer = sv->alt;
		}
	}
	else {
		// Set +/- strand pair
		r1 = read->ref.dir ? pair : read;
		r2 = read->ref.dir ? read : pair;
		// clip
		if (stk::headClip(&r1->cigars, par->min_clip)) addClipRead(r1, DIRECTION::HEAD);
		if (stk::tailClip(&r2->cigars, par->min_clip)) addClipRead(r2, DIRECTION::TAIL);
		// set pos
		sv = &svars.next();
		sv->pos[0] = r1->ref; sv->pos[0].dir = false;
		sv->pos[1] = r2->ref; sv->pos[1].dir = false;
		// inter clip
		if (r1->cigars[-1].option == scigar::SCLIP) {
			if (0 < r1->cigars[-1].length) {
				sv->alt << r1->raw().substring(r1->seq.size() - r1->cigars[-1].length, r1->cigars[-1].length) << "N";
			}
		}
		else sv->type |= (UNCLEAR_TAIL << 8);
		if (r2->cigars[0].option == scigar::SCLIP) {
			if (0 < r2->cigars[0].length) {
				sv->alt << "N" << r2->raw().substring(0, r2->cigars[0].length);
			}
		}
		else sv->type |= (UNCLEAR_HEAD << 8);
		sv->alt.replace("NN", "N");
		// set type
		if (sv->isStrict()) {
			sv->categorize();
		}
		else {
			if (r1->ref.idx == r2->ref.idx) {
				// Tandem duplication
				if (r2->ref.end < r1->ref.begin) sv->type |= DUPLICATION;
				// Deletion
				else if (r1->ref.end + par->varp.svp.max_gap < r2->ref.begin) sv->type |= DELETION;
				// Insertion
				else if (sv->alt.size()) sv->type |= INSERTION;
			}
			else sv->type |= TRANSLOCATION;
		}
		if (stk::tailClip(&r1->cigars, par->min_clip) ||
			stk::headClip(&r2->cigars, par->min_clip)) {
			_buffer = sv->alt;
		}
	}
	// common
	if (_buffer.empty() && (sv->type&0xFF) == 0) {
		svars.resize(svars.size() - 1);
		return;
	}
	if (r1->flag & sbam::FIRST_SEGMENT) sv->read[0] = 1;
	else sv->read[1] = 1;
	sv->qual = stk::pairQual(read, pair, par);
	if (sv->qual < par->min_err_prob) sv->qual = par->min_err_prob;

	// add to realign query (* process all candidates even _buffer is empty)
	trie.addQuery(_buffer);
	if (par->clip_buffer <= svars.size()) realign();
}
void stk::SRAnalysis::realign() {
	// Check candidates
	if (svars.empty()) return;
	// Init.
	stk::ReadCounter counter(summary, par);
	Array<AlignPair*> candidates;
	// Local search
	trie.complete();
	search.search(par->reference, trie);

	//
	auto rnum = par->reference.size();
	auto var = svars.data();
	//
	sforin(q, 0, svars.size()) {
		candidates.clear();
		auto aligns = search.aligns[q];
		sforin(r, 0, rnum) {
			sfor(*aligns) candidates.add($.ptr());
			++aligns;
		}
		//
		if (candidates.empty()) {
			// Paired 
			if (var->pos[0].idx != -1 && var->pos[1].idx != -1) {
				//
				if (var->alt.size()) {
					if ((var->type&0xFF) == 0) var->type |= INSERTION;
					var->alt.replace("N", "|");
					summary->addVariant(*var);
				}
			}
			//
		}
		else {
			// Single
			if (var->pos[0].idx == -1 || var->pos[1].idx == -1) {
				auto best = _selectPair(&candidates, var);
				if (best) {

					/* DEBUG 
					SPrint(best->ref.idx, ":", best->ref.begin, "-", best->ref.end);
					SPrint(best->query.begin, "-", best->query.end, "/", trie.queries[2 * q].size());
					SPrint(best->cigars.toString());

					String rseq = par->reference.raw(best->ref);
					String qseq;
					qseq.resize(trie.queries[2 * q].size());
					sbio::sdna::decode(&trie.queries[2 * q][0], best->query.begin, best->query.length(true), (subyte*)&qseq[0]);
					

					SPrint(best->alref(rseq));
					SPrint(best->match());
					SPrint(best->alque(qseq));
					
					SWrite(var->qual);
					*/

					_makeVariant1(best, var, &trie.queries[2 * q], summary, par);

					/* DEBUG
					SPrint(">>", var->qual);
					*/

					if (par->async_load) par->slock.lock();
					if (best->ref.dir) best->cigars.reverse();
					counter.count(&best->ref, &best->cigars);
					if (par->async_load) par->slock.unlock();
				}
			}
			// Paired
			else {
				auto best = _selectLinker(&candidates, var, par);
				_makeVariant2(best, var, summary, par);
				if (best.first) counter.count(&best.first->ref, &best.first->cigars);
				if (best.second) counter.count(&best.second->ref, &best.second->cigars);
			}
		}
		++var;
	}
	reset();
}

void stk::SRAnalysis::setData(NGSData* d) { summary = d; }
void stk::SRAnalysis::setParam(Param* p) {
	par = p;
	svars.reserve(par->clip_buffer + 1);
	trie.setParam(&par->seqp);
	trie.reserve(par->clip_buffer * 2 + 1);
	search.setParam(&par->seqp);
	// Multi thread search
	if (!par->async_load) search.setThreads(&par->threads);
}
void stk::SRAnalysis::reset() {
	svars.clear();
	trie.reset();
	search.reset();
}