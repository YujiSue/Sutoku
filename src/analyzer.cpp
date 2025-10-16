#include "analyzer.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sbio;
using namespace slib::sbio::sutil;
using namespace stk;
//
typedef int vfilter(sbio::VarList &vl, const stk::Param& par, const SDictionary& pref);
// 
Analyzer::Analyzer() {}
// Constructor
Analyzer::Analyzer(SDictionary &pref) : Analyzer() { setParam(pref); }
// Constructor
Analyzer::~Analyzer() {}

// Export parameter template file (JSON)
Response exportParameterTemplate(stk::Param &par) {
	try {
		Response res;
		if (par.output.empty()) {
			if (par.outdir.empty()) par.output = sfs::joinPath(ssys::current(), "param.json");
			else par.output = sfs::joinPath(par.outdir, "param.json");
		}
		res.output = par.output;
		par.logger.log("Save to '" + res.output + "'");
		par.save(res.output);
		par.logger.log("Completed.");
		return res;
	}
	catch (Exception ex) {
		par.logger.log(ex);
		return Response(ex);
	}
}

// Export read infoformation (SAM format)
Response getReadInfo(stk::Analyzer *an) {
	try {
		// Init.
		Response res;
		auto& par = an->par;
		BamFile bam;
		BamReader reader(an);
		IOStream ostream;
		
		// Proc. each BAM
		sfor(par.inputs) {
			par.logger.log("Open '" + $_ + "'.");
			
			/* Open BAM */
			bam.open($_);

			/* Set output stream */
			SFile f;
			if (par.outdir.empty())
				f.open(sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".read.txt"), MAKE);
			else f.open(sfs::joinPath(an->par.outdir, sfs::fileName($_, false) + ".read.txt"), MAKE);
			ostream.setFileOStream(f);

			/* Log */
			par.logger.log("Export reads to '" + f.path() + "'.");

			/* Writeout read info */
			reader.readinfo(&bam, ostream);

			/* Log */
			par.logger.log("Exported.");
		}
		/* Log */
		par.logger.log("Completed.");
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}
/**
* Get depth info.
*/
Response getDepth(stk::Analyzer* an) {
	try {
		// Init.
		Response res;
		auto& par = an->par;
		NGSData data;
		CNAnalysis cna(an);
		String out, sep;

		// Set separator char.
		if (par.oformat == "auto") par.oformat = "csv";
		if (par.oformat == "tsv") sep = "\t";
		else sep = ",";

		// Make header
		stringarray header;
		/* Non-targeted */
		if (par.target.empty()) { sfor(par.reference) header.add($_.name); }
		/* Targeted */
		else { sfori(par.reference) { sfor(par.target[i]) header.add(par.reference[i].name + ":" + S($_.begin + 1) + "-" + S($_.end + 1)); } }
		// Proc. each BSM
		sfor(par.inputs) {
			/* Log */
			par.logger.log("Open summary data '" + $_ + "'.");
			/* Load data */
			data.load($_);
			data.print(); // For confirmation
			
			/*  */
			cna.setData(&data);
			if (par.outdir.empty())
				out = sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".depth." + par.oformat);
			else out = sfs::joinPath(par.outdir, sfs::fileName($_, false) + ".depth." + par.oformat);
			SFile f(out, slib::sio::MAKE);
			par.logger.log("Export depth to '" + out + "'.");
			
			f << toString(header, sep) << LF; f.flush();
			//
			smath::Vector<svecf> values(header.size());
			cna.depth(values);
			Matrix<float> mat;
			smath::toMat(mat, values, -1.f);
			f << toString(mat, par.oformat);
			par.logger.log("Finished.");
			cna.reset();
		}
		/* Log */
		par.logger.log("Completed.");
		// Return res.
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}
/**
* Get copy num info.
*/
Response getCopy(stk::Analyzer* an) {
	// Init.
	Response res;
	auto& par = an->par;
	NGSData data;
	CNAnalysis cna(an);

	String out, sep;
	if (par.oformat == "auto") par.oformat = "csv";
	if (par.oformat == "tsv") sep = "\t";
	else sep = ",";
	
	stringarray header;
	if (par.target.empty()) {
		sfor(par.reference) header.add($_.name);
	}
	else {
		sfori(par.reference) {
			sfor(par.target[i]) header.add(par.reference[i].name + ":" + S($_.begin + 1) + "-" + S($_.end + 1));
		}
	}

	
	sfor(par.inputs) {
		par.logger.log("Open summary data '" + $_ + "'.");
		data.load($_);
		data.print();
		cna.setData(&data, (par.control.isLoaded() ? &par.control : nullptr));
		if (par.output.empty()) {
			if (par.outdir.empty())
				out = sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".cp.csv");
			else out = sfs::joinPath(par.outdir, sfs::fileName($_, false) + ".cp.csv");
		}
		else out = par.output;
		SFile f(out, slib::sio::MAKE);
		par.logger.log("Export copy to '" + out + "'.");
		smath::Vector<svecf> values(header.size());
		cna.copynum(values);
		f << toString(header, sep) << NL; f.flush();
		

		
		Matrix<float> mat;
		smath::toMat(mat, values);
		f << toString(mat, par.oformat);
		
		par.logger.log("Finished.");
		cna.reset();
	}
	par.logger.log("Completed.");
	return res;
}
/**
* Export split read list
*/ 
Response getSRead(stk::Analyzer* an) {
	try {
		// Init.
		Response res;
		auto& par = an->par;
		NGSData data;
		stk::SRAnalysis sra(an);
		
		// Proc. each BSM
		sfor(par.inputs) {
			/* Log */
			par.logger.log("Open summary data '" + $_ + "'.");

			/* Load BSM data */
			data.load($_);
			data.print(); // Confirmation

			/* Make file to save */
			SFile f;
			if (par.outdir.empty()) f.open(sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".sreads.txt"), MAKE);
			else f.open(sfs::joinPath(par.outdir, sfs::fileName($_, false) + ".sreads.txt"), MAKE);

			/* Log */
			par.logger.log("Export split/chimeric reads to '" + f.path() + "'.");

			/* Export */
			IOStream fs(f, OSTREAM | FILEIO);
			sra.splitreads(&data, fs);

			/* Log */
			par.logger.log("Exported.");
		}

		/* Log */
		par.logger.log("Completed.");

		/* Return res. */
		return res;
	}
	// Error proc.
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}

/**
* Make summary of results. (= Detect split/chimeric reads and count depth)
*/
Response summerize(stk::Analyzer* an) {
	/*
	auto test = "TTGCTTCACCCATCACACGATAACTTTTCAAGTGCATTTTCTATGGAGAACCCAATAATAATCTTTATTGAGAACAACTGTTCGTCGGAAGAACGCCGAACACGCGGCACACGCGTCCACTCCGAAATGAATACCCGTCGAATCTGTACCACAAATGGCGCAAGCAGGAACCTGAATA";
	ubytearray qtest(strlen(test));
	sdna::encode((const subyte*)test, 0, qtest.size(), &qtest[0]);


	AlignPair ap;
	ap.ref.idx = 3;
	ap.ref.begin = 8228663;
	ap.query.begin = 0;
	ap.cigars = "61M";
	ap.ref.end = ap.ref.begin + ap.cigars.refSize() - 1;
	ap.query.end = ap.cigars.queSize() - 1;

	
	sbio::AlignExtend aext(&an->par.seqp);
	aext.extend(&an->par.reference[3], &qtest, &ap);

	auto rs = an->par.reference[3].raw(ap.ref);
	SPrint(ap.cigars.toString());
	SPrint(ap.alref(rs));
	SPrint(ap.match());
	SPrint(ap.alque(S(test).substring(0, ap.query.end + 1)));



	SPrint(an->par.inputs);

	NGSData ex(an->par.inputs[0]+".bsm");

	stk::checkVariants(&ex.variants[18], &an->par);



	*/



	


	// Init.
	Response res;
	auto& par = an->par;
	BamFile bam;
	NGSData data(&par.reference, par.seqtype, par.depth_bin);
	BamReader br(an);
	par.logger.log("Initialization");

	// Status init.
	par.status.setState(stk::INITIALIZE);

	// Log
	par.logger.log(sstr::rfill("Bin size ...", ' ', 20) + SNumber(par.depth_bin).toString());
	par.logger.log(sstr::rfill("PCR duplicate ...", ' ', 20) + (par.ignore_dp ? "Ignore" : "Count"));
	par.logger.log(sstr::rfill("Detect SV ...", ' ', 20) + SNumber(par.detect_sv).toString());
	par.logger.log(sstr::rfill("Thread count ...", ' ', 20) + SNumber(par.max_thread).toString());

	// Proc. each BAM
	sfor(par.inputs) {
		try {
			/* Log */
			par.logger.log("Open '" + $_ + "'.");

			/* Open BAM and summary init. */
			bam.open($_);
			data.setSource(bam);
			if (bam.info.ref_num != par.reference.size()) {
				par.logger.log("Reference mismatch.");
				continue;
			}

			/* Log */
			par.logger.log("Started to make summary data.");
			
			/* Start to make summary */
			br.summarize(&bam, &data);
			
			/* Check finished normaly */
			if (par.status.state == stk::FINISHED) {
				par.logger.log("Save result to '" + $_ + ".bsm'.");
				/** Write out BSM **/
				data.save($_ + ".bsm");
				par.logger.log("Completed.");
			}
			else par.logger.log("Failed.");
			
			/* Reset */
			bam.close();
			an->reset();
		}
		catch (Exception ex) {
			an->par.logger.log(ex);
			return Response(ex);
		}
	}
	return res;
}

Response vsearch(stk::Analyzer* an) {
	Response res;
	auto& par = an->par;
	VarSearch vs(an);
	NGSData data;
	String out;
	par.status.state = stk::INITIALIZE;
	sfor(par.inputs) {
		try {
			par.logger.log("Load '" + $_ + "'.");
			data.load($_);
			if (data.summary.refnum != par.reference.size()) {
				par.logger.log("Reference mismatch.");
				continue;
			}
			// Show summary
			data.print();
			// Detection
			vs.detect(&data);
			//
			if (par.status.state == stk::FINISHED) {
				if (par.oformat == "auto") par.oformat = "tsv";
				if (par.outdir.empty())
					out = sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".variants." + par.oformat);
				else out = sfs::joinPath(par.outdir, sfs::fileName($_, false) + ".variants." + par.oformat);				//
				par.logger.log("Save result to '" + out + "'.");
				vs.variants.save(out);
				par.logger.log("Completed.");
			}
			else par.logger.log("Failed.");
			data.reset();
			vs.reset();
			an->reset();
		}
		catch (Exception ex) {
			an->par.logger.log(ex);
			return Response(ex);
		}
	}
	return res;
}
Response smvs(stk::Analyzer* an) {
	Response res;
	auto& par = an->par;
	BamFile bam;
	NGSData data(&par.reference, par.seqtype, par.depth_bin);
	String out;
	BamReader br(an);
	VarSearch vs(an);
	par.status.state = stk::INITIALIZE;
	sfor(par.inputs) {
		try {
			// Bam file open
			par.logger.log("Open '" + $_ + "'.");
			bam.open($_);
			if (bam.info.ref_num != par.reference.size()) {
				par.logger.log("Reference mismatch.");
				continue;
			}
			// Run proc. to make summary
			par.logger.log("Started to make summary data.");
			//
			br.summarize(&bam, &data);
			// Check finished normaly
			if (par.status.state == stk::FINISHED) {
				par.logger.log("Save result to '" + $_ + ".bsm'.");
				data.save($_ + ".bsm");
				par.logger.log("Completed.");
			}
			else par.logger.log("Failed.");
			// Close
			bam.close();
			// 
			par.logger.log("Started to search variants."); 
			vs.detect(&data);
			if (par.status.state == stk::FINISHED) {
				if (par.oformat == "auto") par.oformat = "txt";
				if (par.outdir.empty())
					out = sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + ".variants." + par.oformat);
				else out = sfs::joinPath(par.outdir, sfs::fileName($_, false) + ".variants." + par.oformat);
				par.logger.log("Save result to '" + out + "'.");
				vs.variants.save(out);
				par.logger.log("Completed.");
			}
			else par.logger.log("Failed.");
			data.reset();
		}
		catch (Exception ex) {
			an->par.logger.log(ex);
			return Response(ex);
		}
	}
	return res;
}
/**
* Integrate BSM files (e.g. for preparation of control data)
*/
Response integrateSummaries(stk::Analyzer* an) {
	try {
		// Init.
		Response res;
		auto& par = an->par;



		if (par.outdir.empty()) par.outdir = sfs::splitPath(par.inputs[0]).first;
		String out = sfs::joinPath(par.outdir, "integrated.bsm");
		
		// 
		if (par.inputs.size() == 1) {
			par.logger.log("Copy summary to '" + out + "'.");
			sfs::copy(par.inputs[0], out, OVERWRITE);
			return res;
		}
		// 
		NGSData ori, tmp;
		par.logger.log("Open summary '" + par.inputs[0] + "'.");
		ori.load(par.inputs[0]);
		ori.print();
		//
		sforin(it, par.inputs.begin()+1, par.inputs.end()) {
			tmp.reset();
			par.logger.log("Open summary '" + $_ + "'.");
			tmp.load($_);
			tmp.print();
			par.logger.log("Integration.");
			stk::integrate(&ori, &tmp, &par);
		}

		// Log.
		par.logger.log("All files have been integrated.");
		ori.print();

		// Writeout integrated summary
		par.logger.log("Save result to '" + out + "'.");
		ori.save(out);
		
		// Return res.
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}
Response subtractSV(stk::Analyzer* an) {
	try {
		Response res;
		auto& par = an->par;
		/*
		NGSData data;
		String out;
		sfor(par.inputs) {
			data.load($_);
			data.subtract(par.control, par.varp);
			if (par.output.empty()) {
				if (par.outdir.empty())
					out = sfs::joinPath(sfs::splitPath($_).first, sfs::fileName($_, false) + "_subtract.bsm");
				else out = sfs::joinPath(par.outdir, sfs::fileName($_, false) + "_subtract.bsm");
			}
			else out = par.output;
			data.save(out);
		}
		*/
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}

//
Response verifyVariants(stk::Analyzer* an) {
	try {
		Response res;
		auto& par = an->par;
		if (par.oformat == "auto") par.oformat = sfs::extension(par.inputs[0]);
		String out;
		VarFilter filter(&par.reference, &par.annotdb, &par.varp);
		sfor(par.inputs) {
			par.logger.log("Open variant list '" + $_ + "'.");
			VarList vl;
			vl.setReference(&par.reference);
			vl.load($_, &par.reference);
			// Annotation
			if (par.annotation) {
				par.logger.log("Annotation started.");
				stk::annotate(vl, &par);
				par.logger.log("Completed.");
			}
			//
			filter.filter(vl);
			//
			sforeach(vf, par.vfilters) {
				//SPrint(vf.toString());
				par.logger.log(S("Run filter '") << vf["plugin"] << "'.");
				sapp::SPlugIn<sbio::VarList&, const stk::Param&, const SDictionary&> plugin(vf["plugin"], "vfilter");
				auto fres = plugin.exec(vl, par, vf["args"]);
			}
			if (par.export_filtered) {
				vl.sort([](const sgenvar& t1, const sgenvar& t2) {
					if (t1->flag == NOT_USE_FLAG || t1->flag == UNAVAILABLE_FLAG) {
						if (t2->flag == NOT_USE_FLAG || t2->flag == UNAVAILABLE_FLAG) return (*t1) < (*t2);
						else return false;
					}
					else if (t2->flag == NOT_USE_FLAG || t2->flag == UNAVAILABLE_FLAG) return true;
					else return (*t1) < (*t2);
					});
			}
			else vl.tidyUp();

			//
			if (par.outdir.empty())
				out = sfs::joinPath(sfs::splitPath(par.inputs[0]).first, sfs::fileName(par.inputs[0], false) + "_verified." + par.oformat);
			else out = sfs::joinPath(par.outdir, sfs::fileName(par.inputs[0], false) + "_verified." + par.oformat);
			par.logger.log("Save verified list to '" + out + "'.");
			//
			vl.save(out, "cols=ID,Chr,Pos,Len,Ref,Alt/Ins,Type,Qual,Freq,Genotype,Cov,Allele Cov,Gene,Site,Mutation,Substitution,Effect,Filter");

		}
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}

// 
Response mergeVariants(stk::Analyzer* an) {
	try {
		Response res;
		auto& par = an->par;
		String out;
		VarFilter filter(&par.reference, &par.annotdb);
		VarList merge, vl;
		merge.setReference(&par.reference);
		sfor(par.inputs) {
			par.logger.log("Open variant list '" + $_ + "'.");
			vl.load($_);
			par.logger.log("Merge variants.");
			filter.merge(merge, vl);
			par.logger.log("Completed.");
			vl.clearAll();
		}
		// Annotation
		if (par.annotation) {
			par.logger.log("Annotation started.");
			stk::annotate(merge, &par);
			par.logger.log("Completed.");
		}
		//
		if (par.outdir.empty())
			out = sfs::joinPath(sfs::splitPath(par.inputs[0]).first, sfs::fileName(par.inputs[0], false) + "_merge." + par.oformat);
		else out = sfs::joinPath(par.outdir, sfs::fileName(par.inputs[0], false) + "_merge." + par.oformat);
		par.logger.log("Save merged list to '" + out + "'.");
		merge.save(out);
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}
//
Response uniqueVariants(stk::Analyzer* an) {
	try {
		Response res;
		auto& par = an->par;
		if (par.oformat == "auto") par.oformat = sfs::extension(par.inputs[0]);
		String out;
		VarFilter filter(&par.reference, &par.annotdb, &par.varp);
		VarList uni;
		sfor(par.inputs) {
			par.logger.log("Open variant list '" + $_ + "'.");
			uni.clearAll();
			uni.setReference(&par.reference);
			uni.load($_, &par.reference);
			filter.subtract(uni, par.vcontrol);
			// Annotation
			if (par.annotation) {
				par.logger.log("Annotation started.");
				stk::annotate(uni, &par);
				par.logger.log("Completed.");
			}
			//
			filter.filter(uni);
			//
			sforeach(vf, par.vfilters) {
				//SPrint(vf.toString());
				par.logger.log(S("Run filter '") << vf["plugin"] << "'.");
				sapp::SPlugIn<sbio::VarList&, const stk::Param&, const SDictionary &> plugin(vf["plugin"], "vfilter");
				auto fres = plugin.exec(uni, par, vf["args"]);
			}
			//
			if (par.export_filtered) {
				uni.sort([](const sgenvar& t1, const sgenvar& t2) {
					if (t1->flag == NOT_USE_FLAG || t1->flag == UNAVAILABLE_FLAG) {
						if (t2->flag == NOT_USE_FLAG || t2->flag == UNAVAILABLE_FLAG) return (*t1) < (*t2);
						else return false;
					}
					else if (t2->flag == NOT_USE_FLAG || t2->flag == UNAVAILABLE_FLAG) return true;
					else return (*t1) < (*t2);
					});
			}
			else uni.tidyUp();

			//
			if (par.outdir.empty())	
				out = sfs::joinPath(sfs::splitPath(par.inputs[0]).first, sfs::fileName(par.inputs[0], false) + "_unique." + par.oformat);
			else out = sfs::joinPath(par.outdir, sfs::fileName(par.inputs[0], false) + "_unique." + par.oformat);
			par.logger.log("Save unique list to '" + out + "'.");
			//
			uni.save(out, "cols=ID,Chr,Pos,Len,Ref,Alt/Ins,Type,Qual,Freq,Genotype,Cov,Allele Cov,Gene,Site,Mutation,Substitution,Effect,Filter");
			
		}
		return res;
	}
	catch (Exception ex) {
		an->par.logger.log(ex);
		return Response(ex);
	}
}
//
Response commonVariants(stk::Analyzer* an) {
	try {
		Response res;
		auto& par = an->par;
		String out;
		VarFilter filter(&par.reference, &par.annotdb, &par.varp);
		VarList com, vl;
		com.setReference(&par.reference);
		com.load(par.inputs[0], &par.reference);
		sforin(vit, par.inputs.begin() + 1, par.inputs.end()) {
			par.logger.log("Open variant list '" + (*vit) + "'.");
			vl.load(*vit, &par.reference);
			filter.common(com, vl);

		}
		// Annotation
		if (par.annotation) {
			par.logger.log("Annotation started.");
			stk::annotate(com, &par);
			par.logger.log("Completed.");
		}
		//
		filter.filter(com);
		//
		sforeach(vf, par.vfilters) {
			//SPrint(vf.toString());
			par.logger.log(S("Run filter '") << vf["plugin"] << "'.");
			sapp::SPlugIn<sbio::VarList&, const stk::Param&, const SDictionary&> plugin(vf["plugin"], "vfilter");
			auto fres = plugin.exec(com, par, vf["args"]);
		}
		//
		if (par.export_filtered) {
			com.sort([](const sgenvar& t1, const sgenvar& t2) {
				if (t1->flag == NOT_USE_FLAG || t1->flag == UNAVAILABLE_FLAG) {
					if (t2->flag == NOT_USE_FLAG || t2->flag == UNAVAILABLE_FLAG) return (*t1) < (*t2);
					else return false;
				}
				else if (t2->flag == NOT_USE_FLAG || t2->flag == UNAVAILABLE_FLAG) return true;
				else return (*t1) < (*t2);
				});
		}
		else com.tidyUp();

		//
		if (par.outdir.empty())
			out = sfs::joinPath(sfs::splitPath(par.inputs[0]).first, sfs::fileName(par.inputs[0], false) + "_common." + par.oformat);
		else out = sfs::joinPath(par.outdir, sfs::fileName(par.inputs[0], false) + "_common." + par.oformat);
		com.save(out);
		return res;
	}
	catch (Exception ex) {
		ex.print();
		return Response(ex);
	}
}

// Hub function
Response Analyzer::analyze() {
	if (par.command == "template") return exportParameterTemplate(par);
	if (par.inputs.empty()) {
		par.logger.log("No input.");
		return Response(sapp::INSUFFICIENT_ARGS_ERROR, "No input.", "Input file(s) is/are required. Please refer help for details.");
	}
	if (par.command == "readinfo") return getReadInfo(this);
	else if (par.command == "depth") return getDepth(this);
	else if (par.command == "copynum") return getCopy(this);
	else if (par.command == "splitread") return getSRead(this);
	else if (par.command == "summary") return summerize(this);
	else if (par.command == "vsearch") return vsearch(this);
	else if (par.command == "analyze") return smvs(this);
	else if (par.command == "integrate") return integrateSummaries(this);
	else if (par.command == "subtract") return subtractSV(this);
	else if (par.command == "verify") return verifyVariants(this);
	else if (par.command == "merge") return mergeVariants(this);
	else if (par.command == "unique") return uniqueVariants(this);
	else if (par.command == "common") return commonVariants(this);
	else {
		par.logger.log("'" + par.command + "' is NOT supported command under this version.");
		return Response();
	}
}
// Setter
void Analyzer::setParam(SDictionary &pref) { par.set(pref); }
// 
void Analyzer::reset() {
	if (par.threads.isWorking()) par.threads.complete();
	par.status.reset();
}
