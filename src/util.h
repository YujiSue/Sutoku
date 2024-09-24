#ifndef STK_UTIL_H
#define STK_UTIL_H
#include "sutil/sthread.h"
#include "smath/vector.h"
#include "sbioinfo/align.h"
#include "sbioinfo/bam.h"
#include "sapp/sapp.h"
using namespace slib;
using namespace slib::sbio;
using namespace slib::smath;
using namespace slib::sapp;
namespace stk {
	constexpr subyte INITIALIZE = 0x01;
	constexpr subyte RUNNING = 0x02;
	constexpr subyte ERRORED = 0x03;
	constexpr subyte FINISHED = 0x04;
	class Param;
	class Status {
	private:
		SLock _lock;

	public:
		/*
		#ifdef SUB_PROCESS
			SharedMemory smem;
		#endif
		*/
		int state;
		String message;
		size_t total_task;
		svecs current_task;

	public:
		Status();
		~Status();
		void reset();
		void setState(const int st);
		void setTask(const size_t total, const size_t count);
		void setMessage(const char* s);
		double progress();
		//void share();
	};

	extern void showProgress(Status* status);
	extern void closeThread(std::thread* thread, Status* status);
	extern void interruptProc(stk::Status* status, SLogger* logger, Exception* ex);
	extern bool headClip(CigarArray* cigars, size_t sz);
	extern bool tailClip(CigarArray* cigars, size_t sz);
	extern double pairQual(sbam::ReadInfo* read, sbam::ReadInfo* pair, stk::Param* par);
	extern void countLowCover(int r, NGSData* data, suinteger* count);
	extern void checkVariants(Array<SVar>* variants, stk::Param* par);
	extern void integrate(NGSData* data1, NGSData* data2, stk::Param* par);
	extern void integrate(Array<SVar>* variants1, Array<SVar>* variants2, stk::Param* par);
	extern void subtract(Array<SVar> *variants1, Array<SVar>* variants2, stk::Param* par);
	extern void annotate(VarList& variants, stk::Param* par);
}

#define ASSERT_STATE if(status->state!=stk::FINISHED){return;}

#endif