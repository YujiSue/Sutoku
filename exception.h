#ifndef STK_EXCEPTION_H
#define STK_EXCEPTION_H

//#include "sbasic/exception.h"
//using namespace slib;
/*
constexpr suint STK_NO_ERROR = 0x0000;

constexpr suint STK_NO_READ = 0x0001;

constexpr suint STK_NO_BAM_INDEX = 0x0010;

constexpr suint STK_COMPLETED = 0xFFFF;
*/

#include "sapp/sappex.h"
using namespace slib;
using namespace slib::sapp;
namespace stk {
	class AnalyzerException : public SAppException {};
	/*
#define NO_INPUT_ERR 0x0001

#define NO_READ_ERR 0x0011

	class AnalyzerException : public SException {
	public:
		AnalyzerException(const char* f, sint l, const char* func, sint e = 0, const char* target = nullptr, const char* note = nullptr);
		~AnalyzerException();
	};

	struct Status {
		suint state;
		String msg;
		//SLock lock;
		suinteger total_task;
		//svecs current_task;
#ifdef SUB_PROCESS
		SharedMemory smem;
#endif
		Status();
		~Status();

		void init();
		void setState(suint st);
		void setTask(suinteger total, sint count);
		void setMessage(const char* s);
		double progress();
		void share();
	};
	*/
}

#endif