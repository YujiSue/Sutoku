#include "profile.h"
#include "analyzer.h"
using namespace slib;
using namespace slib::sio;
using namespace slib::sutil;
using namespace slib::sapp;

class Sutoku : public SCuiApp {
public:
	Sutoku() : SCuiApp(app_profile, prof_format) {}
	~Sutoku() {}
	int exec() {
		stk::Analyzer analyzer(preference);
		auto res = analyzer.analyze();
		return res.code;
	}
};
RUN_CUI_APP(Sutoku)
