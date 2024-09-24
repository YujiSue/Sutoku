set(CODES main.cpp analyzer.cpp param.cpp util.cpp bamreader.cpp sranalysis.cpp cnanalysis.cpp varcheck.cpp varsearch.cpp)
add_executable(sutoku 
	main.cpp
	analyzer.cpp
	param.cpp
	util.cpp
	bamreader.cpp
	sranalysis.cpp
	cnanalysis.cpp
	varcheck.cpp
	varsearch.cpp)
install(TARGETS sutoku)