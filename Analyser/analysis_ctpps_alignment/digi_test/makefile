all: distributions

distributions : distributions.cc input_files.h
	g++ -O3 -g -Wall -Wextra -Wno-attributes --std=c++11\
			`root-config --libs` `root-config --cflags` \
			-I"$(CMSSW_RELEASE_BASE)/src" \
			-L"$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc530" \
			-lDataFormatsFWLite \
			-lDataFormatsCommon \
			distributions.cc -o distributions
