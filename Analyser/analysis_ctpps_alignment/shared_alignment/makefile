all: collect_alignments

collect_alignments: collect_alignments.cc ../shared_track.h ../shared_alignment.h
	g++ -O3 -g -Wall -Wextra --std=c++11\
		`root-config --libs` `root-config --cflags` \
		-I"$(CMSSW_RELEASE_BASE)/src" \
		-L"$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc530" \
		-lDataFormatsFWLite \
		-lDataFormatsCommon \
		collect_alignments.cc -o collect_alignments
