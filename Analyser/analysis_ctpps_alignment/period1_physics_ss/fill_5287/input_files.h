#include <string>
#include <vector>

//----------------------------------------------------------------------------------------------------

std::vector<std::string> input_files;

std::string input_ntuple_name;

void InitInputFiles()
{
	input_ntuple_name = "TotemNtuple";

	input_files.clear();

	std::string prefix;
	char buf[100];
   
	prefix = "root://eostotem.cern.ch//eos/totem/data/ctpps/reduction/version3/280349/";
	for (int idx : { 0, 1, 2, 9, 10, 11, 12, 13, 14, 15 })
	{
		sprintf(buf, "%i", idx);
		input_files.push_back(prefix + buf + "/reco.root");
	}
   
	prefix = "root://eostotem.cern.ch//eos/totem/data/ctpps/reduction/version3/280364/";
	for (int idx : { 0, 1, 2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19 })
	{
		sprintf(buf, "%i", idx);
		input_files.push_back(prefix + buf + "/reco.root");
	}
}
