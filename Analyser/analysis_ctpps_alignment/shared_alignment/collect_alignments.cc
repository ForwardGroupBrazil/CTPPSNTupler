#include "../shared_track.h"
#include "../shared_alignment.h"

int main()
{
	// define input
	vector<string> input_dirs = {
		"period1_physics_margin/fill_4947",
		"period1_physics_margin/fill_4953",
		"period1_physics_margin/fill_4961",
		"period1_physics_margin/fill_4964",
		"period1_physics_margin/fill_4976",

		"period1_physics/fill_4964",
		"period1_physics/fill_4985",
		"period1_physics/fill_4988",
		"period1_physics/fill_4990",
		"period1_physics/fill_5005",
		"period1_physics/fill_5013",
		"period1_physics/fill_5017",
		"period1_physics/fill_5020",
		"period1_physics/fill_5021",
		"period1_physics/fill_5024",
		"period1_physics/fill_5026",
		"period1_physics/fill_5027",
		"period1_physics/fill_5028",
		"period1_physics/fill_5029",
		"period1_physics/fill_5030",
		"period1_physics/fill_5038",
		"period1_physics/fill_5043",
		"period1_physics/fill_5045",
		"period1_physics/fill_5048",
		"period1_physics/fill_5052",
		"period1_physics/fill_5261",
		"period1_physics/fill_5264",
		"period1_physics/fill_5265",
		"period1_physics/fill_5266",
		"period1_physics/fill_5267",
		"period1_physics/fill_5274",
		"period1_physics/fill_5275",
		"period1_physics/fill_5276",
		"period1_physics/fill_5277",
		"period1_physics/fill_5279",
		"period1_physics/fill_5287",
		"period1_physics/fill_5288",
	};
	
	// collection
	AlignmentResultsCollection output;

	// process each input directory
	for (const auto dir : input_dirs)
	{
		// get horizontal alignment
		AlignmentResultsCollection input_x;
		string fn_x = "../" + dir + "/process_alignments.out";
		int error_x = input_x.Load(fn_x);
		if (error_x != 0)
		{
			printf("ERROR: can't load alignment results from file '%s'.\n", fn_x.c_str());
			continue;
		}

		auto it_x = input_x.find("method x");
		if (it_x == input_x.end())
		{
			printf("ERROR: can't get result 'method x' from file '%s'.\n", fn_x.c_str());
			continue;
		}

		// get vertical alignment
		AlignmentResultsCollection input_y;
		string fn_y = "../" + dir + "/y_alignment.out";
		int error_y = input_y.Load(fn_y);
		if (error_y != 0)
		{
			printf("ERROR: can't load alignment results from file '%s'.\n", fn_y.c_str());
			continue;
		}

		auto it_y = input_y.find("y alignment");
		if (it_y == input_y.end())
		{
			printf("ERROR: can't get result 'method x' from file '%s'.\n", fn_y.c_str());
			continue;
		}

		// combine horizontal and verctical alignment
		set<unsigned int> rpIds;

		for (const auto p : it_x->second)
			rpIds.insert(p.first);

		for (const auto p : it_y->second)
			rpIds.insert(p.first);

		AlignmentResults ars_combined;

		for (unsigned int rpId : rpIds)
		{
			auto ari_x = it_x->second.find(rpId);
			if (ari_x == it_x->second.end())
			{
				printf("ERROR: can't get x alignment result for RP %i, dir '%s'.\n", rpId, dir.c_str());
				continue;
			}

			auto ari_y = it_y->second.find(rpId);
			if (ari_y == it_y->second.end())
			{
				printf("ERROR: can't get y alignment result for RP %i, dir '%s'.\n", rpId, dir.c_str());
				continue;
			}

			const AlignmentResult &ar_x = ari_x->second;
			const AlignmentResult &ar_y = ari_y->second;

			ars_combined[rpId] = AlignmentResult(ar_x.sh_x, 150E-3, ar_y.sh_y, 100E-3);
		}

		output[dir + "/2017_01_17"] = ars_combined;
	}

	output.Write("collect_alignments.out");

	return 0;
}
