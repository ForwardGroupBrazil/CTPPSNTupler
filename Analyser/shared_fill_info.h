struct FillInfo
{
  unsigned int fillNumber;
  unsigned int runMin, runMax;
  string alignmentTag;

  FillInfo(unsigned int _fn=0, unsigned int _rmin=0, unsigned int _rmax=0, const string &at="") :
    fillNumber(_fn), runMin(_rmin), runMax(_rmax), alignmentTag(at)
    {
    }
};

//----------------------------------------------------------------------------------------------------

struct FillInfoCollection : public vector<FillInfo> 
{
  FillInfo FindByRun(unsigned int run) const
  {
    for (const auto it : *this)
    {
      if (it.runMin <= run && it.runMax >= run)
	return it;
    }

    printf("ERROR in FillInfoCollection::FindByRun > run %u not found in the collection.\n", run);
    exit(1);

    return FillInfo();
  }
};

//----------------------------------------------------------------------------------------------------

FillInfoCollection fillInfoCollection;

//----------------------------------------------------------------------------------------------------

void InitFillInfoCollection()
{
  fillInfoCollection.push_back(FillInfo(4947, 273725, 273730, "analysis_ctpps_alignment/period1_physics_margin/fill_4947/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4953, 274094, 274094, "analysis_ctpps_alignment/period1_physics_margin/fill_4953/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4961, 274198, 274200, "analysis_ctpps_alignment/period1_physics_margin/fill_4961/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4964, 274240, 274241, "analysis_ctpps_alignment/period1_physics_margin/fill_4964/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4964, 274244, 274244, "analysis_ctpps_alignment/period1_physics/fill_4964/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4976, 274282, 274286, "analysis_ctpps_alignment/period1_physics_margin/fill_4976/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4985, 274387, 274388, "analysis_ctpps_alignment/period1_physics/fill_4985/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4988, 274420, 274422, "analysis_ctpps_alignment/period1_physics/fill_4988/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(4990, 274440, 274443, "analysis_ctpps_alignment/period1_physics/fill_4990/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5005, 274954, 274959, "analysis_ctpps_alignment/period1_physics/fill_5005/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5013, 274966, 274971, "analysis_ctpps_alignment/period1_physics/fill_5013/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5017, 274998, 275001, "analysis_ctpps_alignment/period1_physics/fill_5017/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5020, 275059, 275074, "analysis_ctpps_alignment/period1_physics/fill_5020/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5021, 275124, 275125, "analysis_ctpps_alignment/period1_physics/fill_5021/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5024, 275282, 275293, "analysis_ctpps_alignment/period1_physics/fill_5024/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5026, 275309, 275311, "analysis_ctpps_alignment/period1_physics/fill_5026/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5027, 275319, 275338, "analysis_ctpps_alignment/period1_physics/fill_5027/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5028, 275344, 275345, "analysis_ctpps_alignment/period1_physics/fill_5028/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5029, 275370, 275371, "analysis_ctpps_alignment/period1_physics/fill_5029/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5030, 275375, 275376, "analysis_ctpps_alignment/period1_physics/fill_5030/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5038, 275656, 275659, "analysis_ctpps_alignment/period1_physics/fill_5038/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5043, 275757, 275783, "analysis_ctpps_alignment/period1_physics/fill_5043/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5045, 275828, 275847, "analysis_ctpps_alignment/period1_physics/fill_5045/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5048, 275886, 275890, "analysis_ctpps_alignment/period1_physics/fill_5048/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5052, 275911, 275931, "analysis_ctpps_alignment/period1_physics/fill_5052/2017_01_17"));

  fillInfoCollection.push_back(FillInfo(5261, 279760, 279767, "analysis_ctpps_alignment/period1_physics/fill_5261/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5264, 279794, 279794, "analysis_ctpps_alignment/period1_physics/fill_5264/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5265, 279823, 279823, "analysis_ctpps_alignment/period1_physics/fill_5265/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5266, 279841, 279841, "analysis_ctpps_alignment/period1_physics/fill_5266/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5267, 279844, 279865, "analysis_ctpps_alignment/period1_physics/fill_5267/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5274, 279931, 279931, "analysis_ctpps_alignment/period1_physics/fill_5274/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5275, 279966, 279966, "analysis_ctpps_alignment/period1_physics/fill_5275/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5276, 279975, 279975, "analysis_ctpps_alignment/period1_physics/fill_5276/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5277, 279993, 280024, "analysis_ctpps_alignment/period1_physics/fill_5277/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5279, 280187, 280194, "analysis_ctpps_alignment/period1_physics/fill_5279/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5287, 280327, 280364, "analysis_ctpps_alignment/period1_physics/fill_5287/2017_01_17"));
  fillInfoCollection.push_back(FillInfo(5288, 280383, 280385, "analysis_ctpps_alignment/period1_physics/fill_5288/2017_01_17"));

  fillInfoCollection.push_back(FillInfo(5393, 282730, 282735, "TODO"));
  fillInfoCollection.push_back(FillInfo(5401, 282920, 282924, "TODO"));
  fillInfoCollection.push_back(FillInfo(5405, 283040, 283043, "TODO"));
  fillInfoCollection.push_back(FillInfo(5406, 283049, 283067, "TODO"));
  fillInfoCollection.push_back(FillInfo(5418, 283305, 283308, "TODO"));
  fillInfoCollection.push_back(FillInfo(5421, 283353, 283359, "TODO"));
  fillInfoCollection.push_back(FillInfo(5423, 283407, 283416, "TODO"));
  fillInfoCollection.push_back(FillInfo(5424, 283453, 283453, "TODO"));
  fillInfoCollection.push_back(FillInfo(5427, 283478, 283481, "TODO"));
  fillInfoCollection.push_back(FillInfo(5433, 283548, 283560, "TODO"));
  fillInfoCollection.push_back(FillInfo(5437, 283672, 283685, "TODO"));
  fillInfoCollection.push_back(FillInfo(5439, 283820, 283835, "TODO"));
  fillInfoCollection.push_back(FillInfo(5441, 283863, 283865, "TODO"));
  fillInfoCollection.push_back(FillInfo(5442, 283876, 283878, "TODO"));
  fillInfoCollection.push_back(FillInfo(5443, 283884, 283885, "TODO"));
  fillInfoCollection.push_back(FillInfo(5446, 283933, 283934, "TODO"));
  fillInfoCollection.push_back(FillInfo(5448, 283946, 283964, "TODO"));
  fillInfoCollection.push_back(FillInfo(5450, 284006, 284014, "TODO"));
  fillInfoCollection.push_back(FillInfo(5451, 284025, 284044, "TODO"));
}
